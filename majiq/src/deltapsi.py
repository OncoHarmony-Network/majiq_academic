import h5py
import scipy.misc
import sys
import multiprocessing as mp
import traceback
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
import majiq.src.utils as majiq_utils
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
from majiq.src.io_utils import dump_bin_file, load_bin_file
from majiq.src.psi import deltapsi_posterior, gen_prior_matrix
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper
from voila.api import Voila

from voila.constants import ANALYSIS_DELTAPSI
import collections


def deltapsi(args):
    return pipeline_run(DeltaPsi(args))


def deltapsi_quantification(list_of_lsv, chnk, process_conf, logger):

    logger = majiq_utils.get_logger("%s/%s.majiq.log" % (process_conf.outDir, chnk),
                                    silent=process_conf.silent, debug=process_conf.debug)

    logger.info("Quantifying LSVs PSI.. %s" % chnk)
    num_exp = [len(process_conf.files1), len(process_conf.files1)]

    f_list = [None, None]

    f_list[0] = majiq_io.get_extract_lsv_list(list_of_lsv, process_conf.files1)
    f_list[1] = majiq_io.get_extract_lsv_list(list_of_lsv, process_conf.files2)

    prior_matrix = np.array(load_bin_file(get_prior_matrix_filename(process_conf.outDir,
                                                                    process_conf.names)))

    for lidx, lsv_id in enumerate(list_of_lsv):
        if lidx % 50 == 0:
            print("Event %d ..." % lidx)
            sys.stdout.flush()
        lsv_samples = [None, None]

        for grp_idx in range(2):
            lsv_samples[grp_idx] = f_list[grp_idx][lidx].coverage
            lsv_type = f_list[grp_idx][lidx].type

        psi1, psi2 = [np.array(xx) for xx in lsv_samples]
        msamples = psi1.shape[2]
        del lsv_samples

        post_matrix, posterior_psi1, posterior_psi2, mu_psi1, mu_psi2 = deltapsi_posterior(psi1, psi2, prior_matrix,
                                                                                           msamples, num_exp,
                                                                                           process_conf.nbins,
                                                                                           lsv_type)

        qm = QueueMessage(QUEUE_MESSAGE_DELTAPSI_RESULT, (post_matrix, posterior_psi1, posterior_psi2,
                                                          mu_psi1, mu_psi2, lsv_id), chnk)
        process_conf.queue.put(qm, block=True)



prior_conf = collections.namedtuple('conf', 'iter plotpath breakiter names binsize')

class DeltaPsi(BasicPipeline):

    def run(self):
        self.deltapsi()

    def deltapsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """

        majiq_utils.create_if_not_exists(self.logger_path)
        logger = majiq_utils.get_logger("%s/deltapsi_majiq.log" % self.logger_path, silent=self.silent,
                                             debug=self.debug)

        logger.info("Majiq deltapsi v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("GROUP1: %s" % self.files1)
        logger.info("GROUP2: %s" % self.files2)
        self.nbins = 20

        files = [self.files1, self.files2]

        self.queue = mp.Queue()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]

        lsv_dict1, lsv_types1, lsv_summarized1, meta1, lsv_dict_graph1 = majiq_io.extract_lsv_summary(self.files1)
        list_of_lsv_group1 = majiq_filter.merge_files_hdf5(lsv_dict1, lsv_summarized1, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=logger)

        print(len(list_of_lsv_group1))

        lsv_dict2, lsv_types2, lsv_summarized2, meta2, lsv_dict_graph2 = majiq_io.extract_lsv_summary(self.files2)
        list_of_lsv_group2 = majiq_filter.merge_files_hdf5(lsv_dict2, lsv_summarized2, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=logger)
        print (len(list_of_lsv_group2))
        list_of_lsv = list(set(list_of_lsv_group1).intersection(set(list_of_lsv_group2)))

        psi_space, prior_matrix = gen_prior_matrix(lsv_dict1, lsv_summarized1, lsv_dict2, lsv_summarized2, lsv_types1,
                                                   self.outDir, prior_conf(self.iter, self.plotpath, self.breakiter,
                                                                           self.names, self.binsize),
                                                   numbins=self.nbins, defaultprior=self.default_prior,
                                                   minpercent=self.min_exp, logger=logger)

        logger.info("Saving prior matrix for %s..." % self.names)
        dump_bin_file(prior_matrix, get_prior_matrix_filename(self.outDir, self.names))

        lchnksize = max(len(list_of_lsv)/self.nthreads, 1) + 1
        print(len(list_of_lsv))
        weights = [None, None]
        for grp_idx, fil in enumerate(files):
            #TODO: FIX WEIGHTS, maybe just return rho
            weights[grp_idx] = self.calc_weights(self.weights[grp_idx], fil, list_of_lsv, self.lock, lchnksize, self.queue,
                                                 self.names[grp_idx], logger=logger)

        self.weights = weights
        if len(list_of_lsv) > 0:

            pool = mp.Pool(processes=self.nthreads, initializer=process_conf, initargs=[deltapsi_quantification, self],
                           maxtasksperchild=1)
            [xx.acquire() for xx in self.lock]
            pool.map_async(process_wrapper, majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            with Voila(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
                out_h5p.add_genome(meta1['genome'])
                out_h5p.set_analysis_type(ANALYSIS_DELTAPSI)
                out_h5p.add_experiments(group_name=self.names[0], experiment_names=meta1['experiments'])
                out_h5p.add_experiments(group_name=self.names[1], experiment_names=meta2['experiments'])
                # out_h5p.add_metainfo(meta1['genome'], group1=self.names[0], experiments1=meta1['experiments'],
                #                      group2=self.names[1], experiments2=meta2['experiments'])

                in_h5p = h5py.File(files[0][0], 'r')
                queue_manager(in_h5p, out_h5p, self.lock, self.queue, num_chunks=self.nthreads, logger=logger,
                              list_of_lsv_graphics=lsv_dict_graph1)
                in_h5p.close()
            pool.join()

        logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                     self.names[1],
                                                                                                     self.outDir))
        logger.info("Alakazam! Done.")
