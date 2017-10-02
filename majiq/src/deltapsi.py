import h5py
import sys
import collections
import multiprocessing as mp

import majiq.src.io as majiq_io
import majiq.src.deprecated_io as majiq_deprio
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.psi import deltapsi_posterior, gen_prior_matrix
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper, chunks
import majiq.src.logger as majiq_logger

from voila.api import Voila
from voila.constants import ANALYSIS_DELTAPSI


def deltapsi(args):
    return pipeline_run(DeltaPsi(args))


def deltapsi_quantification(list_of_lsv, chnk, process_conf, logger):

    logger.info("Quantifying LSVs PSI.. %s" % chnk)
    num_exp = [len(process_conf.files1), len(process_conf.files1)]

    f_list = [None, None]

    f_list[0] = majiq_deprio.get_extract_lsv_list(list_of_lsv, process_conf.files1)
    f_list[1] = majiq_deprio.get_extract_lsv_list(list_of_lsv, process_conf.files2)

    prior_matrix = np.array(majiq_io.load_bin_file(get_prior_matrix_filename(process_conf.outDir,
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

        majiq_logger.create_if_not_exists(self.logger_path)
        logger = majiq_logger.get_logger("%s/deltapsi_majiq.log" % self.logger_path, silent=self.silent,
                                        debug=self.debug)

        logger.info("Majiq deltapsi v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("GROUP1: %s" % self.files1)
        logger.info("GROUP2: %s" % self.files2)
        self.nbins = 20

        files = [self.files1, self.files2]

        self.queue = mp.Queue()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]

        lsv_empirical_psi1 = {}
        lsv_empirical_psi2 = {}
        meta1 = majiq_io.read_meta_info(self.files1)
        list_of_lsv1, lsv_dict_graph = majiq_io.extract_lsv_summary(self.files1, epsi=lsv_empirical_psi1,
                                                                    minnonzero=self.minpos, min_reads=self.minreads,
                                                                    percent=self.min_exp, logger=logger)
        logger.info("Group %s: %s LSVs" %(self.names[0], len(list_of_lsv1)))

        meta2 = majiq_io.read_meta_info(self.files2)
        list_of_lsv2, lsv_dict_graph = majiq_io.extract_lsv_summary(self.files2, epsi=lsv_empirical_psi2,
                                                                    minnonzero=self.minpos, min_reads=self.minreads,
                                                                    percent=self.min_exp, logger=logger)
        logger.info("Group %s: %s LSVs" %(self.names[1], len(list_of_lsv1)))

        list_of_lsv = list(set(list_of_lsv1).intersection(set(list_of_lsv2)))
        logger.info("Number quantifiable LSVs: %s" % len(list_of_lsv))

        psi_space, prior_matrix = gen_prior_matrix(lsv_dict_graph, lsv_empirical_psi1, lsv_empirical_psi2,
                                                   self.outDir, names=self.names, breakiter=self.breakiter,
                                                   plotpath=self.plotpath, iter=self.iter, binsize=self.binsize,
                                                   numbins=self.nbins, defaultprior=self.default_prior,
                                                   minpercent=self.min_exp, logger=logger)

        logger.info("Saving prior matrix for %s..." % self.names)
        majiq_io.dump_bin_file(prior_matrix, get_prior_matrix_filename(self.outDir, self.names))

        lchnksize = max(len(list_of_lsv)/self.nthreads, 1) + 1

        weights = [None, None]
        for grp_idx, fil in enumerate(files):
            #TODO: FIX WEIGHTS, maybe just return rho
            weights[grp_idx] = self.calc_weights(self.weights[grp_idx], fil, list_of_lsv, lchnksize,
                                                 self.names[grp_idx], logger=logger)

        self.weights = weights
        if len(list_of_lsv) > 0:

            pool = mp.Pool(processes=self.nthreads, initializer=process_conf, initargs=[deltapsi_quantification, self],
                           maxtasksperchild=1)
            [xx.acquire() for xx in self.lock]
            pool.map_async(process_wrapper, chunks(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            with Voila(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
                out_h5p.add_genome(meta1['genome'])
                out_h5p.set_analysis_type(ANALYSIS_DELTAPSI)
                out_h5p.add_experiments(group_name=self.names[0], experiment_names=meta1['experiments'])
                out_h5p.add_experiments(group_name=self.names[1], experiment_names=meta2['experiments'])

                in_h5p = h5py.File(files[0][0], 'r')
                queue_manager(in_h5p, out_h5p, self.lock, self.queue, num_chunks=self.nthreads, logger=logger,
                              list_of_lsv_graphics=lsv_dict_graph)
                in_h5p.close()
            pool.join()

        logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                     self.names[1],
                                                                                                     self.outDir))
        logger.info("Alakazam! Done.")
