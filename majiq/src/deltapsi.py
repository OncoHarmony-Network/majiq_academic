import collections
import multiprocessing as mp
import sys

import majiq.src.io as majiq_io
from majiq.src.psi import deltapsi_posterior, gen_prior_matrix

import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper, chunks
from voila.api import Matrix
from voila.constants import ANALYSIS_DELTAPSI
import os

def deltapsi(args):
    return pipeline_run(DeltaPsi(args))


def deltapsi_quantification(list_of_lsv, chnk, conf, logger):
    logger.info("Quantifying LSVs Delta PSI.. %s" % chnk)
    num_exp = [len(conf.files1), len(conf.files1)]

    f_list = [None, None]

    f_list[0] = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files1)
    f_list[1] = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files2)

    prior_matrix = np.array(majiq_io.load_bin_file(get_prior_matrix_filename(conf.outDir, conf.names)))

    for lidx, lsv_id in enumerate(list_of_lsv):
        if lidx % 50 == 0:
            print("Event %d ..." % lidx)
            sys.stdout.flush()
        if f_list[0][lidx].coverage.shape[1] < 2:
            continue
        post_dpsi, post_psi1, post_psi2, mu_psi1, mu_psi2 = deltapsi_posterior(f_list[0][lidx].coverage,
                                                                               f_list[1][lidx].coverage, prior_matrix,
                                                                               f_list[0][lidx].coverage.shape[2],
                                                                               num_exp, conf.nbins,
                                                                               conf.lsv_type_dict[lsv_id])
        qm = QueueMessage(QUEUE_MESSAGE_DELTAPSI_RESULT, (post_dpsi, post_psi1, post_psi2,
                                                          mu_psi1, mu_psi2, lsv_id), chnk)
        conf.queue.put(qm, block=True)


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

        weights = [None, None]
        self.queue = mp.Queue()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]

        manager = mp.Manager()
        self.lsv_type_dict = manager.dict()

        lsv_empirical_psi1 = {}
        list_of_lsv1 = majiq_io.extract_lsv_summary(self.files1, epsi=lsv_empirical_psi1, types_dict=self.lsv_type_dict,
                                                    minnonzero=self.minpos, min_reads=self.minreads,
                                                    percent=self.min_exp, logger=logger)
        weights[0] = self.calc_weights(self.weights[0], self.files1, list_of_lsv1, self.names[0], logger=self.logger)
        logger.info("Group %s: %s LSVs" % (self.names[0], len(list_of_lsv1)))

        lsv_empirical_psi2 = {}
        list_of_lsv2 = majiq_io.extract_lsv_summary(self.files2, epsi=lsv_empirical_psi2, types_dict=self.lsv_type_dict,
                                                    minnonzero=self.minpos, min_reads=self.minreads,
                                                    percent=self.min_exp, logger=logger)
        weights[1] = self.calc_weights(self.weights[1], self.files2, list_of_lsv2, self.names[1], logger=self.logger)

        logger.info("Group %s: %s LSVs" % (self.names[1], len(list_of_lsv1)))

        list_of_lsv = list(set(list_of_lsv1).intersection(set(list_of_lsv2)))
        logger.info("Number quantifiable LSVs: %s" % len(list_of_lsv))
        # assert meta1['m_samples'] == meta2['m_samples'], \
        #     "Groups have different number of bootstrap samples(%s,%s)" % (meta1['m_samples'], meta2['m_samples'])

        psi_space, prior_matrix = gen_prior_matrix(self.lsv_type_dict, lsv_empirical_psi1, lsv_empirical_psi2,
                                                   self.outDir, names=self.names, breakiter=self.breakiter,
                                                   plotpath=self.plotpath, iter=self.iter, binsize=self.binsize,
                                                   numbins=self.nbins, defaultprior=self.default_prior,
                                                   minpercent=self.min_exp, logger=logger)

        logger.info("Saving prior matrix for %s..." % self.names)
        majiq_io.dump_bin_file(prior_matrix, get_prior_matrix_filename(self.outDir, self.names))

        self.weights = weights
        if len(list_of_lsv) > 0:
            nthreads = min(self.nthreads, len(list_of_lsv))
            pool = mp.Pool(processes=nthreads, initializer=process_conf, initargs=[deltapsi_quantification, self],
                           maxtasksperchild=1)
            [xx.acquire() for xx in self.lock]

            pool.imap_unordered(process_wrapper, chunks(list_of_lsv, nthreads))
            pool.close()

            with Matrix(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
                out_h5p.analysis_type = ANALYSIS_DELTAPSI
                out_h5p.group_names = self.names
                exps1 = [os.path.splitext(os.path.basename(xx))[0] for xx in self.files1]
                exps2 = [os.path.splitext(os.path.basename(xx))[0] for xx in self.files2]
                out_h5p.experiment_names = [exps1, exps2]
                queue_manager(out_h5p, self.lock, self.queue, num_chunks=nthreads, logger=logger,
                              lsv_type=self.lsv_type_dict)

            pool.join()

        logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                     self.names[1],
                                                                                                     self.outDir))
        logger.info("Alakazam! Done.")
