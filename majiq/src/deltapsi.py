import collections
import multiprocessing as mp
import sys

import majiq.src.io as majiq_io
import psutil
from majiq.src.psi import deltapsi_posterior, gen_prior_matrix

import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper, chunks
from voila import constants
from voila.api import Matrix
from voila.constants import ANALYSIS_DELTAPSI


def deltapsi(args):
    return pipeline_run(DeltaPsi(args))


def deltapsi_quantification(list_of_lsv, chnk, conf, logger):
    logger.info("Quantifying LSVs Delta PSI.. %s" % chnk)
    num_exp = [len(conf.files1), len(conf.files1)]

    f_list = [None, None]

    f_list[0] = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files1)
    f_list[1] = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files2)

    if conf.weights[0] is None:
        weights1 = majiq_io.load_weights(list_of_lsv, conf.outDir, conf.names[0])
        weights2 = majiq_io.load_weights(list_of_lsv, conf.outDir, conf.names[1])
    else:
        weights1 = {xx: conf.weights[0] for xx in list_of_lsv}
        weights2 = {xx: conf.weights[1] for xx in list_of_lsv}

    prior_matrix = np.array(majiq_io.load_bin_file(get_prior_matrix_filename(conf.outDir, conf.names)))
    with Matrix(get_quantifier_voila_filename(conf.outDir, conf.names, deltapsi=True), lock=conf.lock_out,
                mode='a') as out_h5p:

        for lidx, lsv_id in enumerate(list_of_lsv):

            if lidx % 50 == 0:
                print("Event %d ..." % lidx)
                sys.stdout.flush()
            if f_list[0][lsv_id].coverage.shape[1] < 2:
                continue
            lsv_type = conf.lsv_type_dict[lsv_id]
            boots1 = f_list[0][lsv_id].coverage * weights1[lsv_id][:, None, None]
            boots2 = f_list[1][lsv_id].coverage * weights2[lsv_id][:, None, None]

            post_dpsi, post_psi1, post_psi2, mu_psi1, mu_psi2 = deltapsi_posterior(boots1, boots2, prior_matrix,
                                                                                   boots1.shape[2],
                                                                                   num_exp, conf.nbins, lsv_type)

            out_h5p.delta_psi(lsv_id).add(lsv_type=lsv_type, bins=post_dpsi, group_bins=[post_psi1, post_psi2],
                                          group_means=[mu_psi1, mu_psi2], junctions=conf.junc_info[lsv_id])

        # qm = QueueMessage(QUEUE_MESSAGE_DELTAPSI_RESULT, (post_dpsi, post_psi1, post_psi2,
        #                                                   mu_psi1, mu_psi2, lsv_id), chnk)
        # conf.queue.put(qm, block=True)


prior_conf = collections.namedtuple('conf', 'iter plotpath breakiter names binsize')


class DeltaPsi(BasicPipeline):

    def store_results(self, output, results, msg_type, extra={}):

        lsv_type = self.lsv_type_dict[results[5]]
        output.delta_psi(results[5]).add(lsv_type=lsv_type, bins=results[0],
                                         group_bins=[results[1], results[2]],
                                         group_means=[results[3], results[4]],
                                         junctions=extra['junc_info'][results[5]])

    def run(self):
        self.deltapsi()

    def deltapsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """

        majiq_logger.create_if_not_exists(self.outDir)
        logger = majiq_logger.get_logger("%s/deltapsi_majiq.log" % self.outDir, silent=self.silent,
                                         debug=self.debug)

        logger.info("Majiq deltapsi v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("GROUP1: %s" % self.files1)
        logger.info("GROUP2: %s" % self.files2)
        self.dpsi = True
        self.nbins = 20
        manager = mp.Manager()
        self.lsv_type_dict = manager.dict()
        self.junc_info = manager.dict()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]
        self.lock_out = mp.Lock()
        self.queue = manager.Queue()

        weights = [None, None]

        lsv_empirical_psi1 = {}
        # junc_info = {}
        list_of_lsv1, exps1 = majiq_io.extract_lsv_summary(self.files1, epsi=lsv_empirical_psi1,
                                                           types_dict=self.lsv_type_dict,
                                                           minnonzero=self.minpos, min_reads=self.minreads,
                                                           junc_info=self.junc_info, percent=self.min_exp, logger=logger)
        weights[0] = self.calc_weights(self.weights[0], list_of_lsv1, name=self.names[0], file_list=self.files1,
                                       logger=logger)
        logger.info("Group %s: %s LSVs" % (self.names[0], len(list_of_lsv1)))

        lsv_empirical_psi2 = {}
        list_of_lsv2, exps2 = majiq_io.extract_lsv_summary(self.files2, epsi=lsv_empirical_psi2,
                                                           types_dict=self.lsv_type_dict,
                                                           minnonzero=self.minpos, min_reads=self.minreads,
                                                           junc_info=self.junc_info, percent=self.min_exp, logger=logger)
        weights[1] = self.calc_weights(self.weights[1], list_of_lsv2, name=self.names[1], file_list=self.files2,
                                       logger=logger)

        logger.info("Group %s: %s LSVs" % (self.names[1], len(list_of_lsv1)))

        list_of_lsv = list(set(list_of_lsv1).intersection(set(list_of_lsv2)))
        logger.info("Number quantifiable LSVs: %s" % len(list_of_lsv))

        psi_space, prior_matrix = gen_prior_matrix(self.lsv_type_dict, lsv_empirical_psi1, lsv_empirical_psi2,
                                                   self.outDir, names=self.names, plotpath=self.plotpath,
                                                   iter=self.iter, binsize=self.binsize,
                                                   numbins=self.nbins, defaultprior=self.default_prior,
                                                   minpercent=self.min_exp, logger=logger)

        logger.info("Saving prior matrix for %s..." % self.names)
        majiq_io.dump_bin_file(prior_matrix, get_prior_matrix_filename(self.outDir, self.names))

        self.weights = weights
        if len(list_of_lsv) > 0:
            nthreads = min(self.nthreads, len(list_of_lsv))

            with Matrix(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
                out_h5p.file_version = constants.VOILA_FILE_VERSION
                out_h5p.analysis_type = ANALYSIS_DELTAPSI
                out_h5p.group_names = self.names
                out_h5p.prior = prior_matrix
                out_h5p.experiment_names = [exps1, exps2]
                # queue_manager(out_h5p, self.lock, self.queue, num_chunks=nthreads, func=self.store_results,
                #               logger=logger, junc_info=junc_info)

            pool = mp.Pool(processes=nthreads, initializer=process_conf, initargs=[deltapsi_quantification, self],
                           maxtasksperchild=1)
            [xx.acquire() for xx in self.lock]

            pool.imap_unordered(process_wrapper, chunks(list_of_lsv, nthreads))
            pool.close()
            pool.join()

        if self.mem_profile:
            mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
            logger.info("Max Memory used %.2f MB" % mem_allocated)

        logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                     self.names[1],
                                                                                                     self.outDir))
        logger.info("Alakazam! Done.")
