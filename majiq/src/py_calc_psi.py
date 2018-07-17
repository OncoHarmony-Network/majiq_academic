import multiprocessing as mp
import sys

import majiq.src.io as majiq_io
import psutil

import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper, chunks
from majiq.src.psi import psi_posterior
from voila import constants
from voila.api import Matrix
from voila.constants import ANALYSIS_PSI


################################
# PSI calculation pipeline     #
################################


def calcpsi(args):
    return pipeline_run(CalcPsi(args))


def psi_quantification(list_of_lsv, chnk, conf, logger):
    logger.info("Quantifying LSVs PSI.. %s" % chnk)
    f_list = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files)

    if conf.weights is None:
        weights = majiq_io.load_weights(list_of_lsv, conf.outDir, conf.name)
    else:
        weights = {xx: conf.weights for xx in list_of_lsv}

    for lidx, lsv_id in enumerate(list_of_lsv):
        if lidx % 50 == 0:
            print("Event %d ..." % lidx)
            sys.stdout.flush()

        psi = f_list[lsv_id].coverage * weights[lsv_id][:, None, None]

        mu_psi, post_psi = psi_posterior(psi, psi.shape[2], len(conf.files), conf.nbins, conf.lsv_type_dict[lsv_id])
        qm = QueueMessage(QUEUE_MESSAGE_PSI_RESULT, (post_psi, mu_psi, lsv_id), chnk)
        conf.queue.put(qm, block=True)


class CalcPsi(BasicPipeline):

    def store_results(self, output, results, msg_type, extra={}):

        lsv_type = self.lsv_type_dict[results[2]]
        output.psi(results[2]).add(lsv_type=lsv_type, bins=results[0], means=results[1],
                                   junctions=extra['junc_info'][results[2]])

    def run(self):
        self.calcpsi()

    def calcpsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """
        majiq_logger.create_if_not_exists(self.outDir)

        logger = majiq_logger.get_logger("%s/psi_majiq.log" % self.outDir, silent=self.silent, debug=self.debug)

        logger.info("Majiq psi v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("Running Psi ...")
        logger.info("GROUP: %s" % self.files)

        manager = mp.Manager()
        self.lsv_type_dict = manager.dict()

        self.nbins = 40
        self.queue = manager.Queue()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]
        junc_info = {}
        list_of_lsv, exps = majiq_io.extract_lsv_summary(self.files, types_dict=self.lsv_type_dict,
                                                         minnonzero=self.minpos, min_reads=self.minreads,
                                                         percent=self.min_exp, junc_info=junc_info, logger=logger)

        self.weights = self.calc_weights(self.weights, list_of_lsv, name=self.name, file_list=self.files, logger=logger)

        if len(list_of_lsv) > 0:
            nthreads = min(self.nthreads, len(list_of_lsv))
            pool = mp.Pool(processes=nthreads, initializer=process_conf, initargs=[psi_quantification, self],
                           maxtasksperchild=1)
            [xx.acquire() for xx in self.lock]

            pool.map_async(process_wrapper, chunks(list_of_lsv, nthreads))
            pool.close()
            with Matrix(get_quantifier_voila_filename(self.outDir, self.name), 'w') as out_h5p:
                out_h5p.file_version = constants.VOILA_FILE_VERSION
                out_h5p.analysis_type = ANALYSIS_PSI
                out_h5p.experiment_names = [exps]
                out_h5p.group_names = [self.name]
                queue_manager(out_h5p, self.lock, self.queue, num_chunks=nthreads, func=self.store_results,
                              logger=logger, junc_info=junc_info)

            pool.join()

        if self.mem_profile:
            mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
            logger.info("Max Memory used %.2f MB" % mem_allocated)
        logger.info("PSI calculation for %s ended succesfully! "
                    "Result can be found at %s" % (self.name, self.outDir))
