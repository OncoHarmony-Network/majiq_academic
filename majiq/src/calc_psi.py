import multiprocessing as mp
import sys
import traceback

import h5py

import majiq.src.io as majiq_io
from majiq.src.psi import psi_posterior

from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper, chunks
from voila.api import Matrix
from voila.constants import ANALYSIS_PSI
import majiq.src.logger as majiq_logger


################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    return pipeline_run(CalcPsi(args))


def psi_quantification(list_of_lsv, chnk, conf, logger):

    logger.info("Quantifying LSVs PSI.. %s" % chnk)
    f_list = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files, conf.m_samples)
    for lidx, lsv_id in enumerate(list_of_lsv):
        if lidx % 50 == 0:
            print("Event %d ..." % lidx)
            sys.stdout.flush()

        psi = np.array(f_list[lidx].coverage) * conf.weights[:, None, None]
        mu_psi, post_psi = psi_posterior(psi, psi.shape[2], len(conf.files), conf.nbins,
                                         conf.lsv_type_dict[lidx])
        qm = QueueMessage(QUEUE_MESSAGE_PSI_RESULT, (post_psi, mu_psi, lsv_id), chnk)
        conf.queue.put(qm, block=True)


class CalcPsi(BasicPipeline):
    def run(self):
        self.calcpsi()

    def calcpsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """
        majiq_logger.create_if_not_exists(self.logger_path)

        logger = majiq_logger.get_logger("%s/psi_majiq.log" % self.logger_path, silent=self.silent, debug=self.debug)

        logger.info("")
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("Running Psi ...")
        logger.info("GROUP: %s" % self.files)

        manager = mp.Manager()
        self.lsv_type_dict = manager.dict()

        self.nbins = 40
        self.queue = mp.Queue()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]

        meta = majiq_io.read_meta_info(self.files)
        self.m_samples = meta['m_samples']

        list_of_lsv = majiq_io.extract_lsv_summary(self.files, minnonzero=self.minpos, types_dict=self.lsv_types,
                                                   min_reads=self.minreads, percent=self.min_exp,
                                                   logger=logger)
        nthreads = min(self.nthreads, len(list_of_lsv))
        weights = self.calc_weights(self.weights, self.files, list_of_lsv, self.lock, self.queue, self.name)
        self.weights = weights

        if len(list_of_lsv) > 0:
            pool = mp.Pool(processes=nthreads, initializer=process_conf, initargs=[psi_quantification, self],
                           maxtasksperchild=1)
            [xx.acquire() for xx in self.lock]

            pool.map_async(process_wrapper, chunks(list_of_lsv, nthreads))
            pool.close()
            with Matrix(get_quantifier_voila_filename(self.outDir, self.name), 'w') as out_h5p:
                out_h5p.genome = meta['genome']
                out_h5p.analysis_type = ANALYSIS_PSI
                out_h5p.add_experiments(group_name=self.name, experiment_names=meta['experiments'])

                queue_manager(out_h5p, self.lock, self.queue, num_chunks=nthreads, logger=logger)

            pool.join()

        logger.info("PSI calculation for %s ended succesfully! "
                    "Result can be found at %s" % (self.name, self.outDir))
        logger.info("Alakazam! Done.")
