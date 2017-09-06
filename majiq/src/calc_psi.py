import multiprocessing as mp
import sys
import traceback

import h5py

import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.utils as majiq_utils
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper
from majiq.src.psi import psi_posterior
from voila.api import Voila
from voila.constants import ANALYSIS_PSI



################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    return pipeline_run(CalcPsi(args))


def psi_quantification(list_of_lsv, chnk, process_conf, logger):

    logger.info("Quantifying LSVs PSI.. %s" % chnk)
    f_list = majiq_io.get_extract_lsv_list(list_of_lsv, process_conf.files)
    for lidx, lsv_id in enumerate(list_of_lsv):
        if lidx % 50 == 0:
            print("Event %d ..." % lidx)
            sys.stdout.flush()

        psi = np.array(f_list[lidx].coverage) * process_conf.weights[:, None, None]
        mu_psi, post_psi = psi_posterior(psi, psi.shape[2], len(process_conf.files), process_conf.nbins,
                                         f_list[lidx].type)
        qm = QueueMessage(QUEUE_MESSAGE_PSI_RESULT, (post_psi, mu_psi, lsv_id), chnk)
        process_conf.queue.put(qm, block=True)


class CalcPsi(BasicPipeline):
    def run(self):
        self.calcpsi()

    def calcpsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """
        majiq_utils.create_if_not_exists(self.logger_path)

        logger = majiq_utils.get_logger("%s/psi_majiq.log" % self.logger_path, silent=self.silent, debug=self.debug)

        logger.info("")
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("Running Psi ...")
        logger.info("GROUP: %s" % self.files)
        self.nbins = 40

        self.queue = mp.Queue()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]

        lsv_dict, lsv_types, lsv_summarized, meta, lsv_dict_graph = majiq_io.extract_lsv_summary(self.files)

        list_of_lsv = majiq_filter.merge_files_hdf5(lsv_dict=lsv_dict, lsv_summarized=lsv_summarized,
                                                    minnonzero=self.minpos, min_reads=self.minreads,
                                                    percent=self.min_exp, logger=logger)
        lchnksize = max(len(list_of_lsv) / self.nthreads, 1) + 1
        weights = self.calc_weights(self.weights, self.files, list_of_lsv, self.lock, lchnksize, self.queue, self.name)
        self.weights = weights
        if len(list_of_lsv) > 0:
            pool = mp.Pool(processes=self.nthreads, initializer=process_conf, initargs=[psi_quantification, self],
                           maxtasksperchild=1)
            [xx.acquire() for xx in self.lock]

            pool.map_async(process_wrapper, majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            with Voila(get_quantifier_voila_filename(self.outDir, self.name), 'w') as out_h5p:
                out_h5p.add_genome(meta['genome'])
                out_h5p.set_analysis_type(ANALYSIS_PSI)
                out_h5p.add_experiments(group_name=self.name, experiment_names=meta['experiments'])

                in_h5p = h5py.File(self.files[0], 'r')
                queue_manager(in_h5p, out_h5p, self.lock, self.queue, num_chunks=self.nthreads,
                              list_of_lsv_graphics=lsv_dict_graph, logger=logger)
                in_h5p.close()

            pool.join()

        logger.info("PSI calculation for %s ended succesfully! "
                    "Result can be found at %s" % (self.name, self.outDir))
        logger.info("Alakazam! Done.")
