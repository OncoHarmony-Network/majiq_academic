import multiprocessing as mp
import majiq.src.io as majiq_io
import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *

################################
# PSI calculation pipeline     #
################################

def calc_weights(args):
    return pipeline_run(Weights(args))


class Weights(BasicPipeline):

    def run(self):
        self.weights_pip()

    def weights_pip(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """
        majiq_logger.create_if_not_exists(self.logger_path)
        logger = majiq_logger.get_logger("%s/weights_majiq.log" % self.logger_path, silent=self.silent,
                                        debug=self.debug)

        logger.info("Majiq weights v%s" % VERSION)
        logger.info("Command: %s" % self)
        logger.info("Running weights ...")
        logger.info("GROUP: %s" % self.files)

        self.nbins = 40
        self.weights = 'auto'

        manager = mp.Manager()
        self.lsv_type_dict = manager.dict()

        self.lock = [mp.Lock() for xx in range(self.nthreads)]
        self.queue = mp.Queue()
        junc_info = {}
        list_of_lsv = majiq_io.extract_lsv_summary(self.files, minnonzero=self.minpos, types_dict=self.lsv_type_dict,
                                                   min_reads=self.minreads, percent=self.min_exp, junc_info=junc_info,
                                                   logger=logger)

        weights = self.calc_weights(self.weights, list_of_lsv, name=self.name, file_list=self.files, logger=logger)

        logger.info("Weights for %s are %s" %(self.name, weights))
        logger.info("Weights calculation ended succesfully!")
        logger.info("Alakazam! Done.")

