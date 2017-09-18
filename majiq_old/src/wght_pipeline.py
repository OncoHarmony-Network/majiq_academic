import multiprocessing as mp
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.utils as majiq_utils
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run


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
        majiq_utils.create_if_not_exists(self.logger_path)
        self.logger = majiq_utils.get_logger("%s/psi_majiq.log" % self.logger_path, silent=self.silent,
                                             debug=self.debug)

        self.logger.info("")
        self.logger.info("Command: %s" % self)
        self.logger.info("Running Psi ...")
        self.logger.info("GROUP: %s" % self.files)
        self.nbins = 40
        self.weights = 'auto'

        lock_arr = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()
        lsv_dict, lsv_types, lsv_summarized, meta = majiq_io.extract_lsv_summary(self.files)

        list_of_lsv = majiq_filter.merge_files_hdf5(lsv_dict=lsv_dict, lsv_summarized=lsv_summarized,
                                                    minnonzero=self.minpos, min_reads=self.minreads,
                                                    percent=self.min_exp, logger=self.logger)

        lchnksize = max(len(list_of_lsv)/self.nthreads, 1) + 1
        weights = self.calc_weights(self.weights, self.files, list_of_lsv, lock_arr, lchnksize, q,
                                    self.name, store=False)

        self.logger.info("Weights for %s are %s" %(self.name, weights))
        self.logger.info("Weights calculation ended succesfully!")
        self.logger.info("Alakazam! Done.")

