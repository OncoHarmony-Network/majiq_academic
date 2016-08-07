from majiq.src.utils import create_if_not_exists
import old_majiq.src.io as majiq_io
from majiq.src.polyfitnb import fit_nb
import abc
import numpy as np
from numpy.ma import masked_less


# ###############################
# Data loading and Boilerplate #
################################


def get_clean_raw_reads(matched_info, matched_lsv, outdir, names, num_exp):
    res = []
    for eidx in xrange(num_exp):
        for ldx, lsv in enumerate(matched_info):
            jlist = masked_less(matched_lsv[ldx][eidx], 0)

            num = jlist.sum(axis=1)
            res.append([lsv[1], num.data])
    #with open('%s/clean_reads.%s.pkl' % (outdir, names), 'wb') as fp:
    majiq_io.dump_bin_file(res, '%s/clean_reads.%s.pkl' % (outdir, names))


def pipeline_run(pipeline):
    """ Exception catching for all the pipelines """
    try:
        return pipeline.run()
    except KeyboardInterrupt:
        if pipeline.logger:
            pipeline.logger.info("MAJIQ manually interrupted. Avada kedavra...")


class BasicPipeline:
    def __init__(self, args):
        """Basic configuration shared by all pipelines"""

        #trick to dump argparse arguments into self
        self.__dict__.update(args.__dict__)
        create_if_not_exists(self.output)
        if self.plotpath:
            create_if_not_exists(self.plotpath)
        self.logger_path = self.logger
        if not self.logger_path:
            self.logger_path = self.output

        self.nthreads = args.nthreads
        self.psi_paths = []
        self.nchunks = self.nthreads #if args.nchunks == -1 else args.nchunks
        try:
            self.replica_len = [len(self.files1), len(self.files2)]
        except AttributeError:
            pass

    def fitfunc(self, const_junctions):
        """Calculate the Negative Binomial function to sample from using the Constitutive events"""
        if self.debug:
            self.logger.debug("Skipping fitfunc because --debug!")
            return np.poly1d([1, 0])
        else:
            self.logger.debug("Fitting NB function with constitutive events...")
            return fit_nb(const_junctions, "%s/nbfit" % self.output, self.plotpath, logger=self.logger)


    @abc.abstractmethod
    def run(self):
        """This is the entry point for all pipelines"""
        return