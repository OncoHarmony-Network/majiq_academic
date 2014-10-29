import abc
import pickle
from multiprocessing import Pool
from grimoire.utils.utils import create_if_not_exists, get_logger
from analysis.polyfitnb import fit_nb
import numpy as np
import analysis.filter as majiq_filter
import analysis.io as majiq_io
import analysis.psi as majiq_psi
import os
import builder as majiq_builder
from numpy.ma import masked_less
import pipe as pipe
# ###############################
# Data loading and Boilerplate #
################################


def get_clean_raw_reads(matched_info, matched_lsv, outdir, names, num_exp):
    res = []
    for eidx in xrange(num_exp):
        for ldx, lsv in enumerate(matched_info):
            jlist = masked_less(matched_lsv[ldx][eidx], 0)

            num = jlist.sum()
            res.append([lsv[1], num])

        with open('%s/clean_reads.%s%d.pkl' % (outdir, names, eidx), 'wb') as fp:
            pickle.dump(res, fp)


def _pipeline_run(pipeline):
    """ Exception catching for all the pipelines """
    try:
        return pipeline.run()
    except KeyboardInterrupt:
        if pipeline.logger:
            pipeline.logger.info("MAJIQ manually interrupted. Avada kedavra...")


def builder(args):
    majiq_builder.main(args)


class BasicPipeline:
    def __init__(self, args):
        """Basic configuration shared by all pipelines"""
        #trick to dump argparse arguments into self
        self.__dict__.update(args.__dict__)
        create_if_not_exists(self.output)
        if self.plotpath:
            create_if_not_exists(self.plotpath)
        logger_path = self.logger
        if not logger_path:
            logger_path = self.output

        self.logger = get_logger("%s/majiq.log" % logger_path, silent=self.silent, debug=self.debug)
        self.nthreads = args.nthreads
        self.nz = args.nz
        self.psi_paths = []
        try:
            self.replica_len = [len(self.files1), len(self.files2)]
        except AttributeError:
            pass

    @abc.abstractmethod
    def run(self, lsv):
        """This is the entry point for all pipelines"""
        return

    def gc_content_norm(self, lsv_list, const_list):
        """Normalize the matrix using the gc content"""
        self.logger.info("GC content normalization...")
        if self.gcnorm:
            for lidx, lsv in enumerate(lsv_list[0]):
                lsv_list[0][lidx] = np.multiply(lsv, lsv_list[2][lidx])
            const_list[0] = np.multiply(const_list[0], const_list[2])
        return lsv_list, const_list

    def fitfunc(self, const_junctions):
        """Calculate the Negative Binomial function to sample from using the Constitutive events"""
        if self.debug:
            self.logger.debug("Skipping fitfunc because --debug!")
            return np.poly1d([1, 0])
        else:
            self.logger.info("Fitting NB function with constitutive events...")
            return fit_nb(const_junctions, "%s/nbfit" % self.output, self.plotpath, logger=self.logger)

    def mark_stacks(self, lsv_list, fitfunc):
        if self.markstacks >= 0:
            self.logger.info("Marking and masking stacks for...")
            lsv_list = majiq_filter.lsv_mark_stacks(lsv_list, fitfunc, self.markstacks, self.nbdisp, self.logger)

        return lsv_list


################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    return _pipeline_run(CalcPsi(args))


class CalcPsi(BasicPipeline):
    def run(self):
        self.calcpsi()

    def calcpsi(self):
        """
        Given a file path with the junctions, return psi distributions. 
        write_pickle indicates if a .pickle should be saved in disk
        """
        self.logger.info("")
        self.logger.info("Running PSI groups...")
        self.logger.info("GROUP : %s" % self.files)

        num_exp = len(self.files)

        meta_info = [0] * num_exp

        filtered_lsv = [None] * num_exp
        fitfunc = [None] * num_exp
        for ii, fname in enumerate(self.files):
            meta_info[ii], lsv_junc, const = majiq_io.load_data_lsv(fname, self.name, self.logger)

            #fitting the function
            lsv_junc, const = self.gc_content_norm(lsv_junc, const)
            fitfunc[ii] = self.fitfunc(const[0])

            filtered_lsv[ii] = self.mark_stacks(lsv_junc, fitfunc[ii])
        matched_lsv, matched_info = majiq_filter.quantifiable_in_group(filtered_lsv, self.minnonzero, self.minreads,
                                                                       self.logger, 0.10)

        conf = {'minnonzero': self.minnonzero,
                'minreads': self.minreads,
                'm': self.m,
                'k': self.k,
                'discardzeros': self.discardzeros,
                'trimborder': self.trimborder,
                'debug': self.debug,
                'alpha': self.alpha,
                'n': self.n,
                'nbins': 40,
                'nz': self.nz}

        if self.nthreads == 1:
            posterior_matrix, names = pipe.calcpsi(matched_lsv, matched_info, num_exp, conf, fitfunc, self.logger)
        else:
            pool = Pool(processes=self.nthreads)
            csize = len(matched_lsv) / int(self.nthreads)
            self.logger.info("CREATING THREADS %s with <= %s lsv" % (self.nthreads, csize))

            for nt in xrange(self.nthreads):
                lb = nt * csize
                ub = min((nt + 1) * csize, len(matched_lsv))
                if nt == self.nthreads - 1:
                    ub = len(matched_lsv)
                lsv_list = matched_lsv[lb:ub]
                lsv_info = matched_info[lb:ub]
                pool.apply_async(pipe.parallel_lsv_child_calculation, [pipe.calcpsi,
                                                                       [lsv_list, lsv_info, num_exp, conf, fitfunc],
                                                                       lsv_info,
                                                                       '%s/tmp' % os.path.dirname(self.output),
                                                                       self.name,
                                                                       nt])
            pool.close()
            pool.join()

            posterior_matrix = []
            names = []
            self.logger.info("GATHER pickles")
            for nt in xrange(self.nthreads):
                tempfile = open("%s/tmp/%s_th%s.calcpsi.pickle" % (os.path.dirname(self.output), self.name, nt))
                ptempt = pickle.load(tempfile)
                posterior_matrix.extend(ptempt[0])
                names.extend(ptempt[1])

        pickle_path = "%s/%s_psigroup.pickle" % (self.output, self.name)
        pickle.dump([posterior_matrix, names, meta_info], open(pickle_path, 'w'))
        self.logger.info("PSI calculation for %s ended succesfully! Result can be found at %s" % (self.name,
                                                                                                  self.output))
        self.logger.info("Alakazam! Done.")


################################
#          Delta PSI           #
################################
def deltapair(args):
    _pipeline_run(DeltaPair(args))


class DeltaPair(BasicPipeline):
    def run(self):
        self.delta_groups()

    def delta_groups(self):
        self.logger.info("")
        self.logger.info("Running deltagroups new model ...")
        self.logger.info("GROUP 1: %s" % self.files1)
        self.logger.info("GROUP 2: %s" % self.files2)

        exec_id = '%s_%s' % (self.names[0], self.names[1])
        tempfile = '%s/%s_temp_mid_exec.pickle' % (self.output, exec_id)
        num_exp = [len(self.files1), len(self.files2)]
        meta_info = [[0] * num_exp[0], [0] * num_exp[1]]
        if not os.path.exists(tempfile):

            filtered_lsv1 = [None] * num_exp[0]
            fitfunc = [[None] * num_exp[0], [None] * num_exp[1]]
            for ii, fname in enumerate(self.files1):
                meta_info[0][ii], lsv_junc, const = majiq_io.load_data_lsv(fname, self.names[0], self.logger)

                #fitting the function
                lsv_junc, const = self.gc_content_norm(lsv_junc, const)
                fitfunc[0][ii] = self.fitfunc(const[0])
                filtered_lsv1[ii] = self.mark_stacks(lsv_junc, fitfunc[0][ii])
            filtered_lsv1 = majiq_filter.quantifiable_in_group(filtered_lsv1, self.minnonzero, self.minreads,
                                                               self.logger, 0.10)
            self.logger.info("Group1: %s quantifiable in group" % str(len(filtered_lsv1[0])))

            filtered_lsv2 = [None] * num_exp[1]
            for ii, fname in enumerate(self.files2):
                meta_info[1][ii], lsv_junc, const = majiq_io.load_data_lsv(fname, self.names[1], self.logger)

                #fitting the function
                lsv_junc, const = self.gc_content_norm(lsv_junc, const)
                fitfunc[1][ii] = self.fitfunc(const[0])
                filtered_lsv2[ii] = self.mark_stacks(lsv_junc, fitfunc[1][ii])
            filtered_lsv2 = majiq_filter.quantifiable_in_group(filtered_lsv2, self.minnonzero, self.minreads,
                                                               self.logger, 0.10)

            self.logger.info("Group2: %s quantifiable in group" % str(len(filtered_lsv2[0])))
            matched_lsv, matched_info = majiq_filter.lsv_intersection(filtered_lsv1, filtered_lsv2)
            self.logger.info("After intersection:  %d/(%d, %d)" % (len(matched_info), len(filtered_lsv1[0]),
                                                                   len(filtered_lsv2[0])))
            group1, group2 = pipe.combine_for_priormatrix(matched_lsv[0], matched_lsv[1], matched_info, num_exp)
            psi_space, prior_matrix = majiq_psi.gen_prior_matrix(self, group1, group2, self.output, numbins=20,
                                                                 defaultprior=self.default_prior)

            #TEMP 
            tout = open(tempfile, 'w+')
            pickle.dump([meta_info, matched_info, matched_lsv, psi_space, prior_matrix, fitfunc], tout)
            tout.close()
            #END TEMP

        else:
            meta_info, matched_info, matched_lsv, psi_space, prior_matrix, fitfunc = pickle.load(open(tempfile))

        conf = {'minnonzero': self.minnonzero,
                'minreads': self.minreads,
                'm': self.m,
                'k': self.k,
                'discardzeros': self.discardzeros,
                'trimborder': self.trimborder,
                'debug': self.debug,
                'alpha': self.alpha,
                'n': self.n,
                'plotpath': self.plotpath,
                'nz': self.nz,
                'names': self.names}

        get_clean_raw_reads(matched_info, matched_lsv[0], self.output, self.names[0], num_exp[0])
        get_clean_raw_reads(matched_info, matched_lsv[1], self.output, self.names[1], num_exp[1])

        if self.nthreads == 1:
            posterior_matrix, names = pipe.deltapsi(matched_lsv, matched_info, num_exp, conf, prior_matrix, fitfunc,
                                                    psi_space, self.logger)
        else:

            pool = Pool(processes=self.nthreads)
            csize = len(matched_lsv[0]) / int(self.nthreads)
            self.logger.info("CREATING THREADS %s with <= %s lsv" % (self.nthreads, csize))

            for nthrd in xrange(self.nthreads):
                lb = nthrd * csize
                ub = min((nthrd + 1) * csize, len(matched_lsv[0]))
                if nthrd == self.nthreads - 1:
                    ub = len(matched_lsv[0])
                lsv_list = [matched_lsv[0][lb:ub], matched_lsv[1][lb:ub]]
                lsv_info = matched_info[lb:ub]
                pool.apply_async(pipe.parallel_lsv_child_calculation, [pipe.deltapsi,
                                                                       [lsv_list, lsv_info, num_exp, conf, prior_matrix,
                                                                        fitfunc, psi_space],
                                                                       matched_info,
                                                                       '%s/tmp' % os.path.dirname(self.output),
                                                                       '%s_%s' % (self.names[0], self.names[1]),
                                                                       nthrd])
            pool.close()
            pool.join()

            posterior_matrix = []
            names = []
            self.logger.info("GATHER pickles")
            for nthrd in xrange(self.nthreads):
                tempfile = open("%s/tmp/%s_%s_th%s.%s.pickle" % (os.path.dirname(self.output), self.names[0],
                                                                 self.names[1], nthrd, pipe.deltapsi.__name__))
                ptempt = pickle.load(tempfile)
                posterior_matrix.extend(ptempt[0])
                names.extend(ptempt[1])

        pickle_path = "%s/%s_%s.%s.pickle" % (self.output, self.names[0], self.names[1], pipe.deltapsi.__name__)
        pickle.dump([posterior_matrix, names, meta_info], open(pickle_path, 'w'))
        self.logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                          self.names[1],
                                                                                                          self.output))
        self.logger.info("Alakazam! Done.")
