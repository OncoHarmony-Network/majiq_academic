import abc
import pickle
from multiprocessing import Pool, Process
import os
import gc

import numpy as np
from numpy.ma import masked_less

from majiq.src import builder as majiq_builder
from majiq.src.psi import combine_for_priormatrix

from majiq.src.utils.utils import create_if_not_exists, get_logger
from majiq.src.polyfitnb import fit_nb
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.psi as majiq_psi
import majiq.src.pipe as pipe


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

    with open('%s/clean_reads.%s.pkl' % (outdir, names), 'wb') as fp:
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
        self.logger_path = self.logger
        if not self.logger_path:
            self.logger_path = self.output


        self.nthreads = args.nthreads
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
            lsv_list = majiq_filter.lsv_mark_stacks(lsv_list, fitfunc, self.markstacks, self.logger)

        return lsv_list


################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    return _pipeline_run(CalcPsi(args))


class CalcPsi(BasicPipeline):
    def run(self):
        self.calcpsi()

    def pre_psi(self, nchunks, logger=None):

        if logger is None:
            logger = get_logger("%s/majiq.log" % self.output, silent=False)

        self.logger = logger

        num_exp = len(self.files)
        meta_info = [0] * num_exp
        filtered_lsv = [None] * num_exp
        fitfunc = [None] * num_exp
        for ii, fname in enumerate(self.files):
            meta_info[ii], lsv_junc, const = majiq_io.load_data_lsv(fname, self.name, logger)
            #fitting the function
            lsv_junc, const = self.gc_content_norm(lsv_junc, const)
            fitfunc[ii] = self.fitfunc(const[0])
            filtered_lsv[ii] = self.mark_stacks(lsv_junc, fitfunc[ii])

        matched_lsv, matched_info = majiq_filter.quantifiable_in_group(filtered_lsv, self.minpos, self.minreads,
                                                                       logger, 0.10)

        get_clean_raw_reads(matched_info, matched_lsv, self.output, self.name, num_exp)

        csize = len(matched_lsv) / nchunks
        outfdir = '%s/tmp/chunks/' % self.output
        if not os.path.exists(outfdir):
            os.makedirs(outfdir)

        logger.info("Saving meta info for %s..." % self.name)
        tout = open("%s/tmp/%s_metainfo.pickle" % (self.output, self.name), 'w+')
        pickle.dump(meta_info, tout)
        tout.close()

        logger.info("Creating %s chunks with <= %s lsv" % (nchunks, csize))
        for nthrd in xrange(nchunks):
            lb = nthrd * csize
            ub = min((nthrd + 1) * csize, len(matched_lsv))
            if nthrd == nchunks - 1:
                ub = len(matched_lsv)
            lsv_list = matched_lsv[lb:ub]
            lsv_info = matched_info[lb:ub]

            out_file = '%s/chunk_%d.pickle' % (outfdir, nthrd)
            tout = open(out_file, 'w+')
            pickle.dump([lsv_list, lsv_info, num_exp, fitfunc], tout)
            tout.close()

        gc.collect()

    def calcpsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """

        logger = get_logger("%s/majiq.log" % self.logger_path, silent=self.silent, debug=self.debug)
        logger.info("")
        logger.info("Running Psi ...")
        logger.info("GROUP: %s" % self.files)



        conf = {'minnonzero': self.minpos,
                'minreads': self.minreads,
                'm': self.m,
                'k': self.k,
                'discardzeros': self.discardzeros,
                'trimborder': self.trimborder,
                'debug': self.debug,
                'nbins': 40,
                'name': self.name}

        nchunks = self.nthreads
        if self.nthreads > 1:
            p = Process(target=self.pre_psi, args=[nchunks])
            p.start()
            p.join()

            pool = Pool(processes=self.nthreads)
        else:
            self.pre_psi(nchunks)

        for nthrd in xrange(nchunks):
            chunk_fname = '%s/tmp/chunks/chunk_%d.pickle' % (self.output, nthrd)
            if self.nthreads == 1:
                pipe.parallel_lsv_child_calculation(pipe.calcpsi,
                                                    [chunk_fname, conf],
                                                    '%s/tmp' % self.output,
                                                    self.name,
                                                    nthrd)

            else:
                pool.apply_async(pipe.parallel_lsv_child_calculation, [pipe.calcpsi,
                                                                       [chunk_fname, conf],
                                                                       '%s/tmp' % self.output,
                                                                       self.name,
                                                                       nthrd])

        if self.nthreads > 1:
            pool.close()
            pool.join()

        posterior_matrix = []
        names = []
        logger.info("GATHER pickles")
        for nt in xrange(self.nthreads):
            tempfile = open("%s/tmp/%s_th%s.calcpsi.pickle" % (self.output, self.name, nt))
            ptempt = pickle.load(tempfile)
            posterior_matrix.extend(ptempt[0])
            names.extend(ptempt[1])

        logger.info("Getting meta info for %s..." % self.name)
        tin = open("%s/tmp/%s_metainfo.pickle" % (self.output, self.name))
        meta_info = pickle.load(tin)
        tin.close()


        pickle_path = "%s/%s_psigroup.pickle" % (self.output, self.name)
        majiq_io.dump_lsvs_voila(pickle_path, posterior_matrix, names, meta_info)
        # pickle.dump([posterior_matrix, names, meta_info], open(pickle_path, 'w'))
        logger.info("PSI calculation for %s ended succesfully! Result can be found at %s" % (self.name, self.output))
        logger.info("Alakazam! Done.")


def prepare_lsvs(dpsi_obj, conf, nchunks, logger=None):

    if logger is None:
        logger = get_logger("%s/majiq.log" % dpsi_obj.output, silent=False)
    dpsi_obj.logger = logger
    exec_id = '%s_%s' % (dpsi_obj.names[0], dpsi_obj.names[1])
    #tempfile = '%s/%s_temp_mid_exec.pickle' % (dpsi_obj.output, exec_id)
    num_exp = [len(dpsi_obj.files1), len(dpsi_obj.files2)]

    filtered_lsv1 = [None] * num_exp[0]
    fitfunc = [[None] * num_exp[0], [None] * num_exp[1]]
    meta_info = [[0] * num_exp[0], [0] * num_exp[1]]

    for ii, fname in enumerate(dpsi_obj.files1):
        meta_info[0][ii], lsv_junc, const = majiq_io.load_data_lsv(fname, dpsi_obj.names[0], logger)

        #fitting the function
        lsv_junc, const = dpsi_obj.gc_content_norm(lsv_junc, const)
        fitfunc[0][ii] = dpsi_obj.fitfunc(const[0])
        filtered_lsv1[ii] = dpsi_obj.mark_stacks(lsv_junc, fitfunc[0][ii])
    filtered_lsv1 = majiq_filter.quantifiable_in_group(filtered_lsv1, dpsi_obj.minpos, dpsi_obj.minreads,
                                                       logger, 0.10)
    logger.info("Group1: %s quantifiable in group" % str(len(filtered_lsv1[0])))

    filtered_lsv2 = [None] * num_exp[1]
    for ii, fname in enumerate(dpsi_obj.files2):
        meta_info[1][ii], lsv_junc, const = majiq_io.load_data_lsv(fname, dpsi_obj.names[1], logger)

        #fitting the function
        lsv_junc, const = dpsi_obj.gc_content_norm(lsv_junc, const)
        fitfunc[1][ii] = dpsi_obj.fitfunc(const[0])
        filtered_lsv2[ii] = dpsi_obj.mark_stacks(lsv_junc, fitfunc[1][ii])
    filtered_lsv2 = majiq_filter.quantifiable_in_group(filtered_lsv2, dpsi_obj.minpos, dpsi_obj.minreads,
                                                       logger, 0.10)
    logger.info("Group2: %s quantifiable in group" % str(len(filtered_lsv2[0])))

    matched_lsv, matched_info = majiq_filter.lsv_intersection(filtered_lsv1, filtered_lsv2)
    logger.info("After intersection:  %d/(%d, %d)" % (len(matched_info), len(filtered_lsv1[0]),
                                                      len(filtered_lsv2[0])))

    group1, group2 = combine_for_priormatrix(matched_lsv[0], matched_lsv[1], matched_info, num_exp)
    psi_space, prior_matrix = majiq_psi.gen_prior_matrix(dpsi_obj, group1, group2, dpsi_obj.output, numbins=20,
                                                         defaultprior=dpsi_obj.default_prior)

    outfdir = '%s/tmp/chunks/' % dpsi_obj.output
    if not os.path.exists(outfdir):
        os.makedirs(outfdir)

    logger.info("Saving prior matrix for %s..." % dpsi_obj.names)
    tout = open("%s/%s_priormatrix.pickle" % (dpsi_obj.output, exec_id), 'w+')
    pickle.dump(prior_matrix, tout)
    tout.close()

    logger.info("Saving meta info for %s..." % dpsi_obj.names)
    tout = open("%s/tmp/%s_metainfo.pickle" % (dpsi_obj.output, exec_id), 'w+')
    pickle.dump(meta_info, tout)
    tout.close()

    get_clean_raw_reads(matched_info, matched_lsv[0], dpsi_obj.output, dpsi_obj.names[0], num_exp[0])
    get_clean_raw_reads(matched_info, matched_lsv[1], dpsi_obj.output, dpsi_obj.names[1], num_exp[1])

    csize = len(matched_lsv[0]) / nchunks

    logger.info("Creating %s chunks with <= %s lsv" % (nchunks, csize))
    for nthrd in xrange(nchunks):
        lb = nthrd * csize
        ub = min((nthrd + 1) * csize, len(matched_lsv[0]))
        if nthrd == nchunks - 1:
            ub = len(matched_lsv[0])
        lsv_list = [matched_lsv[0][lb:ub], matched_lsv[1][lb:ub]]
        lsv_info = matched_info[lb:ub]


        out_file = '%s/chunk_%d.pickle' % (outfdir, nthrd)
        tout = open(out_file, 'w+')
        pickle.dump([lsv_list, lsv_info, num_exp, conf, fitfunc, psi_space], tout)
        tout.close()

    gc.collect()



################################
#          Delta PSI           #
################################
def deltapair(args):
    _pipeline_run(DeltaPair(args))


class DeltaPair(BasicPipeline):
    def run(self):
        self.delta_groups()



    def delta_groups(self):
        logger = get_logger("%s/majiq.log" % self.logger_path, silent=self.silent, debug=self.debug)
        logger.info("")
        logger.info("Running deltagroups new model ...")
        logger.info("GROUP 1: %s" % self.files1)
        logger.info("GROUP 2: %s" % self.files2)

        exec_id = '%s_%s' % (self.names[0], self.names[1])

        num_exp = [len(self.files1), len(self.files2)]
        nchunks = int(self.nthreads)

        conf = {'minnonzero': self.minpos,
                'minreads': self.minreads,
                'm': self.m,
                'k': self.k,
                'discardzeros': self.discardzeros,
                'trimborder': self.trimborder,
                'debug': self.debug,
                'plotpath': self.plotpath,
                'names': self.names}

        p = Process(target=prepare_lsvs, args=(self, conf, nchunks))
        p.start()
        p.join()

        pool = Pool(processes=self.nthreads)
        for nthrd in xrange(nchunks):
            chunk_fname = '%s/tmp/chunks/chunk_%d.pickle' % (self.output, nthrd)
            delta_prior_path = "%s/%s_priormatrix.pickle" % (self.output, exec_id)
            if self.nthreads == 1:
                pipe.parallel_lsv_child_calculation(pipe.deltapsi,
                                                    [chunk_fname, delta_prior_path],
                                                    '%s/tmp' % self.output,
                                                    '%s_%s' % (self.names[0], self.names[1]),
                                                    nthrd)

            else:
                pool.apply_async(pipe.parallel_lsv_child_calculation, [pipe.deltapsi,
                                                                       [chunk_fname, delta_prior_path],
                                                                       '%s/tmp' % self.output,
                                                                       '%s_%s' % (self.names[0], self.names[1]),
                                                                       nthrd])

        if self.nthreads > 1:
            pool.close()
            pool.join()

        posterior_matrix = []
        names = []
        psi_list1 = []
        psi_list2 = []
        logger.info("GATHER pickles")
        for nthrd in xrange(self.nthreads):
            tempfile = open("%s/tmp/%s_%s_th%s.%s.pickle" % (self.output, self.names[0],
                                                             self.names[1], nthrd, pipe.deltapsi.__name__))
            ptempt = pickle.load(tempfile)
            posterior_matrix.extend(ptempt[0])
            names.extend(ptempt[1])
            psi_list1.extend(ptempt[2])
            psi_list2.extend(ptempt[3])

        logger.info("Getting meta info for %s..." % self.names)
        tin = open("%s/tmp/%s_metainfo.pickle" % (self.output, exec_id))
        meta_info = pickle.load(tin)
        tin.close()

        pickle_path = "%s/%s_%s.%s.pickle" % (self.output, self.names[0], self.names[1], pipe.deltapsi.__name__)
        # pickle.dump([posterior_matrix, names, meta_info, psi_list1, psi_list2], open(pickle_path, 'w'))
        majiq_io.dump_lsvs_voila(pickle_path, posterior_matrix, names, meta_info, psi_list1, psi_list2)
        logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                     self.names[1],
                                                                                                     self.output))

        # if self.pairwise:
        #     num_lsv = len(matched_info)
        #     for ii in np.arange(num_exp[0]):
        #         for jj in np.arange(num_exp[1]):
        #             names = ["%s_%d" % (self.names[0], ii+1), "%s_%d" % (self.names[1], jj+1)]
        #             logger.info("Pairwise deltapsi: %s vs %s" % (names[0], names[1]))
        #             conf['names'] = names
        #             fitting_f = [[fitfunc[0][ii]], [fitfunc[1][jj]]]
        #
        #
        #
        #             lsv_vals1 = [[matched_lsv[0][xx][ii]] for xx in range(num_lsv)]
        #             lsv_vals2 = [[matched_lsv[1][xx][jj]] for xx in range(num_lsv)]
        #
        #             lsv_vals = [lsv_vals1, lsv_vals2]
        #             self.pairwise_deltapsi(lsv_vals, matched_info, meta_info, [1, 1], conf, prior_matrix, fitting_f,
        #                                    psi_space, names)

        logger.info("Alakazam! Done.")

    def pairwise_deltapsi(self, matched_lsv, matched_info, meta_info, num_exp, conf, prior_matrix, fitfunc,
                          psi_space, grpnames):
        if self.nthreads == 1:
            posterior_matrix, names, psi_list1, psi_list2 = pipe.deltapsi(matched_lsv, matched_info, num_exp, conf,
                                                                          prior_matrix, fitfunc, psi_space, self.logger)
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
                                                                       '%s/tmp' % self.output,
                                                                       '%s_%s' % (grpnames[0], grpnames[1]),
                                                                       nthrd])
            pool.close()
            pool.join()

            posterior_matrix = []
            names = []
            psi_list1 = []
            psi_list2 = []
            self.logger.info("GATHER pickles")
            for nthrd in xrange(self.nthreads):
                tempfile = open("%s/tmp/%s_%s_th%s.%s.pickle" % (self.output, grpnames[0],
                                                                 grpnames[1], nthrd, pipe.deltapsi.__name__))
                ptempt = pickle.load(tempfile)
                posterior_matrix.extend(ptempt[0])
                names.extend(ptempt[1])
                psi_list1.extend(ptempt[2])
                psi_list2.extend(ptempt[3])

        pickle_path = "%s/%s_%s.%s.pickle" % (self.output, grpnames[0], grpnames[1], pipe.deltapsi.__name__)
        pickle.dump([posterior_matrix, names, meta_info, psi_list1, psi_list2], open(pickle_path, 'w'))
        self.logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (grpnames[0],
                                                                                                          grpnames[1],
                                                                                                          self.output))
