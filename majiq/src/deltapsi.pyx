import sys
import majiq.src.io as majiq_io
cimport majiq.src.io as majiq_io
import psutil

import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.internals.psi cimport deltapsi_posterior, psi_distr_t
from majiq.src.psi import gen_prior_matrix

from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from cython.parallel import prange

from voila.api import Matrix
from voila.constants import ANALYSIS_DELTAPSI, VOILA_FILE_VERSION
cimport numpy as np
import numpy as np
# import collections

def deltapsi(args):
    return pipeline_run(DeltaPsi(args))

cdef void _core_deltapsi(object self):

    cdef dict junc_info = {}
    cdef dict lsv_type_dict = {}
    cdef object logger
    cdef int nbins = 20
    cdef bint is_ir
    cdef string lsv
    cdef int nways, msamples, i
    cdef list list_of_lsv

    cdef dict out_mupsi_d_1 = {}
    cdef dict out_postpsi_d_1 = {}
    cdef dict out_mupsi_d_2 = {}
    cdef dict out_postpsi_d_2 = {}
    cdef dict out_postdpsi_d = {}
    cdef map[string, vector[psi_distr_t]] cov_dict1
    cdef map[string, vector[psi_distr_t]] cov_dict2
    cdef np.ndarray[np.float32_t, ndim=1, mode="c"] o_mupsi_1
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] o_postpsi_1
    cdef np.ndarray[np.float32_t, ndim=1, mode="c"] o_mupsi_2
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] o_postpsi_2
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] o_postdeltapsi
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] prior_m
    cdef list prior_matrix
    cdef map[string, int] lsv_vec


    majiq_logger.create_if_not_exists(self.outDir)
    logger = majiq_logger.get_logger("%s/deltapsi_majiq.log" % self.outDir, silent=self.silent,
                                     debug=self.debug)

    logger.info("Majiq deltapsi v%s-%s" % (VERSION, get_git_version()))
    logger.info("Command: %s" % " ".join(sys.argv))
    logger.info("GROUP1: %s" % self.files1)
    logger.info("GROUP2: %s" % self.files2)

    # weights = [None, None]

    lsv_empirical_psi1 = {}
    junc_info = {}
    list_of_lsv1, exps1 = majiq_io.extract_lsv_summary(self.files1, epsi=lsv_empirical_psi1,
                                                       types_dict=lsv_type_dict,
                                                       minnonzero=self.minpos, min_reads=self.minreads,
                                                       junc_info=junc_info, percent=self.min_exp, logger=logger)
    # weights[0] = self.calc_weights(self.weights[0], list_of_lsv1, name=self.names[0], file_list=self.files1,
    #                                logger=logger)
    logger.info("Group %s: %s LSVs" % (self.names[0], len(list_of_lsv1)))

    lsv_empirical_psi2 = {}
    list_of_lsv2, exps2 = majiq_io.extract_lsv_summary(self.files2, epsi=lsv_empirical_psi2,
                                                       types_dict=lsv_type_dict,
                                                       minnonzero=self.minpos, min_reads=self.minreads,
                                                       junc_info=junc_info, percent=self.min_exp, logger=logger)
    # weights[1] = self.calc_weights(self.weights[1], list_of_lsv2, name=self.names[1], file_list=self.files2,
    #                                logger=logger)

    logger.info("Group %s: %s LSVs" % (self.names[1], len(list_of_lsv1)))

    list_of_lsv = list(set(list_of_lsv1).intersection(set(list_of_lsv2)))
    logger.info("Number quantifiable LSVs: %s" % len(list_of_lsv))

    psi_space, prior_matrix = gen_prior_matrix(lsv_type_dict, lsv_empirical_psi1, lsv_empirical_psi2,
                                               self.outDir, names=self.names, plotpath=self.plotpath,
                                               iter=self.iter, binsize=self.binsize,
                                               numbins=nbins, defaultprior=self.default_prior,
                                               minpercent=self.min_exp, logger=logger)

    for lidx, lsv in enumerate(list_of_lsv):
        nways = lsv_type_dict[lsv][1]
        out_mupsi_d_1[lsv] = np.zeros(shape=nways, dtype=np.float32)
        out_postpsi_d_1[lsv] = np.zeros(shape=(nways, nbins), dtype=np.float32)
        out_mupsi_d_2[lsv] = np.zeros(shape=nways, dtype=np.float32)
        out_postpsi_d_2[lsv] = np.zeros(shape=(nways, nbins), dtype=np.float32)
        out_postdpsi_d[lsv] = np.zeros(shape=(nways, (nbins*2)-1), dtype=np.float32)
        lsv_vec[lsv] = nways

    # logger.info("Saving prior matrix for %s..." % self.names)
    # majiq_io.dump_bin_file(prior_matrix, get_prior_matrix_filename(self.outDir, self.names))

    # self.weights = weights

    nlsv = len(list_of_lsv)
    if nlsv == 0:
        logger.info("There is no LSVs that passes the filters")
        return

    nthreads = min(self.nthreads, nlsv)

    majiq_io.get_coverage_mat(cov_dict1, lsv_vec, self.files1, "", nthreads)
    majiq_io.get_coverage_mat(cov_dict2, lsv_vec, self.files2, "", nthreads)

    for i in prange(nlsv, nogil=True, num_threads=nthreads):
        with gil:
            lsv = list_of_lsv[i]
            nways = cov_dict1[lsv].size()
            msamples = cov_dict1[lsv][0].size()
            o_mupsi_1 = out_mupsi_d_1[lsv]
            o_postpsi_1 = out_postpsi_d_1[lsv]
            o_mupsi_2 = out_mupsi_d_2[lsv]
            o_postpsi_2 = out_postpsi_d_2[lsv]
            o_postdeltapsi = out_postdpsi_d[lsv]
            is_ir = b'i' in lsv_type_dict[lsv][0]

            # print ('type', lsv_type_dict[lsv_id.decode('utf-8')], prior_matrix[1].dtype, prior_matrix[0].dtype)
            if is_ir:
                prior_m = prior_matrix[1]
            else:
                prior_m = prior_matrix[0]

        deltapsi_posterior(cov_dict1[lsv], cov_dict2[lsv], <np.float32_t *> prior_m.data,
                           <np.float32_t *> o_mupsi_1.data, <np.float32_t *> o_mupsi_2.data,
                           <np.float32_t *> o_postpsi_1.data, <np.float32_t *> o_postpsi_2.data,
                           <np.float32_t *> o_postdeltapsi.data, msamples, nways, nbins, is_ir)

    logger.info('Computation done, saving results....')
    with Matrix(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_DELTAPSI
        out_h5p.group_names = self.names
        out_h5p.prior = prior_matrix
        out_h5p.experiment_names = [exps1, exps2]
        for lsv in list_of_lsv:
            out_h5p.delta_psi(lsv.decode('utf-8')).add(lsv_type=lsv_type_dict[lsv][0].decode('utf-8'),
                                                       bins=out_postdpsi_d[lsv],
                                                       group_bins=[out_postpsi_d_1[lsv], out_postpsi_d_2[lsv]],
                                                       group_means=[out_mupsi_d_1[lsv], out_mupsi_d_2[lsv]],
                                                       junctions=junc_info[lsv])


    if self.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        logger.info("Max Memory used %.2f MB" % mem_allocated)

    logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                 self.names[1],
                                                                                                 self.outDir))


class DeltaPsi(BasicPipeline):

    def store_results(self, output, results, msg_type, extra={}):

        lsv_type = self.lsv_type_dict[results[5]]
        output.delta_psi(results[5]).add(lsv_type=lsv_type, bins=results[0],
                                         group_bins=[results[1], results[2]],
                                         group_means=[results[3], results[4]],
                                         junctions=extra['junc_info'][results[5]])

    def run(self):
        _core_deltapsi(self)

