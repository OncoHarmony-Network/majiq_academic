import sys
import majiq.src.io as majiq_io
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

cdef _core_deltapsi(object self):
    """
    Given a file path with the junctions, return psi distributions.
    write_pickle indicates if a .pickle should be saved in disk
    """

    cdef dict junc_info = {}
    cdef dict lsv_type_dict = {}
    cdef object logger
    cdef int nbins = 40
    cdef bint is_ir
    cdef string lsv_id
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

    cdef vector[string] lsv_vec

    majiq_logger.create_if_not_exists(self.outDir)
    logger = majiq_logger.get_logger("%s/deltapsi_majiq.log" % self.outDir, silent=self.silent,
                                     debug=self.debug)

    logger.info("Majiq deltapsi v%s" % VERSION)
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

    # logger.info("Saving prior matrix for %s..." % self.names)
    # majiq_io.dump_bin_file(prior_matrix, get_prior_matrix_filename(self.outDir, self.names))

    # self.weights = weights

    for lsv in list_of_lsv:
        lsv_vec.push_back(lsv.encode('utf-8'))
    nlsv = len(list_of_lsv)
    if nlsv == 0:
        logger.info("There is no LSVs that passes the filters")
        return

    nthreads = min(self.nthreads, len(list_of_lsv))

    cov_dict1 = majiq_io.get_coverage_lsv(list_of_lsv, self.files1, "")
    cov_dict2 = majiq_io.get_coverage_lsv(list_of_lsv, self.files2, "")

    for i in prange(nlsv, nogil=True, num_threads=nthreads):
        lsv_id = lsv_vec[i]
        nways = cov_dict1[lsv_id].size()
        msamples = cov_dict1[lsv_id][0].size()
        with gil:
            o_mupsi_1 = np.zeros(shape=nways, dtype=np.float32)
            out_mupsi_d_1[lsv_id] = o_mupsi_1
            o_postpsi_1 = np.zeros(shape=(nways, nbins), dtype=np.float32)
            out_postpsi_d_1[lsv_id] = o_postpsi_1

            o_mupsi_2 = np.zeros(shape=nways, dtype=np.float32)
            out_mupsi_d_2[lsv_id] = o_mupsi_2
            o_postpsi_2 = np.zeros(shape=(nways, nbins), dtype=np.float32)
            out_postpsi_d_2[lsv_id] = o_postpsi_2

            o_postdeltapsi = np.zeros(shape=(nways, (nbins*2)-1), dtype=np.float32)
            out_postdpsi_d[lsv_id] = o_postdeltapsi

            print ('type', lsv_type_dict[lsv_id.decode('utf-8')], prior_matrix[1].dtype, prior_matrix[0].dtype)
            is_ir = 'i' in lsv_type_dict[lsv_id.decode('utf-8')]
            if is_ir:
                prior_m = prior_matrix[1]
            else:
                prior_m = prior_matrix[0]

        deltapsi_posterior(cov_dict1[lsv_id], cov_dict2[lsv_id], <np.float32_t *> prior_m.data,
                           <np.float32_t *> o_mupsi_1.data, <np.float32_t *> o_mupsi_1.data,
                           <np.float32_t *> o_postpsi_1.data, <np.float32_t *> o_postpsi_2.data,
                           <np.float32_t *> o_postdeltapsi.data, msamples, nways, nbins, is_ir)

    with Matrix(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_DELTAPSI
        out_h5p.group_names = self.names
        out_h5p.prior = prior_matrix
        out_h5p.experiment_names = [exps1, exps2]
        for lsv in list_of_lsv:
            lsv_id = lsv.encode('utf-8')
            out_h5p.delta_psi(lsv).add(lsv_type=lsv_type_dict[lsv], bins=out_postdpsi_d[lsv_id],
                                       group_bins=[out_postpsi_d_1[lsv_id], out_postpsi_d_2[lsv_id]],
                                       group_means=[out_mupsi_d_1[lsv_id], out_mupsi_d_2[lsv_id]],
                                       junctions=junc_info[lsv])


    if self.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        logger.info("Max Memory used %.2f MB" % mem_allocated)

    logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                 self.names[1],
                                                                                                 self.outDir))



# def deltapsi_quantification(list_of_lsv, chnk, conf, logger):
#     logger.info("Quantifying LSVs Delta PSI.. %s" % chnk)
#     num_exp = [len(conf.files1), len(conf.files1)]
#
#     f_list = [None, None]
#
#     f_list[0] = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files1)
#     f_list[1] = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files2)
#
#     if conf.weights[0] is None:
#         weights1 = majiq_io.load_weights(list_of_lsv, conf.outDir, conf.names[0])
#         weights2 = majiq_io.load_weights(list_of_lsv, conf.outDir, conf.names[1])
#     else:
#         weights1 = {xx: conf.weights[0] for xx in list_of_lsv}
#         weights2 = {xx: conf.weights[1] for xx in list_of_lsv}
#
#     prior_matrix = np.array(majiq_io.load_bin_file(get_prior_matrix_filename(conf.outDir, conf.names)))
#
#     for lidx, lsv_id in enumerate(list_of_lsv):
#         if lidx % 50 == 0:
#             print("Event %d ..." % lidx)
#             sys.stdout.flush()
#         if f_list[0][lsv_id].coverage.shape[1] < 2 or f_list[0][lsv_id].coverage.shape[1] != f_list[1][lsv_id].coverage.shape[1]:
#             logger.info("Skipping Incorrect LSV %s" % lsv_id)
#             continue
#
#         boots1 = f_list[0][lsv_id].coverage * weights1[lsv_id][:, None, None]
#         boots2 = f_list[1][lsv_id].coverage * weights2[lsv_id][:, None, None]
#
#         post_dpsi, post_psi1, post_psi2, mu_psi1, mu_psi2 = deltapsi_posterior(boots1, boots2, prior_matrix,
#                                                                                boots1.shape[2],
#                                                                                num_exp, conf.nbins,
#                                                                                conf.lsv_type_dict[lsv_id])
#         qm = QueueMessage(QUEUE_MESSAGE_DELTAPSI_RESULT, (post_dpsi, post_psi1, post_psi2,
#                                                           mu_psi1, mu_psi2, lsv_id), chnk)
#         conf.queue.put(qm, block=True)
#
#
# prior_conf = collections.namedtuple('conf', 'iter plotpath breakiter names binsize')


class DeltaPsi(BasicPipeline):

    def store_results(self, output, results, msg_type, extra={}):

        lsv_type = self.lsv_type_dict[results[5]]
        output.delta_psi(results[5]).add(lsv_type=lsv_type, bins=results[0],
                                         group_bins=[results[1], results[2]],
                                         group_means=[results[3], results[4]],
                                         junctions=extra['junc_info'][results[5]])

    def run(self):
        _core_deltapsi(self)

