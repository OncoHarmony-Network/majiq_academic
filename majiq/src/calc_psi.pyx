import sys
import majiq.src.io as majiq_io
cimport majiq.src.io as majiq_io
import psutil

import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.internals.psi cimport psi_posterior, psi_distr_t
from majiq.src.internals.grimoire cimport qLSV

from libcpp.string cimport string
from libcpp.map cimport map

from libcpp.vector cimport vector
from cython.parallel import prange, parallel, threadid

from voila.api import Matrix
from voila.constants import ANALYSIS_PSI, VOILA_FILE_VERSION
cimport numpy as np
import numpy as np

################################
# PSI calculation pipeline     #
################################


def calcpsi(args):
    return pipeline_run(CalcPsi(args))


cdef _core_calcpsi(object self):
    """
    Given a file path with the junctions, return psi distributions.
    write_pickle indicates if a .pickle should be saved in disk
    """
    cdef int nlsv
    cdef map[string, int] lsv_vec
    cdef map[string, vector[psi_distr_t]] cov_dict

    cdef object logger
    cdef int nbins = 40
    cdef dict junc_info = {}
    cdef dict lsv_type_dict = {}
    cdef dict out_mupsi_d = {}
    cdef dict out_postpsi_d = {}
    cdef bint is_ir
    cdef string lsv
    cdef int nways, msamples, i, loop_step, cc, nthreads
    cdef np.ndarray[np.float32_t, ndim=1, mode="c"] o_mupsi
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] o_postpsi
    cdef list list_of_lsv
    cdef int nchunks = 500
    cdef qLSV qlsvObj

    majiq_logger.create_if_not_exists(self.outDir)

    logger = majiq_logger.get_logger("%s/psi_majiq.log" % self.outDir, silent=self.silent, debug=self.debug)
    logger.info("Majiq psi v%s-%s" % (VERSION, get_git_version()))
    logger.info("Command: %s" % " ".join(sys.argv))
    logger.info("Running Psi ...")
    logger.info("GROUP: %s" % self.files)

    list_of_lsv, exps = majiq_io.extract_lsv_summary(self.files, types_dict=lsv_type_dict,
                                                     minnonzero=self.minpos, min_reads=self.minreads,
                                                     percent=self.min_exp, junc_info=junc_info, logger=logger)
    nlsv = len(list_of_lsv)
    if nlsv == 0:
        logger.info("There is no LSVs that passes the filters")
        return

    for lidx, lsv in enumerate(list_of_lsv):
        nways = lsv_type_dict[lsv][1]
        out_mupsi_d[lsv] = np.zeros(shape=nways, dtype=np.float32)
        out_postpsi_d[lsv] = np.zeros(shape=(nways, nbins), dtype=np.float32)
        lsv_vec[lsv] = nways

    # self.weights = self.calc_weights(self.weights, list_of_lsv, name=self.name, file_list=self.files, logger=logger)

    logger.info("Group %s: %s LSVs" % (self.name, nlsv))
    loop_step = max(1, int(nlsv/10))
    nthreads = min(self.nthreads, nlsv) + 1

    majiq_io.get_coverage_mat(cov_dict, lsv_vec, self.files, "", nthreads)
    logger.info("COV size: %s, %s" % (nlsv, cov_dict.size()))
    for i in prange(nlsv, nogil=True, num_threads=nthreads):

        with gil:
            lsv = list_of_lsv[i]
            # logger.info("2222 Event %s - %s" %(lsv, lsv_type_dict[lsv][0]))
            if i % loop_step == 0 :
                print ("Event %s/%s" %(i, nlsv))

            nways = cov_dict[lsv].size()
            msamples = cov_dict[lsv][0].size()
            o_mupsi = out_mupsi_d[lsv]
            o_postpsi = out_postpsi_d[lsv]
            is_ir = b'i' in lsv_type_dict[lsv][0]

        psi_posterior(cov_dict[lsv], <np.float32_t *> o_mupsi.data,
                      <np.float32_t *> o_postpsi.data, msamples, nways, nbins, is_ir)

    logger.info('Computation done, saving results....')
    with Matrix(get_quantifier_voila_filename(self.outDir, self.name), 'w') as out_h5p:
        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_PSI
        out_h5p.experiment_names = [exps]
        out_h5p.group_names = [self.name]
        for lsv in list_of_lsv:
            # logger.info("2222 Event %s" %(lsv_id))
            out_h5p.psi(lsv.decode('utf-8')).add(lsv_type=lsv_type_dict[lsv][0].decode('utf-8'), bins=out_postpsi_d[lsv],
                        means=out_mupsi_d[lsv], junctions=junc_info[lsv])

    if self.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        logger.info("Max Memory used %.2f MB" % mem_allocated)
    logger.info("PSI calculation for %s ended succesfully! "
                "Result can be found at %s" % (self.name, self.outDir))

class CalcPsi(BasicPipeline):

    def run(self):
        _core_calcpsi(self)