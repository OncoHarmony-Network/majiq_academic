import sys
import rna_majiq.src.io as majiq_io
cimport rna_majiq.src.io as majiq_io
import psutil

import rna_majiq.src.logger as majiq_logger
from rna_majiq.src.basic_pipeline import BasicPipeline, pipeline_run
import rna_majiq.src.constants as constants
from rna_majiq.src.internals.qLSV cimport dpsiLSV, qLSV
from rna_majiq.src.internals.mtypes cimport psi_distr_t
from rna_majiq.src.internals.psi cimport deltapsi_posterior, get_psi_border
from rna_majiq.src.internals.psi cimport gen_prior_matrix

from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from cython.parallel import prange

from rna_voila.api import Matrix
from rna_voila.constants import ANALYSIS_DELTAPSI, VOILA_FILE_VERSION
cimport numpy as np
import numpy as np
from rna_voila.api.licensing import check_license


def deltapsi(args):
    return pipeline_run(DeltaPsi(args))


cdef parse_dpsi_reads(map[string, qLSV*] lsv_map, list files1, list files2, int nthreads, int minreads, int minnonzero):

    cdef pair[string, qLSV*] p_qlsv
    cdef dpsiLSV* dpsiObj

    majiq_io.get_coverage_mat_lsv(lsv_map, files1, nthreads, True, minreads, minnonzero)
    for p_qlsv in lsv_map:
        dpsiObj = <dpsiLSV*> p_qlsv.second
        dpsiObj.add_condition1()

    majiq_io.get_coverage_mat_lsv(lsv_map, files2, nthreads, True, minreads, minnonzero)
    for p_qlsv in lsv_map:
        dpsiObj = <dpsiLSV*> p_qlsv.second
        dpsiObj.add_condition2()


cdef int _core_deltapsi(object self) except -1:

    cdef dict junc_info = {}
    cdef dict lsv_type_dict = {}
    cdef object logger
    cdef int nbins = 20
    cdef bint is_ir
    cdef string lsv
    cdef int nways, i
    cdef list list_of_lsv
    cdef bint voilafile, tsvfile

    cdef np.ndarray[np.float32_t, ndim=1, mode="c"] mupsi1
    cdef np.ndarray[np.float32_t, ndim=1, mode="c"] mupsi2
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] post_psi1
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] post_psi2
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] post_dpsi
    cdef vector[vector[psi_distr_t]] prior_matrix
    cdef vector[psi_distr_t] prior_m
    cdef psi_distr_t psi_border = psi_distr_t(nbins+1)
    cdef map[string, qLSV*] lsv_map
    cdef dict prior_cond = {'mreads': self.prior_minreads, 'mpos': self.prior_minnonzero}


    majiq_logger.create_if_not_exists(self.outDir)
    logFile = self.logger if self.logger else f"{self.outDir}/deltapsi_majiq.log"
    logger = majiq_logger.get_logger(logFile, silent=self.silent, debug=self.debug)
    logger.info(f"Majiq deltapsi v{constants.VERSION}")
    logger.info("Command: %s" % " ".join(sys.argv))
    check_license(self.license, logger)
    logger.info("GROUP1: %s" % self.files1)
    logger.info("GROUP2: %s" % self.files2)

    voilafile = self.output_type in ['all', 'voila']
    tsvfile = self.output_type in ['all', 'tsv']

    lsv_empirical_psi1 = {}
    junc_info = {}
    list_of_lsv1, exps1 = majiq_io.extract_lsv_summary(self.files1, epsi=lsv_empirical_psi1,
                                                       types_dict=lsv_type_dict,
                                                       minnonzero=self.minpos, min_reads=self.minreads,
                                                       junc_info=junc_info, prior_conf=prior_cond,
                                                       percent=self.min_exp, logger=logger)

    logger.info("Group %s: %s LSVs" % (self.names[0], len(list_of_lsv1)))
    lsv_empirical_psi2 = {}
    list_of_lsv2, exps2 = majiq_io.extract_lsv_summary(self.files2, epsi=lsv_empirical_psi2,
                                                       types_dict=lsv_type_dict,
                                                       minnonzero=self.minpos, min_reads=self.minreads,
                                                       junc_info=junc_info, prior_conf=prior_cond,
                                                       percent=self.min_exp, logger=logger)

    logger.info("Group %s: %s LSVs" % (self.names[1], len(list_of_lsv1)))

    list_of_lsv = list(set(list_of_lsv1).intersection(set(list_of_lsv2)))
    logger.info("Number quantifiable LSVs: %s" % len(list_of_lsv))
    prior_matrix = vector[vector[psi_distr_t]](2, vector[psi_distr_t](nbins, psi_distr_t(nbins, 0)))
    gen_prior_matrix(prior_matrix, lsv_type_dict, lsv_empirical_psi1, lsv_empirical_psi2, self.outDir,
                     names=self.names, iter=self.iter, binsize=self.binsize, numbins=nbins,
                     defaultprior=self.default_prior, minpercent=self.min_exp, logger=logger)


    nlsv = len(list_of_lsv)
    if nlsv == 0:
        logger.info("There is no LSVs that passes the filters")
        return 0


    for lsv in list_of_lsv:
        nways = lsv_type_dict[lsv][1]
        is_ir = b'i' in lsv_type_dict[lsv][0]
        m = new dpsiLSV(nways, nbins, is_ir)
        lsv_map[lsv] = <qLSV*> m


    nthreads = min(self.nthreads, nlsv)
    parse_dpsi_reads(lsv_map, self.files1, self.files2, nthreads, self.minreads, self.minpos)
    get_psi_border(psi_border, nbins)

    logger.info("Start Computing DeltaPSI")
    for i in prange(nlsv, nogil=True, num_threads=nthreads):

        with gil:
            lsv = list_of_lsv[i]

        if lsv_map[lsv].is_ir():
            prior_m = prior_matrix[1]
        else:
            prior_m = prior_matrix[0]

        deltapsi_posterior(<dpsiLSV*> lsv_map[lsv], prior_m, psi_border, nbins)


    logger.info('Computation done, saving results....')
    with Matrix(constants.get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w',
                voila_file=voilafile, voila_tsv=tsvfile) as out_h5p:

        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_DELTAPSI
        out_h5p.group_names = self.names
        out_h5p.experiment_names = [exps1, exps2]

        pmatrix = np.ndarray(shape=(nbins, nbins), dtype=np.float32, order="c")
        pmatrix_ir = np.ndarray(shape=(nbins, nbins), dtype=np.float32, order="c")
        for x in range(nbins):
            for y in range(nbins):
                pmatrix[x, y]    = prior_matrix[0][x][y]
                pmatrix_ir[x, y] = prior_matrix[1][x][y]
        out_h5p.prior = [pmatrix, pmatrix_ir]

        for lsv in list_of_lsv:
            dpsiObj_ptr = <dpsiLSV*> lsv_map[lsv]
            nways = dpsiObj_ptr.get_num_ways()
            mupsi1 = np.ndarray(shape=nways, dtype=np.float32, order="c")
            mupsi2 = np.ndarray(shape=nways, dtype=np.float32, order="c")
            postpsi1 = np.ndarray(shape=(nways, nbins), dtype=np.float32, order="c")
            postpsi2 = np.ndarray(shape=(nways, nbins), dtype=np.float32, order="c")
            postdpsi = np.ndarray(shape=(nways, (nbins*2)-1), dtype=np.float32, order="c")

            for x in range(nways):
                mupsi1[x] = dpsiObj_ptr.mu_psi1[x]
                mupsi2[x] = dpsiObj_ptr.mu_psi2[x]
                for y in range(nbins):
                    postpsi1[x, y] = dpsiObj_ptr.post_psi1[x][y]
                    postpsi2[x, y] = dpsiObj_ptr.post_psi2[x][y]

                for y in range((nbins*2)-1):
                    postdpsi[x, y] = dpsiObj_ptr.post_dpsi[x][y]
            dpsiObj_ptr.clear_all()
            out_h5p.delta_psi(lsv.decode('utf-8')).add(lsv_type=lsv_type_dict[lsv][0].decode('utf-8'),
                                                       bins=postdpsi, group_bins=[postpsi1,postpsi2],
                                                       group_means=[mupsi1, mupsi2], junctions=junc_info[lsv])

    if self.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        logger.info("Max Memory used %.2f MB" % mem_allocated)

    logger.info(f"DeltaPSI calculation for {self.names[0]}{constants.GROUP_NAME_SEP}{self.names[1]} ended successfully! Result can be found at {self.outDir}")


class DeltaPsi(BasicPipeline):

    def store_results(self, output, results, msg_type, extra={}):

        lsv_type = self.lsv_type_dict[results[5]]
        output.delta_psi(results[5]).add(lsv_type=lsv_type, bins=results[0],
                                         group_bins=[results[1], results[2]],
                                         group_means=[results[3], results[4]],
                                         junctions=extra['junc_info'][results[5]])

    def run(self):
        _core_deltapsi(self)

