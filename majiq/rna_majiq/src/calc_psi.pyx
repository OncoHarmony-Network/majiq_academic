import sys
import rna_majiq.src.io as majiq_io
cimport rna_majiq.src.io as majiq_io
import psutil

import rna_majiq.src.logger as majiq_logger
from rna_majiq.src.basic_pipeline import BasicPipeline, pipeline_run
import rna_majiq.src.constants as constants
from rna_majiq.src.internals.mtypes cimport psi_distr_t
from rna_majiq.src.internals.psi cimport psi_posterior, get_psi_border
from rna_majiq.src.internals.qLSV cimport psiLSV, qLSV

from libcpp.string cimport string
from libcpp.map cimport map
from cython.parallel import prange, parallel, threadid

from rna_voila.api import Matrix
from rna_voila.constants import ANALYSIS_PSI, VOILA_FILE_VERSION
cimport numpy as np
import numpy as np
from rna_voila.api.licensing import check_license

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
    cdef object logger
    cdef int nbins = 40
    cdef dict junc_info = {}
    cdef dict lsv_type_dict = {}
    cdef bint is_ir
    cdef string lsv
    cdef int nways, msamples, i, cc, nthreads
    cdef np.ndarray[np.float32_t, ndim=1, mode="c"] mupsi
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] postpsi
    cdef list list_of_lsv
    cdef qLSV qlsvObj
    cdef bint voilafile, tsvfile
    cdef psiLSV* psiObj_ptr
    cdef map[string, qLSV*] lsv_map
    cdef psi_distr_t psi_border = psi_distr_t(nbins+1)

    majiq_logger.create_if_not_exists(self.outDir)

    logFile = self.logger if self.logger else f"{self.outDir}/psi_majiq.log"
    logger = majiq_logger.get_logger(logFile, silent=self.silent, debug=self.debug)
    logger.info(f"Majiq psi v{constants.VERSION}")
    logger.info("Command: %s" % " ".join(sys.argv))
    check_license(self.license, logger)
    logger.info("Running Psi ...")
    logger.info("GROUP: %s" % self.files)

    voilafile = self.output_type in ['all', 'voila']
    tsvfile = self.output_type in ['all', 'tsv']

    list_of_lsv, exps = majiq_io.extract_lsv_summary(self.files, types_dict=lsv_type_dict,
                                                     minnonzero=self.minpos, min_reads=self.minreads,
                                                     percent=self.min_exp, junc_info=junc_info, logger=logger)

    nlsv = len(list_of_lsv)
    if nlsv == 0:
        logger.info("There is no LSVs that passes the filters")
        return

    for lsv in list_of_lsv:
        nways = lsv_type_dict[lsv][1]
        is_ir = b'i' in lsv_type_dict[lsv][0]
        psiObj_ptr = new psiLSV(nways, nbins, is_ir)
        lsv_map[lsv] = <qLSV*> psiObj_ptr


    logger.info("Group %s: %s LSVs" % (self.name, nlsv))
    nthreads = min(self.nthreads, nlsv) + 1
    majiq_io.get_coverage_mat_lsv(lsv_map, self.files, nthreads, True, self.minreads, self.minpos)
    get_psi_border(psi_border, nbins)

    for i in prange(nlsv, nogil=True, num_threads=nthreads):

        with gil:
            lsv = list_of_lsv[i]

        psi_posterior(<psiLSV*> lsv_map[lsv], psi_border, nbins)

    logger.info('Computation done, saving results....')
    with Matrix(constants.get_quantifier_voila_filename(self.outDir, self.name), 'w',
                voila_file=voilafile, voila_tsv=tsvfile) as out_h5p:
        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_PSI
        out_h5p.experiment_names = [exps]
        out_h5p.group_names = [self.name]

        for lsv in list_of_lsv:
            psiObj_ptr = <psiLSV*> lsv_map[lsv]
            nways = psiObj_ptr.get_num_ways()
            mupsi = np.ndarray(shape=nways, dtype=np.float32, order="c")
            postpsi = np.ndarray(shape=(nways, nbins), dtype=np.float32, order="c")
            for x in range(nways):
                mupsi[x] = psiObj_ptr.mu_psi[x]
                for y in range(nbins):
                    postpsi[x, y] = psiObj_ptr.post_psi[x][y]

            psiObj_ptr.clear()
            # logger.info("Event %s" %(lsv_id))
            out_h5p.psi(lsv.decode('utf-8')).add(lsv_type=lsv_type_dict[lsv][0].decode('utf-8'), bins=postpsi,
                        means=mupsi, junctions=junc_info[lsv])

    if self.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        logger.info("Max Memory used %.2f MB" % mem_allocated)
    logger.info("PSI calculation for %s ended successfully! "
                "Result can be found at %s" % (self.name, self.outDir))

class CalcPsi(BasicPipeline):

    def run(self):
        _core_calcpsi(self)
