from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from majiq.src.internals.HetStats cimport HetStats
from majiq.src.internals.qLSV cimport dpsiLSV, hetLSV, qLSV, psiLSV

cimport numpy as np

cdef extern from "psi.hpp":

    ctypedef vector[np.float32_t] psi_distr_t ;
    ctypedef pair[int, int] pair_int_t

    cdef psi_distr_t& get_psi_border(psi_distr_t& psi_border, int nbins) nogil ;
    # cdef psi_distr_t& get_psi_border(int nbins) nogil ;
    cdef void psi_posterior(psiLSV* lsvObj, psi_distr_t& psi_border, int nbins) nogil ;

    cdef void deltapsi_posterior(dpsiLSV* lsvObj, vector[psi_distr_t]& prior_matrix, psi_distr_t& psi_border,
                                 int nbins) nogil ;

    cdef void get_samples_from_psi2(vector[psi_distr_t]& i_psi, np.float32_t* osamps, np.float32_t* o_mu_psi,
                                   np.float32_t* o_postpsi, int psi_samples, int j_offset, psi_distr_t& psi_border, int njunc,
                                   int msamples, int nbins, bint is_ir) nogil ;

    cdef void get_samples_from_psi(float* osamps, hetLSV* lsvObj, int psi_samples, psi_distr_t psi_border,
                                   int nbins, int cidx, int fidx) nogil ;

    cdef void get_samples_from_psi3(vector[psi_distr_t]& i_psi, vector[psi_distr_t]& osamps, psi_distr_t& o_mupsi,
                                   vector[psi_distr_t]& o_postpsi, int psi_samples, int j_offset,
                                   psi_distr_t psi_border, int njunc, int msamples, int nbins, bint is_ir) nogil ;

    cdef void test_calc(vector[psi_distr_t]& oPvals, HetStats* HetStatsObj, hetLSV* lsvObj, int psamples,
                        np.float32_t quant) nogil ;

    cdef void adjustdelta(psi_distr_t& o_mixtpdf, psi_distr_t& emp_dpsi, int num_iter, int nbins) nogil ;

    # cdef void test_calc(np.float32_t* oPvals, HetStats* HetStatsObj, hetLSV* lsvObj, int psamples,
    #                     np.float32_t quant) nogil ;
