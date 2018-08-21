from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from majiq.src.internals.HetStats cimport HetStats

cimport numpy as np

cdef extern from "psi.hpp":

    ctypedef vector[np.float32_t] psi_distr_t ;
    ctypedef pair[int, int] pair_int_t

    cdef psi_distr_t get_psi_border(int nbins) nogil ;
    cdef void psi_posterior(vector[psi_distr_t] & i_psi, np.float32_t* o_mupsi, np.float32_t* o_postpsi,
                            int msamples, int njunc, int nbins, bint is_ir) nogil ;
    cdef void deltapsi_posterior(vector[psi_distr_t]& i_psi1, vector[psi_distr_t]& i_psi2, np.float32_t* prior_matrix,
                                 np.float32_t* o_mu_psi1, np.float32_t* o_mu_psi2, np.float32_t* o_post_psi1,
                                 np.float32_t* o_post_psi2, np.float32_t* o_posterior_dpsi,
                                 int msamples, int njunc, int nbins, bint is_ir) nogil ;
    cdef void get_samples_from_psi(vector[psi_distr_t]& i_psi, float* osamps, float* o_mu_psi, float* o_postpsi1,
                                   int psi_samples, int j_offset, psi_distr_t psi_border, int njunc,
                                   int msamples, int nbins) nogil ;

    cdef void test_calc(float* oPvals, vector[float*] samples1, vector[float*] samples2, HetStats* HetStatsObj,
               int njunc, int psamples) nogil ;
