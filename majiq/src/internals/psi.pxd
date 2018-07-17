from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
cimport numpy as np

cdef extern from "psi.hpp":

    ctypedef vector[np.float32_t] psi_distr_t ;

    cdef void psi_posterior(vector[psi_distr_t] & i_psi, np.float32_t* o_mupsi, np.float32_t* o_postpsi,
                            int msamples, int njunc, int nbins, bint is_ir) nogil ;
    cdef void deltapsi_posterior(vector[psi_distr_t]& i_psi1, vector[psi_distr_t]& i_psi2, np.float32_t* prior_matrix,
                                 np.float32_t* o_mu_psi1, np.float32_t* o_mu_psi2, np.float32_t* o_post_psi1,
                                 np.float32_t* o_post_psi2, np.float32_t* o_posterior_dpsi,
                                 int msamples, int njunc, int nbins, bint is_ir) nogil ;
