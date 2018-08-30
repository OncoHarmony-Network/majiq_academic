from libcpp.vector cimport vector
cimport numpy as np

cdef extern from "qLSV.hpp":

    ctypedef vector[np.float32_t] psi_distr_t ;

    cdef cppclass qLSV:
        qLSV() nogil ;
        qLSV(int nways1, bint is_ir, int msamples) nogil ;
        int get_num_ways() nogil ;
        void add(float* coverage, int msamples) nogil ;

    cdef cppclass hetLSV(qLSV):
        hetLSV() nogil ;
        hetLSV(int nways1, int j_offset1, int max_nfiles1, int nbins1, bint is_ir, int cond1) nogil ;
        int get_junction_index() nogil ;
        vector[vector[psi_distr_t]] mu_psi ;
        vector[vector[psi_distr_t]] post_psi ;
