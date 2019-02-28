from libcpp.vector cimport vector
from majiq.src.internals.mtypes cimport *
cimport numpy as np

cdef extern from "qLSV.hpp":

    cdef cppclass qLSV:
        qLSV() nogil ;
        qLSV(int nways1, bint is_ir) nogil ;
        int get_num_ways() nogil ;
        void add(float* coverage, int msamples) nogil ;
        void reset_samps() nogil ;
        void set_bool( bint flt ) nogil ;
        bint is_enabled() nogil ;
        bint is_ir() nogil ;

    cdef cppclass psiLSV(qLSV):
        psiLSV() nogil ;
        psiLSV(int nways1, int nbins1, bint is_ir) nogil ;
        void clear() nogil ;

        psi_distr_t mu_psi ;
        vector[psi_distr_t] post_psi ;

    cdef cppclass dpsiLSV(qLSV):
        dpsiLSV() nogil ;
        dpsiLSV(int nways1, int nbins1, bint is_ir) nogil ;
        void add_condition1() nogil ;
        void add_condition2() nogil ;
        void clear() nogil ;
        void clear_all() nogil ;

        psi_distr_t mu_psi1 ;
        psi_distr_t mu_psi2 ;
        vector[psi_distr_t] post_psi1 ;
        vector[psi_distr_t] post_psi2 ;
        vector[psi_distr_t] post_dpsi ;

    cdef cppclass hetLSV(qLSV):
        hetLSV() nogil ;
        hetLSV(int nways1, int j_offset1, int max_nfiles1, int nbins1, bint is_ir, int cond1) nogil ;
        int get_junction_index() nogil ;
        void create_condition_samples(int size1, int size2, int psi_samples) nogil ;
        void add_condition1(float * k, int x, int nways, int psi_samples) nogil ;
        void add_condition2(float * k, int x, int nways, int psi_samples) nogil ;
        void clear() nogil ;

        vector[vector[psi_distr_t]] mu_psi ;
        vector[vector[psi_distr_t]] post_psi ;
