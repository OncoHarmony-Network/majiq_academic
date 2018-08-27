from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
from majiq.src.internals.psi cimport psi_distr_t
cimport numpy as np


cdef extern from "io_utils.hpp":
    cdef void get_aggr_coverage(map[string, vector[psi_distr_t]]& output, string lsv_id, np.float32_t* coverage,
                                int njunc, int msamples) nogil ;
