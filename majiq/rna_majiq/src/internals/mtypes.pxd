from libcpp.pair cimport pair
from libcpp.vector cimport vector
cimport numpy as np

ctypedef vector[np.float32_t] psi_distr_t ;
ctypedef pair[int, int] pair_int_t