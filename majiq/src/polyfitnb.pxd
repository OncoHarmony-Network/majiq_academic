cimport numpy as np
from libcpp.vector cimport vector
cdef float fit_nb(vector[np.float32_t *] junctionl, int total_juncs, int eff_length, float nbdisp, object logger)