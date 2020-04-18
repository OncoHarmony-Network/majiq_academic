cimport numpy as np
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
cdef float fit_nb(vector[shared_ptr[vector[np.float32_t]]] junctionl, int total_juncs, int eff_length, float nbdisp, object logger)
