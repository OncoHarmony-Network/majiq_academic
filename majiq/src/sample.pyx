import sys
from random import choice

import numpy as np
cimport numpy as np
import majiq.src.polyfitnb as majiq_fit
from scipy.stats import nbinom, poisson
import cython

"""
Sampling from junctions using a Negative Binomial model.
"""
LIM = 100
PSEUDO = 0.0000000001  # EPSILON is too small for some calculations

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef inline float __calc_nbin_p(float r, float mu):
    p = r / (r + mu)
    return p


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray _sample_over_nb(float one_over_r, float mu, int num_samples):
    if one_over_r > 0:
        r = 1 / one_over_r
        p = __calc_nbin_p(r, mu)
        sampl = nbinom.rvs(r, p, size=num_samples)
    else:
        sampl = poisson.rvs(mu, size=num_samples)
    return sampl

cpdef np.ndarray sample_from_junctions(np.ndarray junction_list, int m, int k, float fitted_one_over_r=0.0):
    """Given the filtered reads, bootstrap samples from every junction
    :param junction_list:
    :param m:
    :param k:
    :param discardzeros:
    :param trimborder:
    :param fitted_one_over_r:
    :return:

    """

    cdef np.ndarray all_samples = np.zeros(shape=(junction_list.shape[0], m), dtype=np.float)
    cdef int npos_mult
    cdef int iternumber
    cdef float sampled_mean
    cdef np.ndarray junction, nb50
    cdef int i

    for i, junction in enumerate(junction_list):
        junction = junction[junction > 0]
        if np.count_nonzero(junction) > 0:
            npos_mult = np.count_nonzero(junction)
            for iternumber in range(m):
                sampled_mean = np.mean([choice(junction) for _ in range(k)])
                nb50 = _sample_over_nb(one_over_r=fitted_one_over_r, mu=sampled_mean, num_samples=k)
                all_samples[i, iternumber] = (np.mean(nb50) + 1) * npos_mult

        # all_samples = (all_samples + 1) * npos_mult

    return all_samples
