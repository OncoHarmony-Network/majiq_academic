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
cdef np.ndarray _sample_over_nb_loop(float one_over_r, float mu, int num_samples):
    if one_over_r > 0:
        r = 1 / one_over_r
        p = r / (r + mu)
        sampl = nbinom.rvs(r, p, size=num_samples)
    else:
        sampl = poisson.rvs(mu, size=num_samples)
    return sampl

cdef np.ndarray _sample_from_junctions_old(np.ndarray junction_list, int m, int k, float fitted_one_over_r=0.0):
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
        npos_mult = np.count_nonzero(junction)
        if npos_mult > 0:
            for iternumber in range(m):
                sampled_mean = np.mean(choice(junction, k))
                nb50 = _sample_over_nb(one_over_r=fitted_one_over_r, mu=sampled_mean, num_samples=k)
                all_samples[i, iternumber] = (np.mean(nb50) + 1) * npos_mult

        # all_samples = (all_samples + 1) * npos_mult

    return all_samples


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray _sample_over_poisson(float r, float mu, int num_samples):
    return poisson.rvs(mu, size=num_samples).mean() +1



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef np.ndarray _sample_over_nb(float r, float mu, int num_samples):
    cdef np.ndarray p
    p = r / (r + mu)
    return nbinom.rvs(r, p, size=num_samples).mean()+1


cdef np.ndarray _sample_from_junctions(np.ndarray junction_list, int m, int k, float fitted_one_over_r=0.0):
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
    cdef float r = 0
    cdef np.ndarray junction, m_samples_means
    cdef int i

    if fitted_one_over_r>0:
        r = 1 / fitted_one_over_r
        func = _sample_over_nb
    else:
        func = _sample_over_poisson

    for i, junction in enumerate(junction_list):
        junction = junction[junction > 0]
        npos_mult = np.count_nonzero(junction)
        if npos_mult > 0:
            m_samples_means  = np.array([np.mean(choice(junction, k)) for _ in range(m)])
            all_samples[i] = np.array([func(r, mu=xx, num_samples=k) * npos_mult for xx in m_samples_means])


    return all_samples



cpdef np.ndarray sample_from_junctions(np.ndarray junction_list, int m, int k, float fitted_one_over_r=0.0):
    return _sample_from_junctions(junction_list, m, k, fitted_one_over_r=fitted_one_over_r)
