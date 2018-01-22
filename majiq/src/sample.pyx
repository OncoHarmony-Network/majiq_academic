from numpy.random  import choice
import numpy as np
cimport numpy as np
from scipy.stats import nbinom, poisson
import cython
from quicksect import IntervalNode, Interval, IntervalTree

"""
Sampling from junctions using a Negative Binomial model.
"""
LIM = 100
PSEUDO = 0.0000000001  # EPSILON is too small for some calculations
# DTYPE_t = np.float
ctypedef np.float64_t DTYPE_t


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef inline float _sample_over_poisson(float r, float mu, int num_samples):
    return poisson.rvs(mu, size=num_samples).mean() + 1



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef inline float _sample_over_nb(float r, float mu, int num_samples):
    cdef float p
    p = r / (r + mu)
    return nbinom.rvs(r, p, size=num_samples).mean() + 1


cdef np.ndarray _sample_from_junctions(np.ndarray[DTYPE_t, ndim=2] junction_list, int m, int k,
                                       float fitted_one_over_r=0.0):
    """Given the filtered reads, bootstrap samples from every junction
    :param junction_list:
    :param m:
    :param k:
    :param discardzeros:
    :param trimborder:
    :param fitted_one_over_r:
    :return:

    """

    cdef np.ndarray[DTYPE_t, ndim=2] all_samples = np.zeros(shape=(junction_list.shape[0], m+2), dtype=np.float)
    cdef int npos_mult
    cdef int iternumber
    cdef float m_samples_means, r = 0
    cdef np.ndarray[DTYPE_t, ndim=1] junction
    cdef int i

    if fitted_one_over_r>0:
        r = 1 / fitted_one_over_r
        func = _sample_over_nb
    else:
        func = _sample_over_poisson


    for i in range(junction_list.shape[0]):
        junction = junction_list[i][junction_list[i] > 0]
        npos_mult = np.count_nonzero(junction)
        if npos_mult > 0:
            for mm in range(m):
                m_samples_means  = choice(junction, k).mean()
                all_samples[i, mm] = func(r, mu=m_samples_means, num_samples=k) * npos_mult
            all_samples[i, m] = junction_list[i].sum()
            all_samples[i, m+1] = npos_mult

    return all_samples



cpdef np.ndarray sample_from_junctions(np.ndarray[DTYPE_t, ndim=2] junction_list, int m, int k,
                                       float fitted_one_over_r=0.0):
    return _sample_from_junctions(junction_list, m, k, fitted_one_over_r=fitted_one_over_r)


cpdef dict create_lt(object all_genes):
    cdef dict td = dict()
    cdef dict xx

    for xx in all_genes.values():
        try:
            td[xx['chromosome']].add(xx['start'], xx['end'], (xx['id'], xx['strand']))
        except KeyError:
            td[xx['chromosome']] = IntervalTree()
            td[xx['chromosome']].add(xx['start'], xx['end'], (xx['id'], xx['strand']))

    return td