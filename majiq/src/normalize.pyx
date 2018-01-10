__author__ = 'jordi@biociphers.org'
import sys
import numpy as np
cimport numpy as np
# from scipy import interpolate
# from scipy.stats.mstats_basic import mquantiles
# from majiq.src.config import Config
from majiq.src.constants import *
from majiq.src.polyfitnb cimport get_negbinom_pval

import numpy.ma as ma
import cython

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef int __mark_stacks_loop(np.ndarray[np.float_t, ndim=2] junctions, fitfunc_r, pvalue_limit, logger=None) except -1:
    cdef int lidx, i, j
    cdef float value, pval, mean_rest
    cdef np.ndarray[np.float_t, ndim=1] copy_junc
    cdef np.ndarray[np.int_t, ndim=1] msk

    if pvalue_limit < 0:
        return 0
    if logger:
        logger.debug("Marking and masking stacks")

    for i in range(junctions.shape[0]):
        if np.count_nonzero(junctions[i]) == 0:
            continue
        msk = junctions[i] <= 0
        for j in range(junctions.shape[1]):
            value = junctions[i, j]
            if value > 0:
                msk[j] = 1
                copy_junc = ma.masked_array(junctions[i], mask=msk)
                mean_rest = 0.5 if copy_junc.mean() is ma.masked else copy_junc.mean() * copy_junc.count()
                pval = get_negbinom_pval(fitfunc_r, mean_rest, value)
                if pval < pvalue_limit:
                    junctions[i, j] = 0


from scipy.stats import nbinom, poisson

cdef int __mark_stacks(np.ndarray[np.float_t, ndim=2] junctions, fitfunc_r, pvalue_limit, logger=None) except -1:

    cdef np.ndarray[np.float_t, ndim=2] pvalues
    cdef np.ndarray[np.float_t, ndim=2] mean_rest
    cdef np.ndarray[np.int_t, ndim=2] denom

    if pvalue_limit < 0:
        return 0
    if logger:
        logger.debug("Marking and masking stacks")
    denom = np.count_nonzero(junctions, axis=1)[:, None] - (junctions > 0)
    with np.errstate(divide='ignore',invalid='ignore'):
        mean_rest = (junctions.sum(axis=1)[:, None] - junctions) / denom
        mean_rest[np.isnan(mean_rest)] = 0.5
    if fitfunc_r >0:
        r = 1/fitfunc_r
        p = r/(mean_rest + r)
        pvalues = 1 - nbinom.cdf(junctions, r, p)
    else:
        pvalues = 1 - poisson.cdf(junctions, mean_rest)
    junctions[pvalues<pvalue_limit] = 0

cpdef mark_stacks(np.ndarray[np.float_t, ndim=2] junctions, fitfunc_r, pvalue_limit, logger=None):
    __mark_stacks(junctions, fitfunc_r, pvalue_limit, logger=logger)



#
# def gc_normalization(output_gc_vals, logger=None):
#     if logger is not None:
#         logger.info("Gc Content normalization")
#     # factor, meanbins = gc_factor_calculation(gc_pairs, nbins=10)
#     # a = np.append(factor[exp_n], factor[exp_n][-1])
#     factor, meanbins = output_gc_vals
#     gc_factor = interpolate.interp1d(meanbins, factor, bounds_error=False, fill_value=1)
#     return np.vectorize(gc_factor)
#
#
# def gc_factor_calculation(gc_pairs, nbins=10):
#
#     local_meanbins = np.zeros(shape=(nbins),   dtype=np.dtype('float'))
#     local_factor = np.zeros(shape=(nbins),   dtype=np.dtype('float'))
#
#     count = gc_pairs['COV']
#     gc = gc_pairs['GC']
#
#     if len(gc) == 0:
#         local_factor = [1.0]*nbins
#         local_meanbins = np.arange(0,1, 0.1)
#         return local_factor, local_meanbins
#
#     count, gc = zip(*sorted(zip(count, gc), key=lambda x: x[1]))
#
#     num_regions = len(count)
#     nperbin = num_regions / nbins
#
#     quant_median = [0.0]*8
#     bins = [0]*nbins
#
#     for ii in range(nbins):
#         lb = ii * nperbin
#         if ii == nbins-1:
#             ub = num_regions
#         else:
#             ub = (ii+1) * nperbin
#
#         a = np.asarray(count[lb:ub])
#         t = np.asarray(gc[lb:ub])
#
#         local_meanbins[ii] = np.mean(t)
#         bins[ii] = mquantiles(a, prob=np.arange(0.1, 0.9, 0.1))
#
#     for qnt in range(8):
#         qnt_bns = np.ndarray(len(bins))
#         for idx, bb in enumerate(bins):
#             qnt_bns[idx] = bb[qnt]
#         quant_median[qnt] = np.mean(qnt_bns)
#
#     for ii in range(nbins):
#         offst = np.zeros(len(quant_median), dtype=np.dtype('float'))
#         for idx, xx in enumerate(quant_median):
#             offst[idx] = float(bins[ii][idx]) / float(xx)
#         local_factor[ii] = 1/np.mean(offst)
#
#     # local_meanbins[exp_n] = mean_bins
#     # local_factor[exp_n] = gc_factor
#
#     return local_factor, local_meanbins
#
#
# def prepare_gc_content(gn):
#     majiq_config = Config()
#     gc_pairs = {'GC': [[] for xx in range(majiq_config.num_experiments)],
#                 'COV': [[] for xx in range(majiq_config.num_experiments)]}
#
#     for ex in gn.get_exon_list():
#         gc_val = ex.get_gc_content()
#         st, end = ex.get_coordinates()
#         if gc_val == 0 or end - st < 30:
#             continue
#         for exp_n in range(majiq_config.num_experiments):
#             cov = ex.get_coverage(exp_n)
#             if cov < 1:
#                 continue
#             gc_pairs['GC'][exp_n].append(gc_val)
#             gc_pairs['COV'][exp_n].append(cov)
#
#     return gc_pairs
