__author__ = 'jordi@biociphers.org'
import sys
from itertools import izip

import numpy as np
from numpy.ma import masked_less
from scipy import interpolate
from scipy.stats.mstats_basic import mquantiles
import majiq.src.config as majiq_config
from majiq.src.constants import *
import majiq.src.checkpolyfitnb as majiqfit

#
# def mark_stacks(hdf5_fp, fitfunc_r, pvalue_limit, logger=None):
#     if pvalue_limit < 0:
#         return lsv_list
#     logger.debug("PCR amplification stacks normalization")
#     minstack = sys.maxint
#     # the minimum value marked as stack
#     numstacks = 0
#
#     for lidx, lsv in enumerate(hdf5_fp['LSVs']):
#         junctions = np.array(hdf5_fp[LSV_JUNCTIONS_DATASET_NAME][lsv['coverage']])
#         for i, junction in enumerate(junctions):
#             if np.count_nonzero(junction) == 0:
#                 continue
#             for j, value   in enumerate(junction):
#                 if value > 0:
#                     # TODO Use masker, and marking stacks will probably be faster.
#                     copy_junc = list(junction)
#                     copy_junc.pop(j)
#                     copy_junc = np.array(copy_junc)
#                     copy_junc = copy_junc[copy_junc > 0]
#                     nzpos = np.count_nonzero(copy_junc)
#
#                     #FINISH TODO
#                     mean_rest = np.mean(copy_junc) * nzpos
#                     pval = majiqfit.get_negbinom_pval(fitfunc_r, mean_rest, value)
#                     if pval < pvalue_limit:
#                         lsv_list[0][lidx][i, j] = -2
#                         minstack = min(minstack, value)
#                         numstacks += 1
#         masked_less(lsv_list[0][lidx], 0)


def mark_stacks(lsv_list, fitfunc_r, pvalue_limit, logger=None):
    import numpy.ma as ma
    if pvalue_limit < 0:
        return lsv_list
    logger.debug("Marking and masking stacks")
    minstack = sys.maxint
    # the minimum value marked as stack
    numstacks = 0
    for lidx, junctions in enumerate(lsv_list[0]):
        for i, junction in enumerate(junctions):
            if np.count_nonzero(junction) == 0:
                continue
            for j, value in enumerate(junction):
                if value > 0:
                    msk = junction <= 0
                    msk[j] = 1
                    copy_junc = ma.masked_array(junction, mask=msk)
                    mean_rest = 0.5 if copy_junc.mean() is ma.masked else copy_junc.mean() * copy_junc.count()
                    pval = majiqfit.get_negbinom_pval(fitfunc_r, mean_rest, value)
                    if pval < pvalue_limit:
                        lsv_list[0][lidx][i, j] = -2
                        minstack = min(minstack, value)
                        numstacks += 1
        masked_less(lsv_list[0][lidx], 0)

    # if logger:
    #     logger.debug("Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"
    #                   % (len(lsv_list[0]), numstacks, pvalue_limit, (float(numstacks) / len(lsv_list[0])) * 100))
    return lsv_list

def mark_stacks_per_junc(junction, fitfunc_r, pvalue_limit, logger=None):
    import numpy.ma as ma
    if logger is not None:
        logger.debug("Marking and masking stacks")
    minstack = sys.maxint
    # the minimum value marked as stack
    numstacks = 0

    for j, value in enumerate(junction):
        if value <= 0:
            continue
        msk = junction <= 0
        msk[j] = 1
        copy_junc = ma.masked_array(junction, mask=msk)
        mean_rest = 0.5 if copy_junc.mean() is ma.masked else copy_junc.mean() * copy_junc.count()
        pval = majiqfit.get_negbinom_pval(fitfunc_r, mean_rest, value)
        if pval < pvalue_limit:
            junction[j] = -2
            minstack = min(minstack, value)
            numstacks += 1
    masked_less(junction, 0)

    # if logger:
    #     logger.debug("Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"
    #                   % (len(lsv_list[0]), numstacks, pvalue_limit, (float(numstacks) / len(lsv_list[0])) * 100))
    return junction



def gc_normalization(output_gc_vals, logger=None):
    if logger is not None:
        logger.info("Gc Content normalization")
    # factor, meanbins = gc_factor_calculation(gc_pairs, nbins=10)
    # a = np.append(factor[exp_n], factor[exp_n][-1])
    factor, meanbins = output_gc_vals
    gc_factor = interpolate.interp1d(meanbins, factor, bounds_error=False, fill_value=1)
    return np.vectorize(gc_factor)


def gc_factor_calculation(gc_pairs, nbins=10):

    local_meanbins = np.zeros(shape=(nbins),   dtype=np.dtype('float'))
    local_factor = np.zeros(shape=(nbins),   dtype=np.dtype('float'))

    count = gc_pairs['COV']
    gc = gc_pairs['GC']

    if len(gc) == 0:
        local_factor = [1.0]*nbins
        local_meanbins = np.arange(0,1, 0.1)
        return local_factor, local_meanbins

    count, gc = izip(*sorted(izip(count, gc), key=lambda x: x[1]))

    num_regions = len(count)
    nperbin = num_regions / nbins

    quant_median = [0.0]*8
    bins = [0]*nbins

    for ii in range(nbins):
        lb = ii * nperbin
        if ii == nbins-1:
            ub = num_regions
        else:
            ub = (ii+1) * nperbin

        a = np.asarray(count[lb:ub])
        t = np.asarray(gc[lb:ub])

        local_meanbins[ii] = np.mean(t)
        bins[ii] = mquantiles(a, prob=np.arange(0.1, 0.9, 0.1))

    for qnt in range(8):
        qnt_bns = np.ndarray(len(bins))
        for idx, bb in enumerate(bins):
            qnt_bns[idx] = bb[qnt]
        quant_median[qnt] = np.mean(qnt_bns)

    for ii in range(nbins):
        offst = np.zeros(len(quant_median), dtype=np.dtype('float'))
        for idx, xx in enumerate(quant_median):
            offst[idx] = float(bins[ii][idx]) / float(xx)
        local_factor[ii] = 1/np.mean(offst)

    # local_meanbins[exp_n] = mean_bins
    # local_factor[exp_n] = gc_factor

    return local_factor, local_meanbins


def prepare_gc_content(gn):
    gc_pairs = {'GC': [[] for xx in xrange(majiq_config.num_experiments)],
                'COV': [[] for xx in xrange(majiq_config.num_experiments)]}

    for ex in gn.get_exon_list():
        gc_val = ex.get_gc_content()
        st, end = ex.get_coordinates()
        if gc_val == 0 or end - st < 30:
            continue
        for exp_n in xrange(majiq_config.num_experiments):
            cov = ex.get_coverage(exp_n)
            if cov < 1:
                continue
            gc_pairs['GC'][exp_n].append(gc_val)
            gc_pairs['COV'][exp_n].append(cov)

    return gc_pairs
