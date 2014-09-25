"""
Functions to filter junction pairs by number of positions covered or number of reads
"""
import sys
from polyfitnb import func2nb
from scipy.stats import nbinom
from numpy.ma import masked_less
import numpy as np


def filter_message(when, value, logger, junc):
    message = "%s (Filter=%s). %s" % (when, value, len(junc[0]))
    if logger:
        if type(logger) == bool:
            print message
        else: 
            logger.info(message)


def lsv_mark_stacks(lsv_list, fitfunc, pvalue_limit, dispersion, logger=None):
    a, b = fitfunc.c
    minstack = sys.maxint #the minimum value marked as stack
    numstacks = 0
    for lidx, junctions in enumerate(lsv_list[0]):
        for i, junction in enumerate(junctions):
            for j, value in enumerate(junction):
                if value > 0:
                    #TODO Use masker, and marking stacks will probably be faster.
                    copy_junc = list(junction)
                    copy_junc.pop(j)
                    copy_junc = np.array(copy_junc)
                    copy_junc = copy_junc[copy_junc > 0]
                    #FINISH TODO
                    mean_rest = np.mean(copy_junc)
                    r, p = func2nb(a, b, mean_rest, dispersion)
                    my_nb = nbinom(r, 1-p)
                    pval = 1-my_nb.cdf(value)
                    if pval < pvalue_limit:
                        lsv_list[0][lidx][i, j] = -2 
                        minstack = min(minstack, value)
                        numstacks += 1
        masked_less(lsv_list[0][lidx], 0)

    if logger:
        logger.info("Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"
                    % (junctions.size, numstacks, pvalue_limit, (float(numstacks)/junctions.size)*100))

#TODO: (Jordi) I don't think this return is necessary.
    return lsv_list


def quantifiable_in_group(list_of_experiments, minnonzero, min_reads, logger, per_exp=0.10):

    nexp = len(list_of_experiments)
    filtr = per_exp*float(nexp)
    filtr = 2
    if nexp <= 2:
        filtr = 1
    
    filt_exp = {}
    for idx, exp in enumerate(list_of_experiments):
        temp = lsv_quantifiable( exp, minnonzero, min_reads, logger)
        for ldx, lsv in enumerate(temp[1]):
            if not lsv[1] in filt_exp:
                filt_exp[lsv[1]] = 0
            filt_exp[lsv[1]] += 1

    tlb = {}
    filtered = []
    filtered_info = []

    for idx, exp in enumerate(list_of_experiments):
        for idx_lsv, lsv in enumerate(exp[1]):
            if not lsv[1] in tlb:
                tlb[lsv[1]] = [-1]*nexp
            tlb[lsv[1]][idx] = idx_lsv

    info = tlb.keys()
    for ii in info:
        if not ii in filt_exp:
            continue
        pres = filt_exp[ii]
        if pres < filtr:
            continue
        lsv = []
        id = list_of_experiments[0][1][tlb[ii][0]]
        for idx, exp in enumerate(list_of_experiments):
            local_indx = tlb[ii][idx]
#            pdb.set_trace()
            lsv.append(exp[0][local_indx])
        filtered.append(lsv)
        filtered_info.append(id)

    return filtered, filtered_info


def lsv_quantifiable(list_lsv_tuple, minnonzero, min_reads, logger=False, fon=[True, True], lsv_type='majiq'):

    filter_message("Before quantifiable_filter", minnonzero, logger, np.array(list_lsv_tuple))
    filtered = []
    filtered_info = []
    k = 2
    if lsv_type == 'majiq':
        for lsvdx, lsv in enumerate(list_lsv_tuple[0]):
            for idx in range(lsv.shape[0]):
                if ((not fon[1] or np.count_nonzero(lsv[idx]) >= minnonzero) and
                        (not fon[0] or lsv[idx].sum() >= min_reads)):
                    filtered.append(lsv)
                    filtered_info.append(list_lsv_tuple[1][lsvdx])
                    break

    elif lsv_type == 'miso':

        for lsvdx, lsv in enumerate(list_lsv_tuple[0]):
            total_count = lsv.sum()
            thresh_reads = min_reads + (lsv.shape[0] - 2)*k
            for idx in range(lsv.shape[0]):
                if ((not fon[1] or np.count_nonzero(lsv[idx]) >= minnonzero) and
                        (not fon[0] or total_count >= thresh_reads)):
                    filtered.append(lsv)
                    filtered_info.append(list_lsv_tuple[1][lsvdx])
                    break
    #    filter_message("After quantifiable_filter", minnonzero, logger, array(filtered))
    return filtered, filtered_info


def lsv_intersection(lsv_list1, lsv_list2):

    lsv_match = [[], []]
    match_info = []

    ids1 = set([xx[1] for xx in lsv_list1[1]])
    ids2 = set([xx[1] for xx in lsv_list2[1]])
    matched_names = ids1.intersection(ids2)

    for ii in matched_names:
        for idx, nm in enumerate(lsv_list1[1]):
            if nm[1] == ii:
                lsv_match[0].append(lsv_list1[0][idx])
                match_info.append(nm)
                break
        for idx, nm in enumerate(lsv_list2[1]):
            if nm[1] == ii:
                lsv_match[1].append(lsv_list2[0][idx])
#                match_info.append(lsv_list2[1][idx])
                break

    return lsv_match, match_info

