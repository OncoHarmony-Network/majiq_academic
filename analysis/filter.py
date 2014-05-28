"""
Functions to filter junction pairs by number of positions covered or number of reads
"""
import sys

from pylab import *
from polyfitnb import func2nb
from scipy.stats import nbinom
from numpy.ma import masked_less


def filter_bulk(matrix_filter, *matrices):
    ret = []
    for m in matrices:
        ret.append(m[matrix_filter])

    return ret


def filter_message(when, value, logger, junc):
    message = "%s (Filter=%s). %s"%(when, value, junc[0].shape)
    if logger:
        if type(logger) == bool:
            print message
        else: 
            logger.info(message)


def lsv_mark_stacks(lsv_list, fitfunc, pvalue_limit, dispersion, logger=False):
    a, b = fitfunc.c
    minstack = sys.maxint #the minimum value marked as stack
    numstacks = 0
    for lidx, junctions in enumerate(lsv_list):
        for i, junction in enumerate(junctions):
            for j, value in enumerate(junction):
                if value > 0:
                    #TODO Use masker, and marking stacks will probably be faster.
                    copy_junc = list(junction)
                    copy_junc.pop(j)
                    copy_junc = array(copy_junc)
                    copy_junc = copy_junc[copy_junc > 0]
                    #FINISH TODO
                    mean_rest = mean(copy_junc)
                    r, p = func2nb(a, b, mean_rest, dispersion)
                    my_nb = nbinom(r, 1-p)
                    pval = 1-my_nb.cdf(value)
                    if pval < pvalue_limit:
                        lsv_list[lidx][i, j] = -2 
                        minstack = min(minstack, value)
                        numstacks += 1
        masked_less(lsv_list[lidx], 0) #remask the stacks

    if logger: logger.info("Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"%(junctions.size, numstacks, pvalue_limit, (float(numstacks)/junctions.size)*100))

#TODO: (Jordi) I don't think this return is necessary.
    return lsv_list


def lsv_quantifiable ( list_lsv_tuple , minnonzero, min_reads, logger=False, fon = [True,True]):

    filter_message("Before quantifiable_filter", minnonzero, logger, array(list_lsv_tuple))
    filtered = []
    filtered_info = []
    for lsvdx, lsv in enumerate(list_lsv_tuple[0]):
        for idx in range(lsv.shape[0]):
            if ((not fon[1] or np.count_nonzero(lsv[idx]) >= minnonzero ) and
                (not fon[0] or lsv[idx].sum() >=  min_reads)):

                filtered.append(lsv)
                filtered_info.append(list_lsv_tuple[1][lsvdx])
                break
#    filter_message("After quantifiable_filter", minnonzero, logger, array(filtered))
    return filtered, filtered_info

def lsv_intersection( lsv_list1, lsv_list2 ):

    lsv_match = [[],[]]
    match_info = [[],[]]

    ids1 = set([xx[1] for xx in lsv_list1[1]])
    ids2 = set([xx[1] for xx in lsv_list2[1]])
    matched_names = ids1.intersection(ids2)

    for ii in matched_names:
        for idx, nm in enumerate(lsv_list1[1]):
            if nm[1] == ii:
                lsv_match[0].append(lsv_list1[0][idx])
                match_info.append(lsv_list1[1][idx])
                break
        for idx, nm in enumerate(lsv_list2[1]):
            if nm[1] == ii:
                lsv_match[1].append(lsv_list2[0][idx])
#                match_info.append(lsv_list2[1][idx])
                break

    return lsv_match, match_info

def discardhigh(max0=0, orfilter=True, logger=False, *junc):
    
    filter_message("Before discardhigh", max0, logger, junc)

    if orfilter:
        j1_filtered = [] 
        j2_filtered = [] 
        for i in range(j1.shape[0]):
            if ((j1[i] > 0).sum() < max0) or ((j2[i] > 0).sum() < max0):
                j1_filtered.append(j1[i])
                j2_filtered.append(j2[i])

        ret = array(j1_filtered), array(j2_filtered)

    else:
        raise NotImplemented

    filter_message("After discardhigh", max0, logger, junc)

    return ret

def discardminreads_and(incexcpairs, minreads=0, logger=False):
    """
    Discard events that have at least one inclusion/exclusion pair that have less than minreads reads at least in exclusion or inclusion 
    """
    all_pairs = []
    ret = []
    for exc, inc in incexcpairs:
        exc, inc = discardminreads(minreads, False, logger, True, None, array(exc), array(inc)) 
        all_pairs.append(inc)
        all_pairs.append(exc)

    #convert to numpy arrays
    all_pairs = [array(x) for x in all_pairs]

    return discardminreads(0, True, logger, False, None, *all_pairs)


def discardminreads(minreads=0, orfilter=True, logger=False, returnempty=False, names=None, *junc):
    """
    Given a collection of N experiments with matching junctions, discard junctions that have less reads than *minreads* flag. 
    With *orfilter*, keeps every junction if at least one passes the filter. Without it, all junctions should pass the filter.
    With *returnempty*, returns empty lines instead of no lines. Useful so we don't lose the match with another discard
    """

    filter_message("Before discardminreads", minreads, logger, junc)
    filtered_juncs = [[] for x in xrange(len(junc))]
    filtered_names = []
    for i in range(junc[0].shape[0]): #for all events
        if orfilter:
            already = False
            for junc_num, junction in enumerate(junc):
                if (junction[i].sum() >= minreads):
                    #if one passes the filter, add all to the resultset and break the loop, since at least one passes the threshold
                    for junc_num2, junction2 in enumerate(junc):
                        filtered_juncs[junc_num2].append(junction2[i])
                        already = True

                    if names:
                        filtered_names.append(names[i])

                    break

            if not already and returnempty: #special case for further filtering with discardminreads_and. The -100 values will be eliminated in a posterior filter (as they are always above 0)
                for junc_num2, junction2 in enumerate(junc):
                    filtered_juncs[junc_num2].append([-100]*junction2[i])                

        else: #AND filter: All should pass the filter
            for junc_num, junction in enumerate(junc):
                all_pass = True
                if not (junction[i].sum() >= minreads): #if one of them doesn't pass the filter, flag and break
                    all_pass = False
                    break

            if all_pass:
                for junc_num2, junction2 in enumerate(junc):
                    filtered_juncs[junc_num2].append(junction2[i]) 
                
                if names:
                    filtered_names.append(names[i])  

            elif returnempty: 
                for junc_num2, junction2 in enumerate(junc):
                    filtered_juncs[junc_num2].append([-100]*len(junction2[i]))


    ret = [array(x) for x in filtered_juncs]
    filter_message("After discardminreads", minreads, logger, ret)
    if names: #if reference names were provided, include them
        ret.append(list(filtered_names))

    return ret


def discardmaxreads(maxreads=0, orfilter=True, logger=False, *junc):
    filter_message("Before discardmaxreads", maxreads, logger, junc)

    filtered_juncs = [[] for x in xrange(len(junc))]
    if orfilter:
        for i in range(junc[0].shape[0]): #for all junctions
            for junc_num, junction in enumerate(junc):
                if (junction[i].sum() >= maxreads):
                    #if one passes the filter, add all to resultset
                    for junc_num, junction2 in enumerate(junc):
                        filtered_juncs[junc_num].append(junction2[i])
                    break
        ret = [array(x) for x in filtered_juncs]
    else:
        raise NotImplemented
    filter_message("After discardmaxreads", maxreads, logger, ret)

    return ret


def discardlow(min0=0, orfilter=True, logger=False, names=None, *junc):
    filter_message("Before discardlow", min0, logger, junc)
    filtered_juncs = [[] for x in xrange(len(junc))]
    filtered_names = []
    if orfilter:
        for i in range(junc[0].shape[0]): #for all junctions
            for junc_num, junction in enumerate(junc):
                if ((junction[i] > 0).sum() > min0):
                    #if one passes the filter, add all to resultset
                    for junc_num2, junction2 in enumerate(junc):
                        filtered_juncs[junc_num2].append(junction2[i])
                        
                    if names:
                        filtered_names.append(names[i])

                    break

        ret = [array(x) for x in filtered_juncs]

    else:    
        raise NotImplemented 

    filter_message("After discardlow", min0, logger, ret)

    if names:
        ret.append(list(filtered_names))

    return ret


def norm_junctions(junctions, gc_factors=None, gcnorm=False, trim=False, logger=False):
    "Junction normalization by GC content and/or trimming borders"
    if gcnorm:
        print "Normalizing by GC factor"
        junctions = junctions*gc_factors

    #number of zeros before
    if trim:
        numzeros = _numzeros(junctions) 
        print "Previous junction length: %s Total number of zeros: %s\n"%(junctions.shape[1], sum(numzeros))
        junctions = junctions[:,trim:junctions.shape[1]-trim] #trim the junctions according to the flag

    return junctions


def mark_stacks(junctions, fitfunc, pvalue_limit, dispersion, logger=False):
    a, b = fitfunc.c
    minstack = sys.maxint #the minimum value marked as stack
    numstacks = 0
    for i, junction in enumerate(junctions):
        for j, value in enumerate(junction):
            if value > 0:
                #TODO Use masker, and marking stacks will probably be faster.
                copy_junc = list(junction)
                copy_junc.pop(j)
                copy_junc = array(copy_junc)
                copy_junc = copy_junc[copy_junc > 0]
                #FINISH TODO
                mean_rest = mean(copy_junc)
                r, p = func2nb(a, b, mean_rest, dispersion)
                my_nb = nbinom(r, 1-p)
                pval = 1-my_nb.cdf(value)
                if pval < pvalue_limit:
                    junctions[i, j] = -2 
                    minstack = min(minstack, value)
                    numstacks += 1

    if logger: logger.info("Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"%(junctions.size, numstacks, pvalue_limit, (float(numstacks)/junctions.size)*100))

    return junctions
