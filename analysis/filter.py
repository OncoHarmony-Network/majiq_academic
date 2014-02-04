"""
Functions to filter junction pairs by number of positions covered or number of reads
"""
import sys

from pylab import *
from polyfitnb import func2nb
from scipy.stats import nbinom


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
            logger.debug(message)


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
        exc, inc = discardminreads(minreads, False, logger, True, array(exc), array(inc)) 
        all_pairs.append(inc)
        all_pairs.append(exc)

    #convert to numpy arrays
    all_pairs = [array(x) for x in all_pairs]

    return discardminreads(0, True, logger, False, *all_pairs)


def discardminreads(minreads=0, orfilter=True, logger=False, returnempty=False, *junc):
    """
    Given a collection of N experiments with matching junctions, discard junctions that have less reads than *minreads* flag. 
    With *orfilter*, keeps every junction if at least one passes the filter. Without it, all junctions should pass the filter.
    With *returnempty*, returns empty lines instead of no lines. Useful so we don't lose the match with another discard
    """

    filter_message("Before discardminreads", minreads, logger, junc)
    filtered_juncs = [[] for x in xrange(len(junc))]

    for i in range(junc[0].shape[0]): #for all events
        if orfilter:
            already = False
            for junc_num, junction in enumerate(junc):
                if (junction[i].sum() >= minreads):
                    #if one passes the filter, add all to the resultset and break the loop, since at least one passes the threshold
                    for junc_num, junction2 in enumerate(junc):
                        filtered_juncs[junc_num].append(junction2[i])
                        already = True
                    break

            if not already and returnempty:
                for junc_num, junction2 in enumerate(junc):
                    filtered_juncs[junc_num].append([-100]*junction2[i])                

        else: #AND filter: All should pass the filter
            for junc_num, junction in enumerate(junc):
                all_pass = True
                if not (junction[i].sum() >= minreads): #if one of them doesn't pass the filter, flag and break
                    all_pass = False
                    break

            if all_pass:
                for junc_num, junction2 in enumerate(junc):
                    filtered_juncs[junc_num].append(junction2[i])   

            elif returnempty: 
                for junc_num, junction2 in enumerate(junc):
                    filtered_juncs[junc_num].append([-100]*len(junction2[i]))


    ret = [array(x) for x in filtered_juncs]
    filter_message("After discardminreads", minreads, logger, ret)
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


def discardlow(min0=0, orfilter=True, logger=False, *junc):
    filter_message("Before discardlow", min0, logger, junc)
    filtered_juncs = [[] for x in xrange(len(junc))]
    if orfilter:
        for i in range(junc[0].shape[0]): #for all junctions
            for junc_num, junction in enumerate(junc):
                if ((junction[i] > 0).sum() > min0):
                    #if one passes the filter, add all to resultset
                    for junc_num, junction2 in enumerate(junc):
                        filtered_juncs[junc_num].append(junction2[i])
                    break

        ret = [array(x) for x in filtered_juncs]

    else:    
        raise NotImplemented 

    filter_message("After discardlow", min0, logger, ret)

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


def mark_stacks(junctions, fitfunc, pvalue_limit, dispersion):
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
                my_nb = nbinom(r, p)
                pval = 1-my_nb.cdf(value)
                if pval < pvalue_limit:
                    junctions[i, j] = -2 
                    minstack = min(minstack, value)
                    numstacks += 1

    print "Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"%(junctions.size, numstacks, pvalue_limit, (float(numstacks)/junctions.size)*100)

    return junctions
