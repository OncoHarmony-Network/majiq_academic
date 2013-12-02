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


def filter_message(when, value, debug, junc):
    if debug: 
        print "%s (Filter=%s). %s"%(when, value, junc[0].shape)


def discardhigh(max0=0, orfilter=False, debug=False, *junc):

    filter_message("Before discardhigh", max0, debug, junc)

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
        #NOTE: For this kind of calculation, numnonzeros2 *MUST* be calculated after the first filter, or else it will not fit the new matrix. Unless you find a way of make combined filters...
        numnonzeros1 = (j1 > 0).sum(axis=1)
        j1, j2= filter_bulk((numnonzeros1 < max0), j1,j2)
        numnonzeros2 = (j2 > 0).sum(axis=1)
        ret = filter_bulk((numnonzeros2 < max0), j1, j2)

    filter_message("After discardhigh", max0, debug, junc)

    return ret


def discardminreads(minreads=0, orfilter=False, debug=False, *junc):

    filter_message("Before discardminreads", minreads, debug, junc)

    filtered_juncs = [[] for x in xrange(len(junc))]
    if orfilter:
        for i in range(junc[0].shape[0]): #for all junctions
            for junc_num, junction in enumerate(junc):
                if (junction[i].sum() >= minreads):
                    #if one passes the filter, add all to resultset
                    for junc_num, junction2 in enumerate(junc):
                        filtered_juncs[junc_num].append(junction2[i])
                    break

        ret = [array(x) for x in filtered_juncs]

    else:
        raise NotImplemented
        """
        readsj1 = j1.sum(axis=1)
        j1, j2 = filter_bulk((readsj1 >= minreads), j1, j2)
        readsj2 = j2.sum(axis=1)
        ret = filter_bulk((readsj2 >= minreads), j1, j2)
        """
    filter_message("After discardminreads", minreads, debug, ret)
    
    return ret


def discardmaxreads(maxreads=0, orfilter=False, debug=False, *junc):

    filter_message("Before discardmaxreads", maxreads, debug, junc)

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
        readsj1 = j1.sum(axis=1)
        j1, j2 = filter_bulk((readsj1 <= maxreads), j1, j2)
        readsj2 = j2.sum(axis=1)
        ret = filter_bulk((readsj2 <= maxreads), j1, j2)

    filter_message("After discardmaxreads", maxreads, debug, ret)

    return ret


def discardlow(min0=0, orfilter=False, debug=False, *junc):
    filter_message("Before discardlow", min0, debug, junc)
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
        """  
        numnonzeros1 = (j1 > 0).sum(axis=1)
        j1, j2 = filter_bulk((numnonzeros1 > min0), j1, j2)
        numnonzeros2 = (j2 > 0).sum(axis=1)
        ret = filter_bulk((numnonzeros2 > min0), j1, j2)
        """
    filter_message("After discardlow", min0, debug, ret)

    return ret


def norm_junctions(junctions, gc_factors=None, gcnorm=False, trim=False, debug=False):
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
                r, p = func2nb(a, b, value, dispersion)
                my_nb = nbinom(r, p)
                pval = 1-my_nb.cdf(value)
                if pval < pvalue_limit:
                    junctions[i, j] = -2
                    minstack = min(minstack, value)
                    numstacks += 1

    print "Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"%(junctions.size, numstacks, pvalue_limit, (float(numstacks)/junctions.size)*100)

    return junctions
