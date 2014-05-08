"""
We want to estimate the variability in sample1 and compare it, per junction with the empirical variability of sample2. The methods we want to compare are:

- Poisson (straw man) mean = variance
- Toronto bootstrapping
- Our improved bootstrapping using the fitted function and the NB

We want to plot the variance of method A vs method B in a scatter plot.

Same scatterplot for the winning method and the different normalizations.

"""
from matplotlib import use
use('Agg', warn=False)

import sys
import os
import pickle
import argparse
from random import choice

from pylab import *
import numpy as np
from matplotlib import rcParams
from scipy.stats import pearsonr
from numpy.random import dirichlet
from numpy.ma import masked_less
from polyfit_ectf import norm_junctions

import pipelines
import analysis.sample
import analysis.polyfitnb


DEBUG = True
TESTBREAK = 1500
LIM = 100
EPSILON = 1./sys.maxint
BORDER = 5 #definition of what a border is

def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)
        savefig("%s/%s.png" % (plotpath, plotname.replace(" ", "_")), bbox_inches='tight', width=1000, height=2000, dpi=300) # No spaces allowed, underscores!
        clf()
    else:
        show()


def calc_nonzeromeanvar(junctions):
    nonzeromean = []
    nonzerovar = []
    for junction in junctions:
        nonzerojunction = junction[junction<=0] #discard also -1 and -2, just in case
        if len(nonzerojunction) > 0:
            nonzeromean.append(mean(nonzerojunction))
            nonzerovar.append(var(nonzerojunction))
        else:
            nonzeromean.append(0)
            nonzerovar.append(0)

    return array(nonzeromean), array(nonzerovar)

def calc_weights(junction):
    #Not used yet. The weight_factors array is wrong
    weight_factors =  array([0.0, 0.67221794713393235, 0.822110928783363, 1.0, 0.96555466736662077, 0.97034079788351413, 0.95888448040949026, 0.93916171702266071,
                    0.9586120159245054, 0.98886869384130682, 0.97031850460351132, 0.97919730108521363, 0.98845338964933505, 0.97762016912650862,
                    0.96000443463677787, 0.97795543996841916, 0.98483897464546255, 0.98430701054559211, 0.97621404631851538, 0.97557091482162561,
                    0.99783624218670419, 0.99420256804654417, 0.99996092005852, 1.0, 0.98891003681022327, 0.98408910298925234, 0.98588911669260937,
                    0.98944197552348001, 0.98861266787997559, 0.98334059809099128, 0.98616818121835459, 0.98356568445706294, 1.0, 1.0,
                    0.99415876734414588, 1.0, 0.99732413178991319, 0.98571657557568526, 0.98637294512249951, 0.98846297241187242, 1.0,
                    0.98857076368303576, 1.0, 0.98474007029306976, 0.98212050612598556, 0.99227062085183826, 0.98716235724225032, 0.98604617629365343,
                    0.96908030440229109, 0.97105918006649872, 0.97297718733803484, 0.98431591864639367, 0.98227616224387038, 0.9961571944449884,
                    0.97565056267585271, 0.96725772937340826, 0.95469906291036666, 0.94761567083759468, 0.96284719373281014, 1.0, 0])

    #TODO correction for the borders
    return weight_factors/len(weight_factors)


def _trimborders(junction):
    "Discard the borders of the junctions unless they have reads"
    #discard left side
    new_junction = []
    for i in range(0, BORDER):
        if junction[i] > 0:
            new_junction.append(junction[i])

    new_junction.extend(junction[BORDER:-BORDER]) #add the middle positions
    #discard right side
    for i in range(len(junction)-1, len(junction)-BORDER-1, -1):
        if junction[i] > 0:
            new_junction.append(junction[i])

    return array(new_junction)


def sample_from_junctions(junctions, m, k, discardzeros=False, nb=False, trimborder=False, parameters=None, dispersion=0.1, fit_func=None, poisson=False):

    if nb:
        return analysis.sample.sample_from_junctions(junctions, m, k, discardzeros=discardzeros, dispersion=dispersion, trimborder=trimborder, fitted_func=fit_func)

    # if parameters:
    #     print "Loading parameters from %s..." % parameters
    #     fitted_func = pickle.load(open('%sfitfunc.pickle' % parameters))
    #     a, b = fitted_func.c

    sampled_means = []
    sampled_var = []
    all_samples = []
    for i, junction in enumerate(junctions):
        # if i % 100 == 0:
        #     print "junction %s..."%i,
        #     sys.stdout.flush()

        junction = junction[junction > -EPSILON]  #discard the -1 (or lower) positions regardless of the dzero treatment

        if trimborder:
            junction = analysis.sample._trimborders(junction) #trim the zeroes from the borders regardless of the discardzeros flag
        if discardzeros:
            junction = junction[junction!=0] #a junction array without the zeroes

        if (poisson):
            mean_poisson = np.mean(junction, axis=0)
            sampled_means.append(mean_poisson)
            sampled_var.append(mean_poisson)
            all_samples.append(junction)
        else:
            if len(junction) == 0:
                sampled_means.append(0)
                sampled_var.append(0)
                all_samples.append([0]*(k*m)) #k*m zeroes
            else:
                samples = []
                for iternumber in xrange(m):
                    junction_samples = []
                    #using the matrix
                    #weights = calc_weights(junction)
                    #junction_samples = multinomial(k, weights, junction) #This function is very slow
                    for numsamples in xrange(k):
                        junction_samples.append(choice(junction))

                    samples.extend(junction_samples)

                #calculate the mean and the variance
                sampled_means.append(mean(samples))
                sampled_var.append(var(samples))
                all_samples.append(samples)

    return array(sampled_means), array(sampled_var), array(all_samples)

def plot_pearsoncorr(var1, var2, my_title, my_xlabel, my_ylabel, plotpath=None, max_value=None):
    if DEBUG:
        var1 = var1[:TESTBREAK]
        var2 = var2[:TESTBREAK]

    var1 = array(var1)
    var2 = array(var2)
    xlabel(my_xlabel)
    ylabel(my_ylabel)

    if not max_value:
        max_value = max(max(var1), max(var2))

    xlim(0, max_value)
    ylim(0, max_value)

    #plot([0, max_value], [0, max_value])
    print "PEARSON",var1, var2
    pear, pvalue = pearsonr(var1, var2)
    r_squared = pear**2

    a, b = polyfit(var1, var2, 1)
    fit_func = poly1d([a, b])
    plot(var1, fit_func(var1), '-r')
    #percentage_under = sum(var1 < var2)/float(len(var1)) 
    text(abs(max_value)*0.1, max_value-abs(max_value)*0.2, r'$R^2$: %.2f (p-value: %.2E)'%(r_squared, pvalue), fontsize=18, bbox={'facecolor':'yellow', 'alpha':0.3, 'pad':10})
    title(my_title)
    print r"%s R^2: %.2f (p-value: %.2E)"%(my_title, r_squared, pvalue)
    plot(var1, var2, '.')
    if plotpath:
        _save_or_show(plotpath, my_title)

def filter_bulk(matrix_filter, *matrices):
    ret = []
    for m in matrices:
        ret.append(m[matrix_filter])

    return ret


def check_junctions_in_replicates(lsv_junc1, lsv_junc2, discard_empty_junctions=False):
    ids1 = set([x[1] for x in lsv_junc1[1]])
    ids2 = set([x[1] for x in lsv_junc2[1]])

    matched_names = ids1.intersection(ids2)
    print len(ids1), len(ids2)
    print len(matched_names)
    replica1 = []
    replica2 = []
    for ii in matched_names:
        for idx, nm in enumerate(lsv_junc1[1]):
            if nm[1] == ii:
                replica1.append(lsv_junc1[0][idx])
                dummy=lsv_junc1[0][idx].shape
                break
        for idx, nm in enumerate(lsv_junc2[1]):
            if nm[1] == ii:
                dummy2 = lsv_junc2[0][idx].shape
                if dummy != dummy2:
                    print "ERRRRORRRRRR", dummy, dummy2
                    replica1 = replica1[:-1]
                else:
                    replica2.append(lsv_junc2[0][idx])

                break

    replica1 = np.concatenate(replica1)
    replica1 = replica1.astype(np.float64)
    replica2 = np.concatenate(replica2)
    replica2 = replica2.astype(np.float64)

    if discard_empty_junctions:
        idx_list = []
        for idx in range(replica1.shape[0]):
            if np.count_nonzero(replica1[idx]) == 0 or np.count_nonzero(replica2[idx]) ==0 : idx_list.append(idx)
        replica1 = np.delete(replica1, idx_list, axis=0)
        replica2 = np.delete(replica2, idx_list, axis=0)

    return replica1, replica2


def discard_emtpy_junctions( replica ):
    idx_list = []
    for idx in range(replica.shape[0]):
        if np.count_nonzero(replica[idx]) == 0 : idx_list.append(idx)
    replica = np.delete(replica, idx_list, axis=0)

    return replica

def load_junctions(filename1, filename2, args, fromlsv=False):

    # Parse LSV files
    lsv_junc1, const1 = pipelines.load_data_lsv(filename1)
    lsv_junc2, const2 = pipelines.load_data_lsv(filename2)

    fit_func1 = polyfitnb.fit_nb(const1, "%s_nbfit" % args.output, args.plotpath, nbdisp=args.dispersion, logger=None, discardb=True)
    fit_func2 = polyfitnb.fit_nb(const2, "%s_nbfit" % args.output, args.plotpath, nbdisp=args.dispersion, logger=None, discardb=True)

    if fromlsv:
        replica1, replica2 = junction_sample.check_junctions_in_replicates(lsv_junc1, lsv_junc2, discard_empty_junctions=True)
    else:
        replica1 = discard_empty_junctions(const1)
        replica2 = discard_empty_junctions(const2)

    return replica1, replica2, fit_func1, fit_func2


def split_junction_pool ( replica1, replica2 ):

    low = []
    mid = []
    high = []

    for idx in range(replica1.shape[0]):
        if np.count_nonzero(replica1[idx]) in range(1,6) :
            low.append(replica1[idx])


def main():
    """Script for initial testing of the MAJIQ algorithms for sampling and initial PSI values generator."""
    # TODO: Many functions from this script will be extracted for general usage in the pipeline.

    samples1 = None
    samples2 = None

    parser = argparse.ArgumentParser()
    parser.add_argument('rep1', help='Path for replica1')
    parser.add_argument('rep2', help='Path for replica2')
    parser.add_argument('--sample', default=False, action='store_true', help='Use sampling')
    parser.add_argument('--nb', default=False, action='store_true', help='Use the negative binomial to sample')
    parser.add_argument('--discardzeros', default=False, action='store_true', help='Discard the zeros from the junctions when computing means and variances')
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration')
    parser.add_argument('--m', default=100, type=int, help='Number of bootstrapping samples')
    parser.add_argument('--alpha', default=0.5, type=int, help='Alpha hyperparameter for the dirichlet distribution')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--junctype', default='rand10k', help='The type of junction to analyze. (Inc, Exc or rand10k for now)')
    parser.add_argument('--norm', default=False, action='store_true', help='Normalize by GC and weight factors')
    parser.add_argument('--poisson', default=False, action='store_true', help='Compare using the poisson distribution (Means=Variance)')
    parser.add_argument('--trimborder', default=False, action='store_true', help='Trim the borders when sampling (keeping the ones with reads)')
    parser.add_argument('--minnonzero', default=0, type=int, help='Minimum number of positive positions to consider the junction')
    parser.add_argument('--maxnonzero', default=sys.maxint, type=int, help='Maximum number of positive positions to consider the junction')
    parser.add_argument('--minreads', default=0, type=int, help='Minimum number of reads combining all positions in a junction to be considered')
    parser.add_argument('--maxreads', default=sys.maxint, type=int, help='Maximum number of reads combining all positions in a junction to be considered')
    parser.add_argument('--meanlim', default=25, type=int, help='Plot limit for the mean plotting (for comparison purposes)')
    parser.add_argument('--varlim', default=1500, type=int, help='Plot limit for the var plotting (for comparison purposes)')
    parser.add_argument('--orfilter', default=False, action='store_true', help='When filtering, select sets of junctions where at least one passes the filter, instead of all passing the filter.')
    parser.add_argument('--output', default=None, help="Path to save the results to.")
    parser.add_argument('--dispersion', default=0.1, type=float, help='Dispersion factor (used in junctions sampling).')
    args = parser.parse_args()

    # Parse LSV files
    lsv_junc1, const1 = pipelines.load_data_lsv(args.rep1)
    lsv_junc2, const2 = pipelines.load_data_lsv(args.rep2)

    # Check that all LSVs have the same num. of junctions in both replicates
    replica1, replica2 = check_junctions_in_replicates(lsv_junc1, lsv_junc2)

    # Collapse all constitutive junctions for sampling models
    all_const = np.concatenate([const1, const2])

    #Get the experiment names
    name1 = os.path.basename(args.rep1).split('.')[-2]
    name2 = os.path.basename(args.rep2).split('.')[-2]
    main_title = "%s VS %s\n"%(name1, name2)

    # TODO -ing: divide script in two sections: normalization and sampling

    if args.sample:
        main_title += " Sampling with repetition"
        my_mean1, my_var1, samples1 = sample_from_junctions(replica1, args.m, args.k, discardzeros=args.discardzeros, dispersion=args.dispersion, trimborder=args.trimborder, fitted_func=fit_func, debug=DEBUG)
        my_mean2, my_var2, samples2 = sample_from_junctions(replica2, args.m, args.k, discardzeros=args.discardzeros, dispersion=args.dispersion, trimborder=args.trimborder, fitted_func=fit_func, debug=DEBUG)
        # my_mean1, my_var1, samples1 = sample_from_junctions(replica1, args.m, args.k, discardzeros=args.discardzeros, nb=args.nb, trimborder=args.trimborder)
        # my_mean2, my_var2, samples2 = sample_from_junctions(replica2, args.m, args.k, discardzeros=args.discardzeros, nb=args.nb, trimborder=args.trimborder)


        if args.poisson:
            my_var1 = my_mean1
            print "Made Variance of %s equal to Mean"%name1
            main_title += " (Poisson)"
        elif args.nb:
            main_title += " (Negative Binomial)"
        else:
            main_title += " (Naive Boostrapping)"

        main_title += '\n'

    else:
        main_title += " Empirical\n"
#        replica1 = masked_less(replica1, 0)
#        replica2 = masked_less(replica2, 0)
        print replica1.shape
        if args.discardzeros:
            my_mean1, my_var1 = calc_nonzeromeanvar(replica1)
            my_mean2, my_var2 = calc_nonzeromeanvar(replica2)
        else:
            my_mean1 = np.mean(replica1, axis=1)
            my_var1 =  np.var(replica1, axis=1)
            my_mean2 = np.mean(replica2, axis=1)
            my_var2 =  np.var(replica2, axis=1)


    #plot the means and the variances
    rcParams.update({'font.size': 18})

    if args.discardzeros:
        main_title += " Discarding 0s"

    if args.trimborder:
        main_title += " Trimming ends"

    fig = figure(figsize=[12, 20])

    suptitle(main_title, fontsize=24)
    subplot(2, 1, 1)
    plot_pearsoncorr(my_mean1, my_mean2, "Means", name1, name2, max_value=args.meanlim)
    subplot(2, 1, 2)
    plot_pearsoncorr(my_var1, my_var2, "Variances", name1, name2,  max_value=args.varlim)
    _save_or_show(args.plotpath, "meanvarcorr")

    #calculate PSI
    if args.output:
        if args.sample or args.miso:
            print "Saving samples..."
            pickle.dump(samples1, open("%s_%s_samples.pickle" % (args.output, name1), 'w'))
            pickle.dump(samples2, open("%s_%s_samples.pickle" % (args.output, name2), 'w'))
            print "Done!"



if __name__ == '__main__':
    main()


## LEGACY CODE from Juan, DEPRECATED

#    replica1_gc = my_mat[name1][args.junctype][0, 0][0, 0]['gc_val']
#    replica2_gc = my_mat[name2][args.junctype][0, 0][0, 0]['gc_val']

#   #Coverage filters
#    if args.maxnonzero:
#        print "Before Maxnonzero %s. replica1: %s replica2: %s"%(args.maxnonzero, replica1.shape, replica2.shape)
#        replica1, replica1_gc, replica2, replica2_gc = discardhigh(replica1, replica1_gc, replica2, replica2_gc, args.maxnonzero, args.orfilter)
#        print "After Maxnonzero %s. replica1: %s replica2: %s"%(args.maxnonzero, replica1.shape, replica2.shape)
#
#    if args.minnonzero:
#        print "Before Minnonzero %s. replica1: %s replica2: %s"%(args.minnonzero, replica1.shape, replica2.shape)
#        replica1, replica1_gc, replica2, replica2_gc = discardlow(replica1, replica1_gc, replica2, replica2_gc, args.minnonzero, args.orfilter)
#        print "After Minnonzero %s. replica1: %s replica2: %s"%(args.minnonzero, replica1.shape, replica2.shape)
#
#    if args.minreads:
#        print "Before Minreads %s. replica1: %s replica2: %s"%(args.minreads, replica1.shape, replica2.shape)
#        replica1, replica1_gc, replica2, replica2_gc = discardminreads(replica1, replica1_gc, replica2, replica2_gc, args.minreads, args.orfilter)
#        print "After Minreads %s. replica1: %s replica2: %s"%(args.minreads, replica1.shape, replica2.shape)
#
#    if args.maxreads:
#        print "Before maxreads %s. replica1: %s replica2: %s"%(args.maxreads, replica1.shape, replica2.shape)
#        replica1, replica1_gc, replica2, replica2_gc = discardmaxreads(replica1, replica1_gc, replica2, replica2_gc, args.maxreads, args.orfilter)
#        print "After maxreads %s. replica1: %s replica2: %s"%(args.maxreads, replica1.shape, replica2.shape)


#    #mask everything under 0
#    if args.norm:
#        replica1 = norm_junctions(replica1, gc_factors=replica1_gc, gcnorm=True)
#        replica2 = norm_junctions(replica2, gc_factors=replica2_gc, gcnorm=True)


# def discardhigh(j1, j1_gc, j2, j2_gc, max0, orfilter=False):
#     if orfilter:
#         j1_filtered = []
#         j1_gc_filtered = []
#         j2_filtered = []
#         j2_gc_filtered = []
#         for i in range(j1.shape[0]):
#             if ((j1[i] > 0).sum() < max0) or ((j2[i] > 0).sum() < max0):
#                 j1_filtered.append(j1[i])
#                 j1_gc_filtered.append(j1_gc[i])
#                 j2_filtered.append(j2[i])
#                 j2_gc_filtered.append(j2_gc[i])
#
#         return array(j1_filtered), array(j1_gc_filtered), array(j2_filtered), array(j2_gc_filtered)
#
#     else:
#         #NOTE: For this kind of calculation, numnonzeros2 *MUST* be calculated after the first filter, or else it will not fit the new matrix. Unless you find a way of make combined filters...
#         numnonzeros1 = (j1 > 0).sum(axis=1)
#         j1, j1_gc, j2, j2_gc = filter_bulk((numnonzeros1 < max0), j1, j1_gc, j2, j2_gc)
#         numnonzeros2 = (j2 > 0).sum(axis=1)
#         return filter_bulk((numnonzeros2 < max0), j1, j1_gc, j2, j2_gc)
#
#
# def discardminreads(j1, j1_gc, j2, j2_gc, minreads, orfilter=False):
#     if orfilter:
#         j1_filtered = []
#         j1_gc_filtered = []
#         j2_filtered = []
#         j2_gc_filtered = []
#         for i in range(j1.shape[0]):
#             if (j1[i].sum() >= minreads) or (j2[i].sum() >= minreads):
#                 j1_filtered.append(j1[i])
#                 j1_gc_filtered.append(j1_gc[i])
#                 j2_filtered.append(j2[i])
#                 j2_gc_filtered.append(j2_gc[i])
#
#         return array(j1_filtered), array(j1_gc_filtered), array(j2_filtered), array(j2_gc_filtered)
#     else:
#         readsj1 = j1.sum(axis=1)
#         j1, j1_gc, j2, j2_gc = filter_bulk((readsj1 >= minreads), j1, j1_gc, j2, j2_gc)
#         readsj2 = j2.sum(axis=1)
#         return filter_bulk((readsj2 >= minreads), j1, j1_gc, j2, j2_gc)
#
# def discardmaxreads(j1, j1_gc, j2, j2_gc, maxreads, orfilter=False):
#     if orfilter:
#         j1_filtered = []
#         j1_gc_filtered = []
#         j2_filtered = []
#         j2_gc_filtered = []
#         for i in range(j1.shape[0]):
#             if (j1[i].sum() <= maxreads) or (j2[i].sum() <= maxreads):
#                 j1_filtered.append(j1[i])
#                 j1_gc_filtered.append(j1_gc[i])
#                 j2_filtered.append(j2[i])
#                 j2_gc_filtered.append(j2_gc[i])
#
#         return array(j1_filtered), array(j1_gc_filtered), array(j2_filtered), array(j2_gc_filtered)
#     else:
#         readsj1 = j1.sum(axis=1)
#         j1, j1_gc, j2, j2_gc = filter_bulk((readsj1 <= maxreads), j1, j1_gc, j2, j2_gc)
#         readsj2 = j2.sum(axis=1)
#         return filter_bulk((readsj2 <= maxreads), j1, j1_gc, j2, j2_gc)
#
#
# def discardlow(j1, j1_gc, j2, j2_gc, min0, orfilter=False):
#     if orfilter:
#         j1_filtered = []
#         j1_gc_filtered = []
#         j2_filtered = []
#         j2_gc_filtered = []
#         for i in range(j1.shape[0]):
#             if ((j1[i] > 0).sum() > min0) or ((j2[i] > 0).sum() > min0):
#                 j1_filtered.append(j1[i])
#                 j1_gc_filtered.append(j1_gc[i])
#                 j2_filtered.append(j2[i])
#                 j2_gc_filtered.append(j2_gc[i])
#
#         return array(j1_filtered), array(j1_gc_filtered), array(j2_filtered), array(j2_gc_filtered)
#
#     else:
#         numnonzeros1 = (j1 > 0).sum(axis=1)
#         j1, j1_gc, j2, j2_gc = filter_bulk((numnonzeros1 > min0), j1, j1_gc, j2, j2_gc)
#         numnonzeros2 = (j2 > 0).sum(axis=1)
#         return filter_bulk((numnonzeros2 > min0), j1, j1_gc, j2, j2_gc)

