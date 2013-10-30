import sys
import os
import pickle
import argparse
from random import choice

from scipy.io import loadmat
from pylab import *
import numpy as np
from matplotlib import rcParams
from scipy.stats import pearsonr
from numpy.random import dirichlet
from numpy.ma import masked_less

from polyfit_ectf import norm_junctions
"""
We want to estimate the variability in sample1 and compare it, per junction with the empirical variability of sample2. The methods we want to compare are:

- Poisson (straw man) mean = variance
- Toronto bootstrapping
- Our improved bootstrapping using the fitted function and the NB

We want to plot the variance of method A vs method B in a scatter plot.

Same scatterplot for the winning method and the different normalizations.

"""
DEBUG = False
TESTBREAK = 800
LIM = 10
EPSILON = 1./sys.maxint
BORDER = 5 #definition of what a border is

def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname.replace(" ", "_")), bbox_inches='tight', width=1000, height=2000, dpi=300) #WNo spaces allowed, underscores!
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

def sample_from_junctions(junctions, m, k, discardzeros=False, nb=False, trimborder=False, parameters=None):
    dispersion = 0.1
    if parameters:
        print "Loading parameters from %s..."%parameters
        #poly_func = pickle.load(open('%sfitfunc.pickle'%args.par1))
        fitted_r = pickle.load(open('%snb_r.pickle'%parameters))
        fitted_p = pickle.load(open('%snb_p.pickle'%parameters))
        fitted_index = pickle.load(open('%snb_index.pickle'%parameters))
        fitted_func = pickle.load(open('%sfitfunc.pickle'%parameters))
        a, b = fitted_func.c

    sampled_mean = []
    sampled_var = []
    for i, junction in enumerate(junctions):
        junction = junction[junction > -EPSILON]  #discard the -1 (or lower) positions regardless of the dzero treatment
        if trimborder: 
            junction = _trimborders(junction) #trim the zeroes from the borders regardless of the discardzeros flag
        if discardzeros:
            junction = junction[junction!=0] #a junction array without the zeroes

        nb_samples = []
        if i % 100 == 0:
            print "junction %s"%i
            if DEBUG and i == TESTBREAK: break

        if len(junction) == 0:
            sampled_mean.append(0)
            sampled_var.append(0)
        else:
            if nb:
                emp_mean = mean(junction)
                #recalculating
                if (a*emp_mean)**2 > (emp_mean+dispersion*emp_mean**2):
                    r_nb = emp_mean**2 / ((a*emp_mean+b)**2 - emp_mean)
                else:
                    r_nb = 1/dispersion

                p_nb = r_nb / (r_nb+emp_mean)

            samples = []
            #sample m times
            for iternumber in xrange(m):
                junction_samples = []
                if nb:
                    junction_samples.extend(negative_binomial(r_nb, p_nb, k))
                else:
                    #for the regular sampling with replacement
                    #weights = calc_weights(junction)
                    #junction_samples = multinomial(k, weights, junction) #This function is very slow                    
                    for numsamples in xrange(k):
                        junction_samples.append(choice(junction))
                      
                samples.extend(junction_samples)

            #calculate the mean and the variance for simple sampling
            sampled_mean.append(mean(samples))
            sampled_var.append(var(samples))

    return array(sampled_mean), array(sampled_var)

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


def calc_psi(*samples):
    alpha = 0.5 
    numsamples = 100
    event_matrix = array(samples).reshape(-1, len(samples))
    psi_matrix = []
    for event in event_matrix:
        event_psi_samples = dirichlet(alpha+event, numsamples)
        #discretize the samples 
        print event_psi_samples

    print psi_matrix


def discardhigh(junctions1, junctions1_gc, junctions2, junctions2_gc, maxnonzero):
    numnonzeros1 = (junctions1 > 0).sum(axis=1)   
    pass_threshold = (numnonzeros1 < maxnonzero)
    junctions1 = junctions1[pass_threshold]
    junctions1_gc = junctions1_gc[pass_threshold]
    junctions2 = junctions2[pass_threshold]
    junctions2_gc = junctions2_gc[pass_threshold]
    #after and *not* before filtering the first replica, we calculate the second threshold
    numnonzeros2 = (junctions2 > 0).sum(axis=1)     
    pass_threshold = (numnonzeros2 < maxnonzero)
    return junctions1[pass_threshold], junctions1_gc[pass_threshold], junctions2[pass_threshold], junctions2_gc[pass_threshold] 


def discardlow(junctions1, junctions1_gc, junctions2, junctions2_gc, minnonzero):
    numnonzeros1 = (junctions1 > 0).sum(axis=1) 
    print numnonzeros1, minnonzero
    pass_threshold = (numnonzeros1 > minnonzero)
    junctions1 = junctions1[pass_threshold]
    junctions1_gc = junctions1_gc[pass_threshold]
    junctions2 = junctions2[pass_threshold]
    junctions2_gc = junctions2_gc[pass_threshold]
    print junctions1.shape, junctions2.shape
    #after and *not* before filtering the first replica, we calculate the second threshold
    numnonzeros2 = (junctions2 > 0).sum(axis=1)     
    pass_threshold = (numnonzeros2 > minnonzero)
    return junctions1[pass_threshold], junctions1_gc[pass_threshold], junctions2[pass_threshold], junctions2_gc[pass_threshold] 

def main():
    import matplotlib
    matplotlib = reload(matplotlib)
    matplotlib.use('GTKAgg')

    parser = argparse.ArgumentParser()
    parser.add_argument('matpath', help='Path with matfile with replica1 and replica2')    
    parser.add_argument('par1', help='Path for parameters of replica1')
    parser.add_argument('par2', help='Path for parameters of replica2')
    parser.add_argument('--sample', default=False, action='store_true', help='Use the negative binomial to sample')    
    parser.add_argument('--nb', default=False, action='store_true', help='Use the negative binomial to sample')
    parser.add_argument('--discardzeros', default=False, action='store_true', help='Discard the zeros from the junctions when computing means and variances')
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration')
    parser.add_argument('--m', default=100, type=int, help='Number of samples') 
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--junctype', default='rand10k', help='The type of junction to analyze. (Inc, Exc or rand10k for now)')
    parser.add_argument('--norm', default=False, action='store_true', help='Normalize by GC and weight factors')
    parser.add_argument('--poisson', default=False, action='store_true', help='Compare using the poisson distribution (Means=Variance)')   
    parser.add_argument('--trimborder', default=False, action='store_true', help='Trim the borders when sampling (keeping the ones with reads)')
    parser.add_argument('--minnonzero', default=0, type=int, help='Minimum number of positive positions to consider the junction')   
    parser.add_argument('--maxnonzero', default=0, type=int, help='Maximum number of positive positions to consider the junction') 
    parser.add_argument('--meanlim', default=25, type=int, help='Plot limit for the mean plotting (for comparison purposes)')
    parser.add_argument('--varlim', default=1500, type=int, help='Plot limit for the var plotting (for comparison purposes)')
    args = parser.parse_args()

    my_mat = loadmat(args.matpath)
    name1, name2 = os.path.basename(args.matpath).split('.')[0].split("_") #Get the experiment names from the mat file
    main_title = "%s VS %s\n"%(name1, name2) #main title of the plot generation

    replica1 = my_mat[name1][args.junctype][0, 0][0, 0]['cov']
    replica2 = my_mat[name2][args.junctype][0, 0][0, 0]['cov']
    replica1_gc = my_mat[name1][args.junctype][0, 0][0, 0]['gc_val']
    replica2_gc = my_mat[name2][args.junctype][0, 0][0, 0]['gc_val']


    if args.maxnonzero:
        print "Before Maxnonzero %s. replica1: %s replica2: %s"%(args.maxnonzero, replica1.shape, replica2.shape)
        replica1, replica1_gc, replica2, replica2_gc = discardhigh(replica1, replica1_gc, replica2, replica2_gc, args.maxnonzero)
        print "After Maxnonzero %s. replica1: %s replica2: %s"%(args.maxnonzero, replica1.shape, replica2.shape)      

    if args.minnonzero:
        print "Before Minnonzero %s. replica1: %s replica2: %s"%(args.minnonzero, replica1.shape, replica2.shape)
        replica1, replica1_gc, replica2, replica2_gc = discardlow(replica1, replica1_gc, replica2, replica2_gc, args.minnonzero)
        print "After Minnonzero %s. replica1: %s replica2: %s"%(args.minnonzero, replica1.shape, replica2.shape)      


    replica1 = masked_less(replica1, 0) 
    replica2 = masked_less(replica2, 0)


    if args.norm:
        replica1 = norm_junctions(replica1, gc_factors=replica1_gc, gcnorm=True)
        replica2 = norm_junctions(replica2, gc_factors=replica2_gc, gcnorm=True)

    if args.sample:
        main_title += " Sampling with repetition"
        my_mean1, my_var1 = sample_from_junctions(replica1, args.m, args.k, discardzeros=args.discardzeros, nb=args.nb, trimborder=args.trimborder, parameters=args.par1)
        my_mean2, my_var2 = sample_from_junctions(replica2, args.m, args.k, discardzeros=args.discardzeros, nb=args.nb, trimborder=args.trimborder, parameters=args.par2)

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
        replica1 = masked_less(replica1, 0) 
        replica2 = masked_less(replica2, 0) 
        if args.discardzeros:
            my_mean1, my_var1 = calc_nonzeromeanvar(replica1)
            my_mean2, my_var2 = calc_nonzeromeanvar(replica2)
        else:
            my_mean1 = replica1.mean(axis=1)
            my_var1 = replica1.var(axis=1)
            my_mean2 = replica2.mean(axis=1)
            my_var2 = replica2.var(axis=1) 

    rcParams.update({'font.size': 18})

    if args.discardzeros:
        main_title += " Discarding 0s"    



    fig = figure(figsize=[12, 20])

    suptitle(main_title, fontsize=24)
    subplot(2,1,1)
    plot_pearsoncorr(my_mean1, my_mean2, "Means", name1, name2, max_value=args.meanlim)
    subplot(2,1,2)
    plot_pearsoncorr(my_var1, my_var2, "Variances", name1, name2,  max_value=args.varlim)
    _save_or_show(args.plotpath, "")
    #calculate PSI
    print "PSI between %s and %s"%(name1, name2)
    calc_psi(my_mean1, my_mean2, my_var2)

if __name__ == '__main__':
    main()

