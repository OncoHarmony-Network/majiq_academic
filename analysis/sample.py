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

import analysis.polyfitnb as polyfitnb

"""
Sampling from junctions using a Negative Binomial model.
"""
LIM = 100
EPSILON = 1./sys.maxint
PSEUDO = 0.0000000001 # EPSILON is too small for some calculations

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


def _trimborders(junction, border):
    "Discard the borders of the junctions unless they have reads"
    #TODO use a masker strategy instead of generating a new array (warning: Assigning values in a numpy array creates a new one, so it is more inneficient than this)
    #discard left side
    new_junction = []
    for i in range(0, border):
        if junction[i] > 0:
            new_junction.append(junction[i])

    new_junction.extend(junction[border:-border]) #add the middle positions, disregard masked positions
    #discard right side
    for i in range(len(junction)-1, len(junction)-border-1, -1):
        if junction[i] > 0:
            new_junction.append(junction[i])

    return array(new_junction)

def remove_masked(junction):
    "For performance: Less values to sample from, faster execution time"
    ret = []
    for value in junction:
        if value > -EPSILON: #zero and bigger than zero
            ret.append(value)

    return array(ret)


def mean_junction(junctions, discardzeros=True):
    """Simple mean of junctions without bootstrapping, but discarding zeroes and flagged stuff"""
    ret = []
    for junc in junctions:
        junc = junc[junc > -EPSILON]  #mask the -1 (or lower) positions regardless of the discardzero treatment
        if discardzeros:
            junc = junc[junc!=0] #a junc array without the zeroes

        if len(junc) == 0:
            ret.append(0) 
        else:
            ret.append(junc.mean())

    return array(ret)

def sample_from_junctions(junctions, m, k, dispersion=0.1, discardzeros=5, trimborder=True, fitted_func=None, debug=False, tracklist=None, names=None):
    "Given the filtered reads, bootstrap samples from every junction"
    a, b = fitted_func.c
    sampled_means = []
    sampled_var = []
    all_samples = []
    for i, junction in enumerate(junctions):
        if debug > 0 and i == debug: break
        if i % 100 == 0:
            print "junction %s..."%i,
            sys.stdout.flush()

        if trimborder: 
            junction = _trimborders(junction, trimborder) #trim the zeroes from the borders regardless of the discardzeros flag

        junction = junction[junction > -EPSILON]  #mask the -1 (or lower) positions regardless of the discardzero treatment

        if discardzeros > 0:
            junction = junction[junction!=0] #a junction array without the zeroes
            sys.stdout.flush()

            if junction.shape[0]< discardzeros:
                z = np.zeros(shape=(discardzeros-junction.shape[0]), dtype=int)
                junction = np.concatenate((junction,z), axis=1) #a junction array without the zeroes

        if len(junction) == 0:
            sampled_means.append(0)
            sampled_var.append(0)
            all_samples.append([0]*(k*m)) #k*m zeroes
        else:
            samples = []
            for iternumber in xrange(m):
                junction_samples = []
                for numsamples in xrange(k):
                    junction_samples.append(choice(junction))

                sampled_mean = mean(junction_samples)
                #recalculating
                r_nb, p_nb = polyfitnb.func2nb( a, b, sampled_mean, dispersion )
#                if (a*sampled_mean)**2 > (sampled_mean+dispersion*sampled_mean**2):
#                    r_nb = sampled_mean**2 / ((a*sampled_mean+b)**2 - sampled_mean)
#                else:
#                    r_nb = 1/dispersion
#
#                p_nb = r_nb / (r_nb+sampled_mean)
                samples.extend(negative_binomial(r_nb, p_nb, k)) 
            
            #calculate the mean and the variance 
            sampled_means.append(mean(samples))
            sampled_var.append(var(samples))
            all_samples.append(samples)
            if names and tracklist:
                if names[i] in tracklist:
                    logger.info("TRACKLIST (%s): %s"%(event_name, junctions))

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



if __name__ == '__main__':
    main()

