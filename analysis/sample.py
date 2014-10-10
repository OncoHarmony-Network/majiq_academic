import sys
from random import choice
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
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
        plt.savefig("%s%s.png"%(plotpath, plotname.replace(" ", "_")), bbox_inches='tight', width=1000, height=2000,
                    dpi=300) #WNo spaces allowed, underscores!
        plt.clf()
    else:
        plt.show()


def calc_nonzeromeanvar(junctions):
    nonzeromean = []
    nonzerovar = []
    for junction in junctions:
        nonzerojunction = junction[junction <= 0]
        #discard also -1 and -2, just in case
        if len(nonzerojunction) > 0:
            nonzeromean.append(np.mean(nonzerojunction))
            nonzerovar.append(np.var(nonzerojunction))
        else:
            nonzeromean.append(0)
            nonzerovar.append(0)

    return np.array(nonzeromean), np.array(nonzerovar)


def _trimborders(junction, border):
    "Discard the borders of the junctions unless they have reads"
    #TODO use a masker strategy instead of generating a new array (warning: Assigning values in a numpy
    # array creates a new one, so it is more inneficient than this)
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

    return np.array(new_junction)


def remove_masked(junction):
    "For performance: Less values to sample from, faster execution time"
    ret = []
    for value in junction:
        if value > -EPSILON:
        #zero and bigger than zero
            ret.append(value)

    return np.array(ret)


def mean_junction(junctions, discardzeros=True):
    """Simple mean of junctions without bootstrapping, but discarding zeroes and flagged stuff"""
    ret = []
    for junc in junctions:
        junc = junc[junc > -EPSILON]
        #mask the -1 (or lower) positions regardless of the discardzero treatment
        if discardzeros:
            junc = junc[junc!=0]
             #a junc array without the zeroes

        if len(junc) == 0:
            ret.append(0) 
        else:
            ret.append(junc.mean())

    return np.array(ret)

global EMPIRICAL_NZ, BINOMIAL_NZ 
EMPIRICAL_NZ = 0
BINOMIAL_NZ = -1


def sample_from_junctions(junction_list, m, k, dispersion=0.1, discardzeros=5, trimborder=True, fitted_func=None,
                          debug=False, tracklist=None, Nz=0, names=None, naive=False):
    "Given the filtered reads, bootstrap samples from every junction"
    a, b = fitted_func.c
    sampled_means = []
    sampled_var = []
    all_samples = []
    
    for i, junction in enumerate(junction_list):
        if 0 < debug == i:
            break
        if i % 100 == 0 and debug > 0:
            print "junction %s..." % i,
            sys.stdout.flush()

        if trimborder: 
            junction = _trimborders(junction, trimborder)
            #trim the zeroes from the borders regardless of the discardzeros flag

        junction = junction[junction > -EPSILON]
        #mask the -1 (or lower) positions regardless of the discardzero treatment
        
        if discardzeros > 0:
            junction = junction[junction != 0]
            #a junction array without the zeroes
            sys.stdout.flush()
            if junction.shape[0] < discardzeros:
                z = np.zeros(shape=(discardzeros-junction.shape[0]), dtype=int)
                junction = np.concatenate((junction, z))
                #a junction array without the zeroes

        if np.count_nonzero(junction) == 0:
            sampled_means.append(0)
            sampled_var.append(0)
            all_samples.append([0]*m)
            #k*m zeroes
        else:

            if Nz == EMPIRICAL_NZ:
                npos_mult = np.count_nonzero(junction)
            elif Nz == BINOMIAL_NZ:
                npos_mult = np.binomial(k, float(np.count_nonzero(junction))/float(k), k*m)
            else:
                npos_mult = Nz

            samples = []
            for iternumber in xrange(m):
                junction_samples = []
                for numsamples in xrange(k):
                    junction_samples.append(choice(junction))

                sampled_mean = np.mean(junction_samples)
                #recalculating
                r_nb, p_nb = polyfitnb.func2nb(a, b, sampled_mean, dispersion)
                nb50 = np.random.negative_binomial(r_nb, p_nb, k)
                smpl = np.mean(nb50)
                samples.append(smpl)
            #calculate the mean and the variance 

            sampled_means.append(np.mean(samples))
            sampled_var.append(np.var(samples))
           
#            samples = [ npos_mult* (x+1) for x in samples]
            samples = npos_mult * (np.array(samples)+1)
#            print i, Nz, samples
            all_samples.append(samples)

    return np.array(sampled_means), np.array(sampled_var), np.array(all_samples)


def plot_pearsoncorr(var1, var2, my_title, my_xlabel, my_ylabel, plotpath=None, max_value=None):

    var1 = np.array(var1)
    var2 = np.array(var2)
    plt.xlabel(my_xlabel)
    plt.ylabel(my_ylabel)

    if not max_value:
        max_value = max(max(var1), max(var2))

    plt.xlim(0, max_value)
    plt.ylim(0, max_value)
    
    #plot([0, max_value], [0, max_value])
    pear, pvalue = pearsonr(var1, var2)
    r_squared = pear**2
    a, b = np.polyfit(var1, var2, 1)
    fit_func = np.poly1d([a, b])
    plt.plot(var1, fit_func(var1), '-r')
    #percentage_under = sum(var1 < var2)/float(len(var1)) 
    plt.text(abs(max_value) * 0.1, max_value-abs(max_value)*0.2, r'$R^2$: %.2f (p-value: %.2E)' % (r_squared, pvalue),
             fontsize=18, bbox={'facecolor': 'yellow', 'alpha': 0.3, 'pad': 10})
    plt.title(my_title)
    print r"%s R^2: %.2f (p-value: %.2E)" % (my_title, r_squared, pvalue)
    plt.plot(var1, var2, '.')
    if plotpath:
        _save_or_show(plotpath, my_title)