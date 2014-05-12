import sys
import os
import argparse
import pickle
from collections import defaultdict

from numpy.ma import masked_less
from scipy.io import loadmat
from scipy.stats import nbinom
from pylab import *


#TODO: every function that translates a, b to r, p should use this function if needed for single values, or the one below for a whole ndarray
def func2nb(a, b, x, dispersion):
    """
    Given a and b from the linear fit, calculate the nb parameters for x and return its r and p parameters.
    """
#    if (a*x**2 > (x+dispersion*x**2) or b == 1 :
    if (a*x**2 + b*x> x ) or b == 1 :
#        r = x**2 / ((a*x+b)**2 - x)
        r = x / ( a*x +b - 1 )
    else:
        r = 1/dispersion

    p =  (r / (r+x))
#    p = 1 - (r / (r+x))
    return r, p


def nb_from_func(poly_func, max_line=1, dispersion=0.1):
    """
    From a linear function that describes the distribution of junctions mean/variance, get the NB values

    We use a combined approach to avoid having values of r < 0. 
    For small values, we fit the r *global* negative binomial with a dispersion parameter (1/r in the NB formula)
    of 0.01, which is a bit bigger variance than a simple Poisson. 

    """
    a, b = poly_func.c
    r = []
    p = []

    points = linspace(0.01, max_line, num=800)
    for x in points:
        r_val, p_val = func2nb(a, b, x, dispersion)
        r.append( r_val )
        p.append( p_val )
#        p.append(r_val / (r_val+x))

    return r, p, points

def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if a flag was given"""
    if plotname:
        title(plotname)

    if plotpath:
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)
        savefig("%s/%s.png"%(plotpath, plotname.replace(" ", "_")), bbox_inches='tight') #WNo spaces allowed, underscores!
        clf()
    else:
        show()

def get_ecdf(pvalues):
    hist, bin_edges = histogram(pvalues, range=[0,1], bins=len(pvalues)/10, density=True)
    return cumsum(hist)/len(bin_edges)

def plot_mappability_zeros(junctions, plotpath, numzeros, plotname):
    """
    Plot the percentage of zero positions in the junction sites against the total sum of the positions in the junction. 
    This should be proportional. Dots in the right mean that there are some zero positions in the junctions with high coverage, meaning that 
    there is probably some mappability bias. 
    """
    junction_sum = junctions.sum(axis=1) #Cumulative positional value per junction
    plot(junction_sum, numzeros/float(junctions.shape[1]), '*')
    _save_or_show(plotpath, plotname)


def score_ecdf(ecdf):
    "Give a score to a ecdf calculating the deviation from the 45 degree line"
    return sum(abs(linspace(0, 1, num=len(ecdf))-ecdf))


def plot_negbinomial_fit(mean_junc, std_junc, fit_function, plotpath, plotname):
    #plot the fit of the line 
    xlabel("Mean")
    ylabel("Std")
    xlim(0, max(mean_junc)*1.1) #adjust the x axis of the plot
    ylim(0, max(std_junc)*1.1) #adjust the y axis of the plot
    plot(mean_junc, std_junc, '*') #the mean and std for e very junction that passes the filter
    plot(mean_junc, fit_function(mean_junc), '-r')
    _save_or_show(plotpath, plotname)


def get_pvalues(junctions, a, b, dispersion):
    pvalues = []
    for i, junction in enumerate(junctions):
        if junction.any():
            junction_value = junction.mean()
            r, p = func2nb(a, b, junction_value, dispersion)
            my_nb = nbinom(r, 1-p)
            pval = 1-my_nb.cdf(junction_value)
        else:
            pval = 1 #if no reads, p-value is 1 

        pvalues.append(pval) 
         
    return pvalues

def adjust_fit(starting_a, b, junctions, precision, previous_score, dispersion, final=False, logger=None):
    previous_a = -1
    if logger: logger.info("Starting from %s with precision %s"%(starting_a, precision))
    for corrected_a in arange(starting_a, 0, -precision): #since we are reducing the "a" from the fit and the problem is too much variability, we expect optimization to be getting the "a" below 
        pvalues = get_pvalues(junctions, corrected_a, b, dispersion)
        ecdf = get_ecdf(pvalues)
        score = score_ecdf(ecdf)
        if previous_score < score: #the best fit are previous_a and previous_score
            if previous_a == -1:
                return corrected_a, score, ecdf, pvalues
            else:
                return previous_a, previous_score, previous_ecdf, previous_pvalues

        previous_a = corrected_a
        previous_score = score
        previous_ecdf = ecdf
        previous_pvalues = pvalues
        pvalues = []    

    if logger: logger.warning("WARNING: Something is wrong, please contact Biociphers!")
    return corrected_a, score, ecdf, pvalues #this return should not be hit

def fit_nb(junctions, outpath, plotpath, gcnorm=True, trim=True, minnonzero=5, plotmapzeros=False, discardb=False, nbdisp=0.1, logger=None, bval=False):
    if logger: logger.info("NBFit: Plots will be drawn in %s..."%plotpath)
    #normalize the junctions
    junctions = masked_less(junctions, 1) #mask everything below zero 
    mean_junc = junctions.mean(axis=1)
    std_junc = junctions.std(axis=1)
    #linear regression, retrieve the a and the b plus 
    a, b = polyfit(mean_junc, std_junc, 1)
    if discardb:
        if logger: logger.debug("Discarded b from the polyfit")
        if bval:
            print "fitnb, b=1"
            b = 1
        else:
            print "fitnb, b=0"
            b = 0

    print "GET_PVALS",b
    pvalues = get_pvalues(junctions, a, b, nbdisp)
    ecdf    = get_ecdf(pvalues)
    xlabel("P-value")
    ylabel("non_corrected ECDF")
    plot(linspace(0, 1, num=len(ecdf)), ecdf)
    plot([0, 1], 'k')
    _save_or_show(plotpath, "NON-Corrected ECDF b_%s"%b)


    if logger: logger.info("Fitting function: y = x*a+b. a=%.5f b=%.5f"%(a, b))
    fit_function = poly1d([a, b])
    nb_r, nb_p, matching_x = nb_from_func(fit_function, max(mean_junc), dispersion=nbdisp) #We calculate both r and p parameters of the negative binomial distribution along with the function
    plot_negbinomial_fit(mean_junc, std_junc, fit_function, plotpath, "Before correction")    
    #pvalue study


    #find the corresponding NB parameters to the junction mean
    score = sys.maxint
    corrected_a = a
    prev_precision = 0
    precision_values = [0.1, 0.01]
    if discardb:
        if logger: logger.debug("Discarded b from the polyfit")
        if bval:
            print "fitnb, b=1"
            b = 1
        else:
            print "fitnb, b=0"
            b = 0

    for i, precision in enumerate(precision_values):
        corrected_a, score, ecdf, pvalues = adjust_fit(corrected_a, b, junctions, precision, score, nbdisp, logger)
        if logger: logger.info("Corrected to %.5f with precision %s. Current score is %.5f"%(corrected_a, precision, score))
        if i+1 != len(precision_values): #go "up" in the scale so we dont miss better solution
            corrected_a += precision-precision_values[i+1]
            pvalues = get_pvalues(junctions, corrected_a, b, nbdisp)
            ecdf    = get_ecdf(pvalues)
            score   = score_ecdf(ecdf)

    if logger: logger.debug("Final Score: %.3f Final function A: %.3f"%(score, corrected_a))
    xlabel("P-value")
    ylabel("ECDF")
    plot(linspace(0, 1, num=len(ecdf)), ecdf)
    plot([0, 1], 'k')
    
    _save_or_show(plotpath, "Corrected ECDF b_%s"%b)

    fit_function = poly1d([corrected_a, b])
    if logger: logger.debug("Calculating the nb_r and nb_p with the new fitted function")
    nb_r, nb_p, matching_x = nb_from_func(fit_function, max(mean_junc), dispersion=nbdisp) 
    plot_negbinomial_fit(mean_junc, std_junc, fit_function, plotpath, "After correction") 

    #Save everything into pickle object to maybe later reuse?
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    pickle.dump(fit_function, open("%s/fitfunc.pickle"%outpath, 'w'))

    return fit_function
