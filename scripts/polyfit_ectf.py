import sys
import os
import argparse
import pickle

from scipy.io import loadmat
from scipy.stats import nbinom
from pylab import *

#### SNIPPETS  #######
#Get the names of the columns in an ndarray
#print my_mat['rand10k'].dtype.names


### Plot the samples from the curve fit #####
        #for i in xrange(len(nb_r)):
        #    z = negative_binomial(nb_r[i], nb_p[i], size=200) #from the NB distribution of each point, we sample a few values and then plot them to see how they fit 
        #    plot(z.mean(), z.std(), "+k", markersize=10)


def nb_from_func(poly_func, max_line):
    """
    From a linear function that describes the distribution of junctions mean/variance, get the NB values

    We use a combined approach to avoid having values of r < 0. 
    For small values, we fit the r *global* negative binomial with a dispersion parameter (1/r in the NB formula)
    of 0.1, which is a bit bigger variance than a simple Poisson. 

    """
    dispersion = 0.1 #for low coverage junctions
    a, b = poly_func.c
    #x = arange(0.01, max_line, 0.01)
    r = []
    p = []
    #points = arange(0.01, max_line, 0.01)
    points = linspace(0.01, max_line, num=800)
    for x in points:
        if (a*x)**2 > (x+dispersion*x**2):
            r_val = x**2 / ((a*x+b)**2 - x)
        else:
            r_val = 1/dispersion

        r.append(r_val)
        p.append(r_val / (r_val+x))

    #r = (x**2) / ((a*x+b)**2 - x) 
    #p = r / (r+x)
    return r, p, points

def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if a flag was given"""
    if plotname:
        title(plotname)

    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname.replace(" ", "_")), bbox_inches='tight') #WNo spaces allowed, underscores!
        clf()
    else:
        show()

def _numzeros(junctions):
    "Obtain the number of zeros for every junction"
    return (junctions == 0).sum(axis=1)

def _numnonzeros(junctions):
    "Obtain the number NON zero positions for every junction"
    return (junctions != 0).sum(axis=1)

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


def get_pvalues(sum_junctions, a, b, nonzeros):
    pvalues = []
    b = 0
    for i, junction in enumerate(sum_junctions):
        junction_value = junction/nonzeros[i]
        r = junction_value**2 /  ((a*junction_value+b)**2 - junction_value)
        p = 1 - (r / (r+junction_value))
        my_nb = nbinom(r, p)
        pval = 1-my_nb.cdf(junction_value)
        #if not isnan(pval):  #only include nonnan palues
        pvalues.append(pval) 
        #print "Discarded NAN pvalue: %s out of %s junctions (%.2f%%)"%(len(sum_junctions)-len(pvalues), len(sum_junctions), (len(sum_junctions)-len(pvalues))/float(len(sum_junctions))*100)
    
    return pvalues

def adjust_fit(starting_a, b, sum_junctions, nonzeros, precision, previous_score, final=False):
    previous_a = -1
    print "Starting from %s with precision %s"%(starting_a, precision)
    for corrected_a in arange(starting_a, 0, -precision): #since we are reducing the "a" from the fit and the problem is too much variability, we expect optimization to be getting the "a" below 
        pvalues = get_pvalues(sum_junctions, corrected_a, b, nonzeros)
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

    print "Do I hit??"
    return corrected_a, score, ecdf, pvalues #I am not sure if this return will be hit at all

def norm_junctions(junctions, gc_factors=None, gcnorm=False, trim=False, minnonzero=False, plotmapzeros=False):
    if gcnorm:
        print "Normalizing by GC factor"
        junctions = junctions*gc_factors

    #number of zeros before
    if trim:
        numzeros = _numzeros(junctions) 
        print "Previous junction length: %s Total number of zeros: %s\n"%(junctions.shape[1], sum(numzeros))
        if plotmapzeros:
            plot_mappability_zeros(junctions, plotpath, numzeros, "before_trim")

        junctions = junctions[:,trim:junctions.shape[1]-trim] #trim the junctions according to the flag
        if plotmapzeros:
            plot_mappability_zeros(junctions, plotpath, "after_trim (%s bases)"%trim)

    if minnonzero:
        print "Filter out the junctions with less than %s positions bigger than 0"%minnonzero
        num_before =  junctions.shape[0]
        numzeros = _numzeros(junctions) 
        pass_threshold = (numzeros < junctions.shape[1]-minnonzero)
        junctions = junctions[pass_threshold]
        print "Before: %s junctions. Now: %s junctions. Discarded %s"%(num_before, junctions.shape[0], num_before-junctions.shape[0])    

    return junctions

def process(path, args):
    #copied from Jordis script, this will have to go at some point
    my_mat = loadmat(path)
    gc_factors = my_mat['rand10k']['gc_val'][0, 0]
    junctions = my_mat['rand10k']['cov'][0, 0] # We do [0,0] because matlab 
    outpath = "%s/%s"%(args.output, os.path.basename(path).replace(".mat", ""))
    plotpath = "%s/%s"%(args.plotpath, os.path.basename(path).replace(".mat", ""))
    print "Results will be written in %s..."%outpath
    print "Plots will be drawn in %s..."%plotpath
    #normalize the junctions
    junctions = norm_junctions(junctions, gc_factors, args.gcnorm, args.trim, args.minnonzero, args.plotmapzeros)
    mean_junc = junctions.mean(axis=1)
    std_junc = junctions.std(axis=1)
    #linear regression, retrieve the a and the b plus 
    a, b = polyfit(mean_junc, std_junc, 1)
    print "Fitting function: y = x*a+b. a=%s b=%s"%(a, b)
    fit_function = poly1d([a, b])
    nb_r, nb_p, matching_x = nb_from_func(fit_function, max(mean_junc)) #We calculate both r and p parameters of the negative binomial distribution along the function
    plot_negbinomial_fit(mean_junc, std_junc, fit_function, plotpath, "Before correction")    
    #pvalue study
    nonzeros = _numnonzeros(junctions) 
    sum_junctions = junctions.sum(axis=1) 
    #find the corresponding NB parameters to the junction mean
    score = sys.maxint
    corrected_a = a
    prev_precision = 0
    precision_values = [0.1, 0.01]
    if args.discardb:
        print "Discarded b from the polyfit"
        b = 0

    for i, precision in enumerate(precision_values):
        corrected_a, score, ecdf, pvalues = adjust_fit(corrected_a, b, sum_junctions, nonzeros, precision, score)
        print "Corrected to %s with precision %s. Current score is %s\n"%(corrected_a, precision, score)
        if i+1 != len(precision_values): #go "up" in the scale so we dont miss better solution
            corrected_a += precision-precision_values[i+1]
            pvalues = get_pvalues(sum_junctions, corrected_a, b, nonzeros)
            ecdf    = get_ecdf(pvalues)
            score   = score_ecdf(ecdf)

    print "Final Score: %.3f Final function A: %.3f"%(score, corrected_a)
    xlabel("P-value")
    ylabel("ECDF")
    plot(linspace(0, 1, num=len(ecdf)), ecdf)
    plot([0, 1], 'k')
    
    _save_or_show(plotpath, "Corrected ECDF")

    fit_function = poly1d([corrected_a, b])
    print "Calculating the nb_r and nb_p with the new fitted function"
    nb_r, nb_p, matching_x = nb_from_func(fit_function, max(mean_junc)) #
    plot_negbinomial_fit(mean_junc, std_junc, fit_function, plotpath, "After correction") 
    #Save everything in pickle objects
    pickle.dump(ecdf, open("%s_ecdf.pickle"%outpath, 'w'))
    pickle.dump(fit_function, open("%s_fitfunc.pickle"%outpath, 'w'))
    pickle.dump(nb_r, open("%s_nb_r.pickle"%outpath, 'w'))
    pickle.dump(nb_p, open("%s_nb_p.pickle"%outpath, 'w'))
    pickle.dump(matching_x, open("%s_nb_index.pickle"%outpath, 'w'))
    pickle.dump(junctions, open("%s_filtered_junctions.pickle"%outpath, 'w'))    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('matfiles', nargs='+', help='The matlab files that we read. Can be used as a glob.')
    parser.add_argument('--trim', default=0, type=int, help='Trim the borders of the junctions because of poor mappability')
    parser.add_argument('--minnonzero', default=0, type=int, help='Minimum number of positions ')
    parser.add_argument('--gcnorm', action='store_true',  default=False, help='Correct by GC content')
    parser.add_argument('--discardb', action='store_true',  default=False, help='Discard the b from the polynomial, since we expect our fit to start from 0, 0')
    parser.add_argument('--weightnorm', action='store_true',  default=False, help='Correct using the weight factor')
    parser.add_argument('--plotmapzeros', action='store_true', default=False, help='Plot the zeros study')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    args = parser.parse_args()
    for path in args.matfiles:
        print "\n\n\nProcessing %s..."%path
        process(path, args)

if __name__ == '__main__':
    main()



