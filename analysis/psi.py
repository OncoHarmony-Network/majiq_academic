import sys
import os
import pickle
import argparse
from random import choice
import operator
from collections import defaultdict

from pylab import *
from scipy.special import gamma #deprecation WARNING comes from this import!!! 
from scipy.stats import pearsonr, binom_test
from numpy.random import dirichlet
from numpy import rollaxis
from matplotlib import rcParams

import analysis.filter as majiq_filter
import analysis.adjustdelta as majiq_delta
import analysis.sample as majiq_sample

"""
Calculate and manipulate PSI and Delta PSI values
"""
BSIZE = 0.025 #TODO To parameters
BINS = arange(0, 1, BSIZE) # The bins for PSI values. With a BSIZE of 0.025, we have 40 BINS
BINS_CENTER = arange(0+BSIZE/2, 1, BSIZE) #The center of the previous BINS. This is used to calculate the mean value of each bin.


def plot_matrix(matrix, my_title, plotname, plotpath):
    clf()
    ax = subplot(1,1,1)
    title(my_title)
    imshow(matrix)
    xlabel(u"PSI i")
    ylabel(u"PSI j")
    ax.set_xticklabels([0, 0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0, 0.25, 0.5, 0.75, 1])

    _save_or_show(plotpath, plotname=plotname)

def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname), bbox_inches='tight') 
        clf()
    else:
        show()


def median_psi(junctions, discardzeros=True):
    "Calculates the median PSI for all events"
    medians = []
    for junction in junctions:
        if discardzeros:
            junction = junction[junction!=0] #a junction array without the zeroes

        medians.append(median(junction))

    return array(medians)

def empirical_delta_psi( lsv_list1, lsv_list2, logger=None):
    """Simple PSI calculation without involving a dirichlet prior, coming from reads from junctions"""

#    if  not logger: logger.info("Calculating PSI for 'best set'...")

    delta_psi = []
    for idx, lsv in enumerate(lsv_list1):
        psi1 = np.zeros(shape=len(lsv), dtype=np.dtype('float'))
        psi2 = np.zeros(shape=len(lsv), dtype=np.dtype('float'))
        for ii, rate in enumerate(lsv):
            val = rate /  np.sum(lsv)
            if isnan(val): val = 0.5
            psi1[ii] = val
        
        for ii, rate in enumerate(lsv_list2[idx]):
            val = rate /  np.sum(lsv_list2[idx])
            if isnan(val): val = 0.5
            psi2[ii] = val

        sys.stdout.flush()

        delta_psi.append( psi1 - psi2 )
 #   if logger: logger.info("Calculating delta PSI for 'best set'...")
    return delta_psi 

def simple_psi(inc, exc):
    """Simple PSI calculation without involving a dirichlet prior, coming from reads from junctions"""
    psi = inc/(exc+inc)
    psi[isnan(psi)] = 0.5 #if NaN, is because exc+inc = 0. If we know nothing, then we don't know if its 0 (exclusion) or 1 (inclusion)
    return psi 

def reads_given_psi_lsv(lsv_junc, psi_space):
    #P(vector_i | PSI_i)
    "We do a simple binomial test to evaluate how probable is the data given a PSI range"
    ret = []
    lsv = lsv_junc.sum(axis=1)
    for idx in range(lsv.shape[0]):
        event = []
        for psi_val in psi_space:
            event.append(binom_test(lsv[idx], lsv.sum(), p=psi_val))

        ret.append(array(event) / sum(event))
        print ret
    return array(ret).reshape(-1, len(psi_space)) 

def reads_given_psi(inc_samples, exc_samples, psi_space):
    #P(vector_i | PSI_i)
    "We do a simple binomial test to evaluate how probable is the data given a PSI range"
    ret = []
    inc = inc_samples.sum(axis=1)
    exc = exc_samples.sum(axis=1)
    for i in xrange(inc.shape[0]):
        event = []
        for psi_val in psi_space:
            event.append(binom_test(inc[i], exc[i]+inc[i], p=psi_val))

        ret.append(array(event) / sum(event))
    
    return array(ret).reshape(-1, len(psi_space)) 


class DirichletCalc:
    def __init__(self):
        self.cache = defaultdict(float)

    def pdf(self, x, alpha):
        k = x
        k.extend(alpha)
        key = " ".join(map(str, k))
        if self.cache[key]: #try to figure out if we already calculated this pdf
            return self.cache[key] 
        else: 
            #formula taken from stackoverflow "How to calculate dirichlet PDF", author based it on Wikipedia's PDF definition
            ret = gamma(sum(alpha)) / reduce(operator.mul, [gamma(a) for a in alpha]) * reduce(operator.mul, [x[i]**(alpha[i]-1.0) for i in xrange(len(alpha))])
            self.cache[key] = ret
            return ret


def recalibrate_delta(deltapsi):
    #TODO make deltaPSI follow the following binning system
    arange(-98.75, 100, 2.5)

def lsv_psi(samples_events, name, alpha, n, debug):
    "Given a set of matching inclusion and exclusion samples, calculate psi, save it in disk, and return the psi-per-juntion matrix"
    
    psi_scores = []
    dircalc = DirichletCalc() 
    for i, lsv in enumerate(samples_events):
        if i % 50 == 0:
            print "event %s..."%i,
            sys.stdout.flush()
        if debug > 0 and i == debug: break
        psi = np.zeros(shape=(lsv.shape[0],BINS.shape[0]), dtype=np.float)
        #if debug: print "Paired samples to dirichlet..."
        #sampling PSI by pairing the samples of the previous step sequentially
        for idx, junc in enumerate(lsv):
            total_acum = 0.
            acum_samples = np.zeros(shape=(BINS.shape[0]))
            aggr = np.zeros(shape=(junc.shape[0]))
            for xidx, xx in enumerate(lsv):
                if idx == xidx : continue
                aggr += xx

            samples    = np.ndarray(shape=(2,junc.shape[0]))
            samples[0,:] = junc
            samples[1,:] = aggr

            for paired_samples in samples.T:

                dir_pdf = [dircalc.pdf([x, 1-x], alpha+paired_samples) for x in BINS_CENTER]
                dir_pdf = np.asarray(dir_pdf)
                acum_samples += dir_pdf
                total_acum += sum(dir_pdf) 

            psi[idx]=acum_samples/total_acum

        #if debug: print "Dividing by total acum..."
        psi_scores.append( psi )
        #print "Junction %s PSI distribution: %s sum_N: %s"%(i, psi_matrix[-1], sum(psi_matrix[-1]))


    return psi_scores

def calc_psi(inc_samples, exc_samples, name, alpha, n, debug, psiparam):
    "Given a set of matching inclusion and exclusion samples, calculate psi, save it in disk, and return the psi-per-juntion matrix"
    
    print inc_samples.shape
    samples = vstack([inc_samples, exc_samples]).reshape(2, inc_samples.shape[0], inc_samples.shape[1])
    psi_scores = calc_dirichlet(alpha, n, samples, debug=debug, psiparam=psiparam)
    if psiparam:
        return psi_scores
    else:
        return psi_scores[:,0] #psi_scores[:,1] is PSE


def mean_psi(psi_events):
    "Calculate the mean for every junction. Used for delta PSI calculation."
    ret = []
    for psi_dist in psi_events:
        #print "PSI_DIST", psi_dist
        #print "sum(PSI_DIST)", sum(psi_dist)
        ret.append(sum(psi_dist*BINS_CENTER))

    return array(ret)


def calc_dirichlet(alpha, n, samples_events, debug=False, psiparam=False):
    "Expects 3 dimensional matrix in samples_events"
    psi_matrix = []
    dircalc = DirichletCalc() 
    if psiparam:
        for i, event_samples in enumerate(rollaxis(samples_events, 1)): #The second dimension of the matrix corresponds to the paired samples per event (3rd dimension) for different experiments (1st dimension)       
            if i % 5 == 0:
                print "event %s..."%i,
                sys.stdout.flush()

            if debug > 0 and i == debug: break
            #if debug: print "Paired samples to dirichlet..."
            #sampling PSI by pairing the samples of the previous step sequentially
            total_acum = 0.
            acum_samples = array([0]*BINS) 
            for h, paired_samples in enumerate(event_samples.T):
                dir_pdf = [dircalc.pdf([x, 1-x], alpha+paired_samples) for x in BINS_CENTER]
                acum_samples += dir_pdf
                total_acum += sum(dir_pdf) 

            #if debug: print "Dividing by total acum..."
            psi_matrix.append(acum_samples/total_acum)
            #print "Junction %s PSI distribution: %s sum_N: %s"%(i, psi_matrix[-1], sum(psi_matrix[-1]))
    else:
        for i, event_samples in enumerate(rollaxis(samples_events, 1)): #we iterate through the second dimension of the matrix, which corresponds to the paired samples per event for different experiments        
            #This is sampling the PSI. Instead of doing this, we want to fit a parametric form.
            if i % 50 == 0:
                print "event %s..."%i,
                sys.stdout.flush()

            if debug > 0 and i == debug: break

            if len(event_samples.shape) == 1: #only one sample (this is only used if bootstrapping of reads is deactivated)
                event_psi_samples = dirichlet(event_samples+alpha, n)
            else:
                event_psi_samples = []
                #sampling PSI by pairing the samples of the previous step sequentially (to gain execution time)
                for paired_samples in event_samples.T:
                    event_psi_samples.extend(dirichlet(paired_samples+alpha, n))

            #discretization step. Get the psi samples and transform them into a histogram-like distribution
            event_psi_discrete = []
            for psi_dist in array(event_psi_samples).transpose():
                my_bins = list(BINS)
                my_bins.extend([1]) #extend because:  If `bins` is a sequence,it defines the bin edges, including the rightmost edge, allowing for non-uniform bin widths (form histogram docs)
                counts, limits = histogram(psi_dist, bins=my_bins) 
                event_psi_discrete.append(counts/float(sum(counts)))

            #print event_psi_discrete[-1], len(event_psi_discrete), len(event_psi_discrete[-1])
            psi_matrix.append(event_psi_discrete)
            #print "Junction %s PSI distribution:"%i, psi_matrix[-1]

    psi_matrix = array(psi_matrix)
    return psi_matrix


def gen_prior_matrix( pip, lsv_exp1, lsv_exp2, output ):

    #Start prior matrix
    pip.logger.info("Calculating prior matrix...")
    numbins = 20 #half the delta bins TODO to parameters

    pip.logger.info("Calculate jefferies matrix...")
    dircalc = DirichletCalc()
    #Adjust prior matrix with Jefferies prior        
    jefferies = []
    psi_space = linspace(0, 1-pip.binsize, num=numbins) + pip.binsize/2
    for i in psi_space:
        jefferies.append([])
        for j in psi_space:
            jefferies[-1].append(dircalc.pdf([i, 1-i, j, 1-j], [pip.alpha, pip.alpha, pip.alpha, pip.alpha]))

    #jefferies = array([dircalc.pdf([x, 1-x], [0.5, 0.5]) for x in psi_space])
    jefferies = array(jefferies)
    jefferies /= sum(jefferies)
    plot_matrix(jefferies, "Jefferies Matrix", "jefferies_matrix", pip.plotpath)
    if pip.synthprior:
        #Use a synthetic matrix to generate the values
        prior_matrix = [] 
        uniform = pip.prioruniform/numbins 
        mydist = norm(loc=0, scale=pip.priorstd)
        norm_space = linspace(-1, 1-pip.binsize, num=numbins*2) + pip.binsize/2
        pdfnorm = mydist.pdf(norm_space)
        newdist = (pdfnorm+uniform)/(pdfnorm+uniform).sum()
        plot(linspace(-1, 1, num=len(list(pdfnorm))), pdfnorm)
        _save_or_show(pip.plotpath, plotname="prior_distribution")
        #generate the matrix
        for i in xrange(numbins):
            prior_matrix.append(list(newdist[numbins-i:(numbins*2)-i]))

        prior_matrix = array(prior_matrix)
        prior_matrix /= sum(prior_matrix) #renormalize so it sums 1
        pip._get_delta_info(newdist, norm_space)
        plot_matrix(prior_matrix, "Prior Matrix (before Jefferies)", "prior_matrix_no_jefferies", pip.plotpath)

    elif not pip.jefferiesprior:
        #Using the empirical data to get the prior matrix
        pip.logger.info('Filtering to obtain "best set"...')

        filtered_lsv1 = majiq_filter.lsv_quantifiable(lsv_exp1, minnonzero=5, min_reads=20, logger=pip.logger)
        filtered_lsv2 = majiq_filter.lsv_quantifiable(lsv_exp2, minnonzero=5, min_reads=20, logger=pip.logger)


#        print "FILTER1",filtered_lsv1[1]
#        print "FILTER2",filtered_lsv2[1]

        ids1 = set([xx[1] for xx in filtered_lsv1[1]])
        ids2 = set([xx[1] for xx in filtered_lsv2[1]])
        matched_names = ids1.intersection(ids2)
        best_set_mean1 = [[],[]]
        best_set_mean2 = [[],[]]
    
        for ii in matched_names:
            for idx, nm in enumerate(filtered_lsv1[1]):
                if nm[1] == ii:
                    best_set_mean1[0].append(majiq_sample.mean_junction(filtered_lsv1[0][idx]))
                    best_set_mean1[1].append(filtered_lsv1[1][idx])
                    break
            for idx, nm in enumerate(filtered_lsv2[1]):
                if nm[1] == ii:
                    best_set_mean2[0].append(majiq_sample.mean_junction(filtered_lsv2[0][idx]))
                    best_set_mean2[1].append(filtered_lsv2[1][idx])
                    break

        pip.logger.info("'Best set' is %s events (out of %s)"%(len(best_set_mean1), len(lsv_exp1)))
        best_delta_psi = empirical_delta_psi(best_set_mean1[0], best_set_mean2[0])
        best_delta_psi = np.concatenate(best_delta_psi)

        pip.logger.info("Parametrizing 'best set'...")
        mixture_pdf = majiq_delta.adjustdelta(best_delta_psi, output, plotpath=pip.plotpath, title=" ".join(pip.names), numiter=pip.iter, breakiter=pip.breakiter, V=pip.V, logger=pip.logger)

        pickle.dump(mixture_pdf, open("%s%s_%s_bestset.pickle"%(output, pip.names[0], pip.names[1]), 'w'))

        prior_matrix = []
        for i in xrange(numbins):
            prior_matrix.extend(mixture_pdf[numbins-i:(numbins*2)-i])
        prior_matrix = array(prior_matrix).reshape(numbins, -1)

        #some info for later analysis
#        pickle.dump(event_names, open("%s%s_%s_eventnames.pickle"%(output, self.names[0], self.names[1]), 'w')) 
        if not pip.jefferiesprior:
            plot_matrix(prior_matrix, "Prior Matrix (before Jefferies)", "prior_matrix_no_jefferies", pip.plotpath)

        #Calculate prior matrix
        pip.logger.info("Adding a Jefferies prior to prior (alpha=%s)..."%(pip.alpha))
        #Normalize prior with jefferies
        if pip.jefferiesprior:
            pip.logger.info("Using the Uniform distribution + Jefferies...")
            prior_matrix = jefferies + (pip.prioruniform/numbins)
        else: 
            prior_matrix *= jefferies 

        prior_matrix /= sum(prior_matrix) #renormalize so it sums 1
        plot_matrix(prior_matrix, "Prior Matrix", "prior_matrix", pip.plotpath)
        pip.logger.info("Saving prior matrix for %s..."%(pip.names))
        pickle.dump(prior_matrix, open("%s%s_%s_priormatrix.pickle"%(output, pip.names[0], pip.names[1]), 'w'))
    
    return psi_space, prior_matrix



#deprecated
def sample_psi(psi_scores):
    """
    Input is all junctions PSI distributions for 1 replica
    """
    samples = []
    for pval, limits in psi_scores:
        event_samples = []
        sample_pos = multinomial(100, pval)
        for p in sample_pos:
            event_samples.append(limits[p])

        samples.append(mean(event_samples))

    return array(samples)


def main():
    """
    Script for initial testing of the MAJIQ algorithms for sampling and initial PSI values generator. 

    TODO: Many functions from this script will be extracted for general usage in the pipeline. 
    """
    parser = argparse.ArgumentParser() 
    parser.add_argument('samples', nargs='+', help='Path for samples of conditions 1 to N')
    parser.add_argument('--n', default=1, type=int, help='Number of PSI samples per sample paired') 
    parser.add_argument('--alpha', default=0.5, type=float, help='Alpha hyperparameter for the dirichlet distribution') 
    parser.add_argument('--output', required=True, help="Path to save the results to.")
    parser.add_argument('--name1', default='Inc')
    parser.add_argument('--name2', default='Exc')
    args = parser.parse_args()

    print "Loading samples..."
    samples = []
    for sample in args.samples:
        samples.append(pickle.load(open(sample)))

    samples = vstack(samples)
    print "Calculating PSI for %s and %s..."%(args.name1, args.name2)
    psi_scores = calc_dirichlet(args.alpha, args.n, samples)  
    pickle.dump(psi_scores, open("%s%s_vs_%s_psivalues.pickle"%(args.output, args.name1, args.name2), 'w'))
    print "Done."



if __name__ == '__main__':
    main()

