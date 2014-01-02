import sys
import os
import pickle
import argparse
from random import choice
from scipy.special import gamma
import operator
from collections import defaultdict

from scipy.io import loadmat
from pylab import *
import numpy as np
from matplotlib import rcParams
from scipy.stats import pearsonr
from numpy.random import dirichlet


"""
Calculate and manipulate PSI and Delta PSI values
"""
BSIZE = 0.025 #TODO To parameters
BINS = arange(0, 1, BSIZE) # The bins for PSI values. With a BSIZE of 0.025, we have 40 BINS
BINS_CENTER = arange(0+BSIZE/2, 1, BSIZE) #The center of the previous BINS. This is used to calculate the mean value of each bin.

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


def calc_psi(inc_samples, exc_samples, name, output, alpha, n, debug, psinosample):
    "Given a set of matching inclusion and exclusion samples, calculate psi, save it in disk, and return the psi-per-juntion matrix"
    samples = vstack([inc_samples, exc_samples]).reshape(2, inc_samples.shape[0], inc_samples.shape[1])
    psi_scores = calc_dirichlet(alpha, n, samples, debug=debug, psinosample=psinosample)
    psi = psi_scores[:,0] #psi_scores[:,1] is PSE
    outpath = "%s_%s_psi.pickle"%(output, name)
    print "Writing PSI result into %s ..."%outpath
    pickle.dump(psi, open(outpath, 'w'))
    return psi


def mean_psi(psi_events):
    "Calculate the mean for every junction. Used for delta PSI calculation."
    ret = []
    for psi_dist in psi_events:
        ret.append(sum(psi_dist*BINS_CENTER[:-1]))

    return array(ret)


def calc_dirichlet(alpha, n, samples_events, debug=False, psinosample=False):
    "Expects 3 dimensional matrix in samples_events"
    psi_matrix = []
    dircalc = DirichletCalc() 
    if psinosample:
        for i, event_samples in enumerate(np.rollaxis(samples_events, 1)): #The second dimension of the matrix corresponds to the paired samples per event (3rd dimension) for different experiments (1st dimension)       
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
                total_acum += sum(dir_pdf) #TODO repasar

            #if debug: print "Dividing by total acum..."
            psi_matrix.append(acum_samples/total_acum)
            #print "Junction %s PSI distribution: %s sum_N: %s"%(i, psi_matrix[-1], sum(psi_matrix[-1]))
    else:
        for i, event_samples in enumerate(np.rollaxis(samples_events, 1)): #we iterate through the second dimension of the matrix, which corresponds to the paired samples per event for different experiments        
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
                counts, limits = histogram(psi_dist, bins=BINS)
                event_psi_discrete.append(counts/float(sum(counts)))

            #print event_psi_discrete[-1], len(event_psi_discrete), len(event_psi_discrete[-1])
            psi_matrix.append(event_psi_discrete)
            #print "Junction %s PSI distribution:"%i, psi_matrix[-1]

    return array(psi_matrix)

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

