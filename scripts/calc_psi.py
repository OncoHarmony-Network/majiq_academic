import sys
import os
import pickle
import argparse
from random import choice
from math import log

from scipy.io import loadmat
from pylab import *
import numpy as np
from matplotlib import rcParams
from scipy.stats import pearsonr
from numpy.random import dirichlet

"""
Calculate PSI values
"""
DEBUG = False
TESTBREAK = 50
BINS = linspace(0, 1, num=100)


def calc_psi(alpha, n, debug, *samples_events):
    "Expects N 2D matrices in samples_events"
    psi_matrix = []
    num_events = samples_events[0].shape[0]
    for event_num, event_samples in enumerate(range(samples_events[0].shape[0])): #we iterate through the second dimension of the matrix, which corresponds to the paired samples per event for different experiments    
        if event_num % 50 == 0:
            print "event %s..."%event_num,
            sys.stdout.flush()

        if DEBUG and event_num == TESTBREAK: break

        if len(event_samples.shape) == 1: #only one sample
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
    print
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

    print samples.shape
    psi_scores = calc_psi(args.alpha, args.n, samples)  
    pickle.dump(psi_scores, open("%s%s_vs_%s_psivalues.pickle"%(args.output, args.name1, args.name2), 'w'))
    print "Done."

if __name__ == '__main__':
    main()






