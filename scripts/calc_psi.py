from matplotlib import use
use('Agg', warn=False)
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
import pipelines
import analysis.filter


"""
Calculate PSI values
"""
DEBUG = False
TESTBREAK = 50
BINS = linspace(0, 1, num=100)


def calc_psi(alpha, n, debug, *samples_events):
    """Expects N 2D matrices in samples_events"""
    psi_matrix = []
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

    parser.add_argument('files', nargs='+', help='Path for replicas of conditions 1 to N')
    parser.add_argument('--lsv', default=True, action='store_true', help='Execute pipeline for lsv')
    parser.add_argument('--tmp', default="/tmp/", help='Path to save the temporary files. [Default: %(default)s]')
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    parser.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    parser.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--debug', type=int, default=0, help="Activate this flag for debugging purposes, activates logger and jumps some processing steps.")
    parser.add_argument('--minreads', default=10, type=int, help='Minimum number of reads combining all positions in an event to be considered. [Default: %(default)s]')
    parser.add_argument('--minnonzero', default=5, type=int, help='Minimum number of start positions with at least 1 read for an event to be considered.')
    parser.add_argument('--markstacks', default=0.0000001, type=float, help='Mark stack positions. Expects a p-value. Use a negative value in order to disable it. [Default: %(default)s]')
    parser.add_argument('--tracklist', nargs='+', help='A list of identifiers to track in detail, for debugging purposes')

    parser.add_argument('--n', default=50, type=int, help='Number of PSI samples per sample paired')
    parser.add_argument('--alpha', default=0.5, type=float, help='Alpha hyperparameter for the dirichlet distribution')

    parser.add_argument('--gcnorm', default=None, help='Flag for GC normalization.')
    parser.add_argument('--nbdisp', default=0.1, type=float, help='Dispersion factor (used in junctions sampling).')
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration')
    parser.add_argument('--m', default=100, type=int, help='Number of samples')
    parser.add_argument('--files-contain-scores', dest="files_contain_scores", action='store_true', default=False, help='Feed the script with psis already computed.')

    args = parser.parse_args()


    methods = {
        'Poisson':                  {'discardzeros': 0, 'trimborder': False,   'nb': False},
        'Naive_Boots':              {'discardzeros': 0, 'trimborder': False,   'nb': False},
        'Naive_Boots_trim_borders': {'discardzeros': 1, 'trimborder': True,    'nb': False},
        'Naive_Boots_no_zeros':     {'discardzeros': 1, 'trimborder': False,   'nb': False},
        'Neg_Binomial':             {'discardzeros': 0, 'trimborder': False,   'nb': True},
        'Majiq':                    {'discardzeros': 1, 'trimborder': True,    'nb': True},
        'Majiq_with_zeros':         {'discardzeros': 0, 'trimborder': True,    'nb': True},
        'Majiq_no_stacks':          {'discardzeros': 1, 'trimborder': True,    'nb': True},
        'Majiq_padding_5':          {'discardzeros': 5, 'trimborder': True,    'nb': True},
        'Majiq_padding_10':         {'discardzeros': 10,'trimborder': True,    'nb': True},
        'Majiq_gc_norm':            {'discardzeros': 1, 'trimborder': True,    'nb': True},

    }

    args.__dict__.update(methods['Majiq_padding_5'])
    if args.files_contain_scores:
        psi_scores = []
        for file in args.files:
            psi_scores.append(pickle.load(open(file, 'r')))
    else:
        psi_scores = pipelines.calcpsi(args)

    print "Number of files analyzed: %d\nshapes: (%d, %s)\t(%d, %s)" % (len(psi_scores), len(psi_scores[0]), str(psi_scores[0][0][0].shape), len(psi_scores[1]), str(psi_scores[1][0][0].shape))
    lsv_match, match_info = analysis.filter.lsv_intersection(psi_scores[0], psi_scores[1])

    from os.path import basename
    output_file_name = "_vs_".join([basename(f) for f in args.files])

    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    pickle.dump((lsv_match, match_info), open("%s/%s_psivalues.pickle"%(args.output, output_file_name), 'w'))
    print "Done. Output file can be found at:\n%s/%s_psivalues.pickle" % (args.output, output_file_name)

if __name__ == '__main__':
    main()






