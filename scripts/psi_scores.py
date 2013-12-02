import sys
import os
import pickle
import argparse
from math import log

from pylab import *
from matplotlib import rcParams
from scipy.spatial.distance import cityblock

"""
Calculate PSI values
"""
DEBUG = False
TESTBREAK = 50
BINS = linspace(0, 1, num=100)

def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname), bbox_inches='tight') 
        clf()
    else:
        show()


def calculate_dkl(p, q):
    """
    Calculates distance between two distributions:

    Dkl(P|Q) = sum_i ln(P(i)/Q(i))*P(i)
    """
    pseudo = 0.001
    p += pseudo
    q += pseudo
    left = log(array(p/q))
    return (left*p).sum(axis=1)

def calculate_l1(p, q):
    return (abs(p - q)).sum(axis=1)

def calculate_ead(psi_samples):
    """
    P(S = |PSI_1 - PSI_2|) = sum_psi_1 sum_psi_2 p(psi_1)*p(psi_2)*|psi_1 - psi_2| 

    Expected Absolute Difference = EAD 
    """
    sample1 = psi_samples[:, 0]
    sample2 = psi_samples[:, 1]
    score = 0
    for event_num in xrange(sample1.shape[0]):
        for i in xrange(sample1.shape[1]):
            for j in xrange(sample2.shape[1]):
                psi_val1 = BINS[i]
                psi_val2 = BINS[j]
                score += sample1[event_num][i]*sample2[event_num][j]*abs(psi_val1 - psi_val2)
            
    return score
    #cellcalc = sample1*sample2*abs(sample1 - sample2)
    #return (cellcalc.sum(axis=1)).sum(axis=0)


def main():
    """
    Script for initial testing of the MAJIQ algorithms for sampling and initial PSI values generator. 

    TODO: Many functions from this script will be extracted for general usage in the pipeline. 
    """
    parser = argparse.ArgumentParser() 
    parser.add_argument('psivalue1', help='Path for psi pickles to evaluate')
    parser.add_argument('psivalue2', help='Path for psi pickles to evaluate')
    parser.add_argument('--name', default="")
    parser.add_argument('--plotpsidist', default=False, help="Plot the PSI distributions for ALL junctions. Slow.")
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    args = parser.parse_args()

    psivalue1 = pickle.load(open(args.psivalue1))
    psivalue2 = pickle.load(open(args.psivalue2))
    
    #add together the inclusion values
    inclusion1 = psivalue1[:, 0] #Inclusion values 
    inclusion2 = psivalue2[:, 0] #Inclusion values 
    """
    print "\nAnalysis of %s\n------------------------------------------------------------"%val
    print "Number of junctions:", psi_scores.shape[0]         
    print "Number of Experiments:", psi_scores.shape[1] 
    print "Number of BINS:", psi_scores.shape[2]
    print
    """
    if args.plotpsidist:
        for i in range(psi_scores.shape[0]):
            print "%s..."%i,
            sys.stdout.flush()
            subplot(2,1,1)
            xlim(-0.1,1.1)
            plot(BINS[:-1], psi_scores[i, 0])
            subplot(2,1,2)
            xlim(-0.1,1.1)
            plot(BINS[:-1], psi_scores[i, 1])
            _save_or_show(args.plotpath, "psi_%s"%(i))

    dkl_junctions = calculate_dkl(inclusion1, inclusion2)
    mean_dkl = mean(dkl_junctions)
    var_dkl = var(dkl_junctions)
    print "DKL: Mean: %.5f Var: %.6f Acum: %.3f"%(mean_dkl, var_dkl, sum(dkl_junctions))

    l1_junctions = calculate_l1(inclusion1, inclusion2)

    print "L1: Mean: %.5f Var: %.6f Acum: %.3f"%(mean(l1_junctions), var(l1_junctions), sum(l1_junctions))



if __name__ == '__main__':
    main()

