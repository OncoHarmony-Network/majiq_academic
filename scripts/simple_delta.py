import sys
import os
import argparse
import pickle

from pylab import *
from numpy.random import choice
from scipy.io import loadmat
from numpy.ma import masked_less
 
from analysis.filter import discardhigh, discardlow, discardminreads, discardmaxreads, norm_junctions

BINS = linspace(0, 1, num=99)



def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname), bbox_inches='tight') 
        clf()
    else:
        show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('matpath', help='Path for .mat with junctions')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--title', default=None, help='') 
    parser.add_argument('--minnonzero', default=0, type=int, help='Minimum number of positive positions to consider the junction')   
    parser.add_argument('--maxnonzero', default=0, type=int, help='Maximum number of positive positions to consider the junction') 
    parser.add_argument('--minreads', default=0, type=int, help='Minimum number of reads combining all positions in a junction to be considered') 
    parser.add_argument('--maxreads', default=0, type=int, help='Maximum number of reads combining all positions in a junction to be considered') 
    parser.add_argument('--orfilter', default=False, action='store_true', help='When filtering, select sets of junctions where at least one passes the filter, instead of all passing the filter.')
    parser.add_argument('--output', default=None, help="Path to save the results to.")
    args = parser.parse_args()
    #load matlab matrices and get names of experiments from file name (TODO:All this should change)
    my_mat = loadmat(args.matpath)
    name1, name2 = os.path.basename(args.matpath).split('.')[0].split("_") #Get the experiment names from the mat file
    main_title = "%s VS %s\n"%(name1, name2) #main title of the plot generation
    inclusion1 = my_mat[name1]["Inc"][0, 0][0, 0]['cov']
    exclusion1 = my_mat[name1]["Exc"][0, 0][0, 0]['cov']
    inclusion2 = my_mat[name2]["Inc"][0, 0][0, 0]['cov']
    exclusion2 = my_mat[name2]["Exc"][0, 0][0, 0]['cov']

    inclusion1_gc = my_mat[name1]["Inc"][0, 0][0, 0]['gc_val']
    exclusion1_gc = my_mat[name1]["Exc"][0, 0][0, 0]['gc_val']
    inclusion2_gc = my_mat[name2]["Inc"][0, 0][0, 0]['gc_val']
    exclusion2_gc = my_mat[name2]["Exc"][0, 0][0, 0]['gc_val']

    inclusion1 = norm_junctions(inclusion1, gc_factors=inclusion1_gc, gcnorm=True)
    inclusion2 = norm_junctions(inclusion2, gc_factors=inclusion2_gc, gcnorm=True)
    exclusion1 = norm_junctions(exclusion1, gc_factors=exclusion1_gc, gcnorm=True)
    exclusion2 = norm_junctions(exclusion2, gc_factors=exclusion2_gc, gcnorm=True)

    if args.maxnonzero:
        exclusion1, inclusion1, exclusion2, inclusion2 = discardhigh(args.maxnonzero, args.orfilter, True, exclusion1, inclusion1, exclusion2, inclusion2)
        
    if args.minnonzero:
        exclusion1, inclusion1, exclusion2, inclusion2 = discardlow(args.minnonzero, args.orfilter, True, exclusion1, inclusion1, exclusion2, inclusion2)

    if args.minreads:
        exclusion1, inclusion1, exclusion2, inclusion2 = discardminreads(args.minreads, args.orfilter, True, exclusion1, inclusion1, exclusion2, inclusion2)

    if args.maxreads:
        exclusion1, inclusion1, exclusion2, inclusion2 = discardmaxreads(args.maxreads, args.orfilter, True, exclusion1, inclusion1, exclusion2, inclusion2)

    pseudo = 0.001
    #mask everything under 0
    inclusion1 = masked_less(inclusion1, 0) 
    inclusion2 = masked_less(inclusion2, 0)
    exclusion1 = masked_less(exclusion1, 0) 
    exclusion2 = masked_less(exclusion2, 0)

    inc_reads1 = inclusion1.mean(axis=1)+pseudo
    inc_reads2 = inclusion2.mean(axis=1)+pseudo
    exc_reads1 = exclusion1.mean(axis=1)+pseudo
    exc_reads2 = exclusion2.mean(axis=1)+pseudo
    #totalreads = reads1 + pseudo * 2
    psi1 = inc_reads1/(exc_reads1+inc_reads1)
    psi2 = inc_reads2/(exc_reads2+inc_reads2)
    print psi1
    delta_psi = psi1 - psi2
    xlim(-1, 1)
    hist(delta_psi, bins = 60, histtype='step')

    if args.title:
        title(args.title)

    xlabel("Delta PSI")

    _save_or_show(args.plotpath, "deltapsi")



if __name__ == '__main__':
    main()

