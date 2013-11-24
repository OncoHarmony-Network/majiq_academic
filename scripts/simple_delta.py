import sys
import os
import argparse
import pickle

from pylab import *
from numpy.random import choice
from scipy.io import loadmat
from numpy.ma import masked_less
 
from junction_sample import discardhigh, discardlow, discardminreads, discardmaxreads, norm_junctions

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
    parser.add_argument('--junctype', default='rand10k', help='The type of junction to analyze. (Inc, Exc or rand10k for now)')
    parser.add_argument('--orfilter', default=False, action='store_true', help='When filtering, select sets of junctions where at least one passes the filter, instead of all passing the filter.')
    parser.add_argument('--output', default=None, help="Path to save the results to.")
    args = parser.parse_args()
    #load matlab matrices and get names of experiments from file name (TODO:All this should change)
    my_mat = loadmat(args.matpath)
    name1, name2 = os.path.basename(args.matpath).split('.')[0].split("_") #Get the experiment names from the mat file
    main_title = "%s VS %s\n"%(name1, name2) #main title of the plot generation
    replica1 = my_mat[name1][args.junctype][0, 0][0, 0]['cov']
    replica2 = my_mat[name2][args.junctype][0, 0][0, 0]['cov']
    replica1_gc = my_mat[name1][args.junctype][0, 0][0, 0]['gc_val']
    replica2_gc = my_mat[name2][args.junctype][0, 0][0, 0]['gc_val']

    pseudo = 0.001
    #mask everything under 0
    replica1 = masked_less(replica1, 0) 
    replica2 = masked_less(replica2, 0)

    replica1 = norm_junctions(replica1, gc_factors=replica1_gc, gcnorm=True)
    replica2 = norm_junctions(replica2, gc_factors=replica2_gc, gcnorm=True)

    reads1 = replica1.mean(axis=1)+pseudo
    reads2 = replica2.mean(axis=1)+pseudo
    totalreads = reads1 + reads2 + pseudo * 2
    psi1 = reads1/totalreads
    psi2 = reads2/totalreads
    print psi1
    print psi2
    delta_psi = psi1 - psi2
    xlim(-1, 1)
    hist(delta_psi, bins = 60, histtype='step')

    if args.title:
        title(args.title)

    xlabel("Delta PSI")

    _save_or_show(args.plotpath, "deltapsi")



if __name__ == '__main__':
    main()

