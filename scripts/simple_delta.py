import sys
import os
import argparse
import pickle

from pylab import *
from numpy.random import choice
from scipy.io import loadmat
from numpy.ma import masked_less
 
from analysis.filter import discardhigh, discardlow, discardminreads, discardmaxreads, norm_junctions, mark_stacks

BINS = linspace(0, 1, num=99)

"""
NOTE: This script should be copied one directory up to actually work. 

"""

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
    parser.add_argument('--markstacks', default=0.001, type=float, help='Mark stack positions. Expects a p-value. Use a negative value in order to disable it. [Default: %(default)s]')  
    parser.add_argument('--orfilter', default=False, action='store_true', help='When filtering, select sets of junctions where at least one passes the filter, instead of all passing the filter.')
    parser.add_argument('--output', default=None, help="Path to save the results to.")
    parser.add_argument('--ONLYSTACKS', action='store_true', help="Useless flag for analysis. Used to test if stacks are worth masking.")
    parser.add_argument('--names', nargs='+', required=True, help="The names that identify each of the experiments. [Default: %(default)s]")    
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

    if args.markstacks >= 0:
        print "Marking and masking stacks..."
        inclusion1 = mark_stacks(inclusion1, poly1d([1, 0]), args.markstacks, 0.01)
        inclusion2 = mark_stacks(inclusion2, poly1d([1, 0]), args.markstacks, 0.01)
        exclusion1 = mark_stacks(exclusion1, poly1d([1, 0]), args.markstacks, 0.01)
        exclusion2 = mark_stacks(exclusion2, poly1d([1, 0]), args.markstacks, 0.01)
        if not args.ONLYSTACKS:
            inclusion1 = masked_less(inclusion1, 0) #remask the stacks
            exclusion2 = masked_less(exclusion2, 0) #remask the stacks
            exclusion1 = masked_less(exclusion1, 0) #remask the stacks
            inclusion2 = masked_less(inclusion2, 0) #remask the stacks

    #START Just for analysis, should go into a script
    if args.ONLYSTACKS:
        print "Before", inclusion1.shape
        only_stacks_inc1 = []
        only_stacks_exc1 = []
        only_stacks_inc2 = []
        only_stacks_exc2 = []        
        stacks_ids = [] #the ones we know are being filtered    
        for i in range(inclusion1.shape[0]): 
            if ((inclusion1[i] == -2).sum() > 0) or ((exclusion1[i] == -2).sum() > 0) or (inclusion2[i] == -2).sum() > 0 or (exclusion2[i] == -2).sum() > 0:
                only_stacks_inc1.append(inclusion1[i])
                only_stacks_exc1.append(exclusion1[i])
                only_stacks_inc2.append(inclusion2[i])
                only_stacks_exc2.append(exclusion2[i])    
                stacks_ids.append(i)

        inclusion1 = array(only_stacks_inc1)
        inclusion1 = masked_less(inclusion1, 0) #remask the stacks
        exclusion1 = array(only_stacks_exc1)
        exclusion1 = masked_less(exclusion1, 0) #remask the stacks
        inclusion2 = array(only_stacks_inc2)
        inclusion2 = masked_less(inclusion2, 0) #remask the stacks            
        exclusion2 = array(only_stacks_exc2)
        exclusion2 = masked_less(exclusion2, 0) #remask the stacks
        print "Saving indexes..."
        pickle.dump(stacks_ids, open("%s_%s_vs_%s_indexesSTACKS.pickle"%(args.output, args.names[0], args.names[1]), 'w'))
        print "After", inclusion1.shape
    #END Just for analysis




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
    pickle.dump(delta_psi, open("%s_deltapsi.pickle"%(args.output), 'w'))


if __name__ == '__main__':
    main()

