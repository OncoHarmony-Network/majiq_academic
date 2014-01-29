import argparse
import pickle

"We want to know how many events have a significant PSI of 1 or 0 in one side, and no reads at all in the other"


from collections import defaultdict

from pylab import *
from numpy.ma import masked_less

from analysis.polyfitnb import fit_nb 
from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, mark_stacks
from analysis.psi import calc_psi, mean_psi
from majiq_analysis import load_data
from analysis.sample import sample_from_junctions

def plot_densities(deltadata, my_title):

    deltadata = nan_to_num(deltadata) #substitute nan with zero, because histogram is a 
    title(my_title)
    xlim(-1, 1)
    xlabel("Delta PSI")
    ylabel("Density")
    values, edges = histogram(deltadata, bins = 40)      
    plot(linspace(-1, 1, num=len(values)), values)    


def _save_or_show(plotpath, name):
    if plotpath:
        savefig("%s%s.png"%(plotpath, name), width=200, height=400, dpi=100)
        clf()
    else:
        show()  

def plot_save(deltas, plotpath, plot_title, plot_name):
    check_deltas(deltas)
    plot_densities(deltas, plot_title)
    _save_or_show(plotpath, plot_name)


def check_deltas(deltas):
    print "All deltas:", len(deltas)
    for i in arange(-0.9, 0, 0.1):
        print "Delta < %.1f:"%i, len(deltas[deltas < i])

    for i in arange(0.1, 1, 0.1):
        print "Delta > %.1f:"%i, len(deltas[deltas > i])


def main():
    """
    Script to test what is happening with the empirical delta PSI that we intend to use as our prior
    """
    parser = argparse.ArgumentParser() 
    parser.add_argument('pair',  help='Path for samples of conditions 1 to N')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--debug', type=int,  default=0)
    parser.add_argument('--minreads', type=int,  default=20, help='')   
    args = parser.parse_args()
    print "Loading all..."
    filename = args.pair.split('/')[-1]
    name1, name2 = filename.split('.')[0].split('_')
    all_junctions = load_data(args.pair, name1, name2)

    print "Preprocessing..."
    inc1, exc1, const1, inc2, exc2, const2 = load_data(args.pair, name1, name2) #loading the paired matrixes
    all_junctions = {"inc1": inc1, "exc1": exc1, "inc2": inc2, "exc2": exc2 }

    for junc_set in all_junctions.keys():
        all_junctions[junc_set] = norm_junctions(all_junctions[junc_set]["junctions"], all_junctions[junc_set]["gc_content"])

    for junc_set in all_junctions.keys():
        all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

    print '\nRaw reads with dirichlet PSI ...'
    pseudo = 0.000000001
    inc_reads1 = all_junctions['inc1'].mean(axis=1)+pseudo
    inc_reads2 = all_junctions['inc2'].mean(axis=1)+pseudo
    exc_reads1 = all_junctions['exc1'].mean(axis=1)+pseudo
    exc_reads2 = all_junctions['exc2'].mean(axis=1)+pseudo
    psi1 = inc_reads1/(exc_reads1+inc_reads1)
    psi2 = inc_reads2/(exc_reads2+inc_reads2)
    for i in range(inc_reads1):
    	if inc_reads1 >  exc_reads1



if __name__ == '__main__':
    main()




