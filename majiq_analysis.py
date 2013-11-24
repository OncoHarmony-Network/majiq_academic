import sys
import os
import argparse
import pickle

from scipy.io import loadmat
from scipy.stats import nbinom
from pylab import *
from numpy.ma import masked_less

from analysis.polyfitnb import fit_nb 
from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads
from analysis.sample import sample_from_junctions
from analysis.psi import calc_psi


def _load_data3(my_mat, name, my_type):
    return {"junctions":my_mat[name][my_type][0, 0][0, 0]['cov'], "gc_content": my_mat[name][my_type][0, 0][0, 0]['gc_val']}


def _load_data2(my_mat, name):
    inc = _load_data3(my_mat, name, 'Inc')
    exc = _load_data3(my_mat, name, 'Exc')
    const = _load_data3(my_mat, name, 'rand10k')
    return inc, exc, const

def load_data(path, name1, name2):
    "Load data from a mat file. Should be deprecated soon"
    my_mat = loadmat(path)
    inc1, exc1, const1 = _load_data2(my_mat, name1)
    inc2, exc2, const2 = _load_data2(my_mat, name2)
    return inc1, exc1, const1, inc2, exc2, const2

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help='The matlab files pairs that we read. Can be used as a glob.')
    parser.add_argument('--trim', default=0, type=int, help='Trim the borders of the junctions because of poor mappability')
    parser.add_argument('--minnonzero', default=0, type=int, help='Minimum number of positions ')
    parser.add_argument('--skipgc', dest="norm", action='store_true', default=True, help='Disable GC content')
    parser.add_argument('--discardb', action='store_true',  default=False, help='Discard the b from the polynomial, since we expect our fit to start from 0, 0')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--tmp', required=True, help='Path to save the temporary files.')
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    parser.add_argument('--debug', type=int,  default=0)
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration')
    parser.add_argument('--m', default=100, type=int, help='Number of bootstrapping samples')  
    parser.add_argument('--trimborder', default=5, type=int, help='Trim the borders when sampling (keeping the ones with reads)')
    parser.add_argument('--maxnonzero', default=0, type=int, help='Maximum number of positive positions to consider the junction') 
    parser.add_argument('--minreads', default=0, type=int, help='Minimum number of reads combining all positions in a junction to be considered') 
    parser.add_argument('--maxreads', default=0, type=int, help='Maximum number of reads combining all positions in a junction to be considered')    
    parser.add_argument('--names', nargs='+', required=True, help="The names that identify each of the experiments.")
    parser.add_argument('--n', default=1, type=int, help='Number of PSI samples per sample paired') 
    parser.add_argument('--alpha', default=0.5, type=float, help='Alpha hyperparameter for the dirichlet distribution') 
    parser.add_argument('--nodiscardzeros', default=True, dest="discardzeros", action='store_false', help='Skip discarding zeroes')
    parser.add_argument('--nbdisp', default=0.01, type=int, help='Dispersion for the fallback NB function')
    parser.add_argument('--orfilter', default=False, action='store_true', help='When filtering, select sets of junctions where at least one passes the filter, instead of all passing the filter.')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        print "Creating directory %s..."%args.output
        os.makedirs(args.output)

    if args.plotpath:
        if not os.path.exists(args.plotpath):
            print "Creating directory %s..."%args.plotpath
            os.makedirs(args.plotpath)

    for path in args.files:
        print "\n\n\nProcessing %s..."%path
        inc1, exc1, const1, inc2, exc2, const2 = load_data(path, args.names[0], args.names[1]) #loading the paired matrixes
        
        print "GC content normalization..."
        inc1 = norm_junctions(inc1["junctions"], inc1["gc_content"])
        exc1 = norm_junctions(exc1["junctions"], exc1["gc_content"])
        const1 = norm_junctions(const1["junctions"], const1["gc_content"])
        inc2 = norm_junctions(inc2["junctions"], inc2["gc_content"])
        exc2 = norm_junctions(exc2["junctions"], exc2["gc_content"])
        const2 = norm_junctions(const2["junctions"], const2["gc_content"])

        print "Masking non unique and stats (n < 0)..."
        inc1 = masked_less(inc1, 0) 
        inc2 = masked_less(inc2, 0)
        exc1 = masked_less(exc1, 0) 
        exc2 = masked_less(exc2, 0)
        const1 = masked_less(const1, 0) 
        const2 = masked_less(const2, 0)

        if args.debug:
            print "Skipping fitfunc because --debug!!"
            fitfunc1 = poly1d([1, 0])
            fitfunc2 = poly1d([1, 0])
        else:
            print "Fitting NB function with constitutive events..."
            fitfunc1 = fit_nb(const1, "%s/const1"%args.output, args.plotpath)
            fitfunc2 = fit_nb(const2, "%s/const2"%args.output, args.plotpath)

        if args.maxnonzero or args.minnonzero or args.minreads or args.maxreads:
            print "Filtering..."

        if args.maxnonzero:
            exc1, inc1, exc2, inc2 = discardhigh(args.maxnonzero, args.orfilter, args.debug, exc1, inc1, exc2, inc2)
            
        if args.minnonzero:
            exc1, inc1, exc2, inc2 = discardlow(args.minnonzero, args.orfilter, args.debug, exc1, inc1, exc2, inc2)

        if args.minreads:
            exc1, inc1, exc2, inc2 = discardminreads(args.minreads, args.orfilter, args.debug, exc1, inc1, exc2, inc2)

        if args.maxreads:
            exc1, inc1, exc2, inc2 = discardmaxreads(args.maxreads, args.orfilter, args.debug, exc1, inc1, exc2, inc2)
              
        print exc1.shape
        print inc1.shape

        mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(exc1, args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc1, debug=args.debug)
        mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(inc1, args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc1, debug=args.debug)      
        
        mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(exc2, args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc2, debug=args.debug)
        mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(inc2, args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc2, debug=args.debug)
        
        pickle.dump(exc_samples1, open("%s_%s_exc_samples.pickle"%(args.output, args.names[0]), 'w'))
        pickle.dump(inc_samples1, open("%s_%s_inc_samples.pickle"%(args.output, args.names[1]), 'w'))
        pickle.dump(exc_samples2, open("%s_%s_exc_samples.pickle"%(args.output, args.names[0]), 'w'))
        pickle.dump(inc_samples1, open("%s_%s_inc_samples.pickle"%(args.output, args.names[1]), 'w'))
        
        print "\nCalculating PSI for %s"%(args.names)
        samples1 = vstack([inc_samples1, exc_samples1]).reshape(2, inc_samples1.shape[0], inc_samples1.shape[1])
        samples2 = vstack([inc_samples2, exc_samples2]).reshape(2, inc_samples2.shape[0], inc_samples2.shape[1])
        psi_scores1 = calc_psi(args.alpha, args.n, samples1, debug=args.debug)
        pickle.dump(psi_scores1, open("%s_%s_inc_vs_exc_psivalues.pickle"%(args.output, args.names[0]), 'w'))
        psi_scores2 = calc_psi(args.alpha, args.n, samples2, debug=args.debug)
        pickle.dump(psi_scores2, open("%s_%s_inc_vs_exc_psivalues.pickle"%(args.output, args.names[1]), 'w'))




if __name__ == '__main__':
    main()



