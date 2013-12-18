import sys
import os
import argparse
try:
    import cPickle as pickle
except:
    import pickle

from scipy.io import loadmat
from scipy.stats import nbinom
from pylab import *
from numpy.ma import masked_less

from analysis.polyfitnb import fit_nb 
from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, mark_stacks
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
    "Load data from a matlab file. Should be deprecated soon"
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
    parser.add_argument('--nodiscardb', dest="discardb", action='store_false',  default=True, help='Skip biscarding the b from the NB polynomial function, since we expect our fit to start from x=0, y=0')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--tmp', required=True, help='Path to save the temporary files.')
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    parser.add_argument('--debug', type=int,  default=0)
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration. [Default: %(default)s]')
    parser.add_argument('--m', default=100, type=int, help='Number of bootstrapping samples. [Default: %(default)s]')  
    parser.add_argument('--trimborder', default=5, type=int, help='Trim the borders when sampling (keeping the ones with reads). [Default: %(default)s]')
    parser.add_argument('--maxnonzero', default=0, type=int, help='Maximum number of positive positions to consider the junction.') 
    parser.add_argument('--minreads', default=0, type=int, help='Minimum number of reads combining all positions in a junction to be considered. [Default: %(default)s]') 
    parser.add_argument('--maxreads', default=0, type=int, help='Maximum number of reads combining all positions in a junction to be considered. [Default: %(default)s]')    
    parser.add_argument('--names', nargs='+', required=True, help="The names that identify each of the experiments. [Default: %(default)s]")
    parser.add_argument('--n', default=1, type=int, help='Number of PSI samples per sample paired. [Default: %(default)s]') 
    parser.add_argument('--alpha', default=0.5, type=float, help='Alpha hyperparameter for the dirichlet distribution. [Default: %(default)s]') 
    parser.add_argument('--markstacks', default=0.001, type=float, help='Mark stack positions. Expects a p-value. Use a negative value in order to disable it. [Default: %(default)s]') 
    parser.add_argument('--nodiscardzeros', default=True, dest="discardzeros", action='store_false', help='Skip discarding zeroes')
    parser.add_argument('--nbdisp', default=0.1, type=int, help='Dispersion for the fallback NB function. [Default: %(default)s]')
    parser.add_argument('--orfilter', default=False, action='store_true', help='When filtering, select sets of junctions where at least one passes the filter, instead of all passing the filter. [Default: %(default)s]')
    parser.add_argument('--ONLYSTACKS', action='store_true', help="Useless flag for analysis. Used to test if stacks are worth masking.")
    args = parser.parse_args()

    #create directories if they dont exist
    if not os.path.exists(args.output):
        print "\nCreating directory %s..."%args.output
        os.makedirs(args.output)

    if args.plotpath:
        if not os.path.exists(args.plotpath):
            print "\nCreating directory %s..."%args.plotpath
            os.makedirs(args.plotpath)

    for path in args.files:
        print "\nProcessing %s..."%path
        inc1, exc1, const1, inc2, exc2, const2 = load_data(path, args.names[0], args.names[1]) #loading the paired matrixes
        all_junctions = {"inc1": inc1, "exc1": exc1, "const1": const1, 
                         "inc2": inc2, "exc2": exc2, "const2": const2 } #TODO Generalize to N junctions

        print "GC content normalization..."
        for junc_set in all_junctions.keys():
            all_junctions[junc_set] = norm_junctions(all_junctions[junc_set]["junctions"], 
                                                     all_junctions[junc_set]["gc_content"])

        print "Masking non unique..."
        for junc_set in all_junctions.keys():
            all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

        if args.debug:
            print "Skipping fitfunc because --debug!!"
            fitfunc1 = poly1d([1, 0])
            fitfunc2 = poly1d([1, 0])
        else:
            print "Fitting NB function with constitutive events..."
            fitfunc1 = fit_nb(all_junctions["const1"], "%s_nbfit1"%args.output, args.plotpath, nbdisp=args.nbdisp)
            fitfunc2 = fit_nb(all_junctions["const2"], "%s_nbfit2"%args.output, args.plotpath, nbdisp=args.nbdisp)

        if args.markstacks >= 0:
            print "Marking and masking stacks..."
            for junc_set in all_junctions.keys():
                if junc_set.find("const") == -1:
                    if junc_set.endswith("1"):
                        f = fitfunc1
                    else:
                        f = fitfunc2

                    print "... %s"%junc_set
                    all_junctions[junc_set] = mark_stacks(all_junctions[junc_set], f, args.markstacks, args.nbdisp)
                    if not args.ONLYSTACKS:
                        all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) #remask the stacks

        #START Just for analysis, should go into a script
        if args.ONLYSTACKS:
            print "Before", all_junctions["inc1"].shape
            only_stacks_inc1 = []
            only_stacks_exc1 = []
            only_stacks_inc2 = []
            only_stacks_exc2 = []        
            stacks_ids = [] #the ones we know are being filtered    
            for i in range(all_junctions["inc1"].shape[0]): 
                if ((all_junctions["inc1"][i] == -2).sum() > 0) or ((all_junctions["exc1"][i] == -2).sum() > 0) or (all_junctions["inc2"][i] == -2).sum() > 0 or (all_junctions["exc2"][i] == -2).sum() > 0:
                    only_stacks_inc1.append(all_junctions["inc1"][i])
                    only_stacks_exc1.append(all_junctions["exc1"][i])
                    only_stacks_inc2.append(all_junctions["inc2"][i])
                    only_stacks_exc2.append(all_junctions["exc2"][i])    
                    stacks_ids.append(i)

            all_junctions["inc1"] = array(only_stacks_inc1)
            all_junctions["inc1"] = masked_less(all_junctions["inc1"], 0) #remask the stacks
            all_junctions["exc1"] = array(only_stacks_exc1)
            all_junctions["exc1"] = masked_less(all_junctions["exc1"], 0) #remask the stacks
            all_junctions["inc2"] = array(only_stacks_inc2)
            all_junctions["inc2"] = masked_less(all_junctions["inc2"], 0) #remask the stacks            
            all_junctions["exc2"] = array(only_stacks_exc2)
            all_junctions["exc2"] = masked_less(all_junctions["exc2"], 0) #remask the stacks
            print "Saving indexes..."
            pickle.dump(stacks_ids, open("%s_%s_vs_%s_indexesSTACKS.pickle"%(args.output, args.names[0], args.names[1]), 'w'))
            print "After", all_junctions["inc1"].shape
        #END Just for analysis


        if args.maxnonzero or args.minnonzero or args.minreads or args.maxreads:
            print "Filtering..."

        if args.maxnonzero:
            all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"] = discardhigh(args.maxnonzero, args.orfilter, args.debug, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"])
            
        if args.minnonzero:
            all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"] = discardlow(args.minnonzero, args.orfilter, args.debug, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"])

        if args.minreads:
            all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"] = discardminreads(args.minreads, args.orfilter, args.debug, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"])

        if args.maxreads:
            all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"] = discardmaxreads(args.maxreads, args.orfilter, args.debug, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"])

        mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(all_junctions["exc1"], args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc1, debug=args.debug)
        mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(all_junctions["inc1"], args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc1, debug=args.debug)      
        
        mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(all_junctions["exc2"], args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc2, debug=args.debug)
        mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(all_junctions["inc2"], args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc2, debug=args.debug)
        
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



