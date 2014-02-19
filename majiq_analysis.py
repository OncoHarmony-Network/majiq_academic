import sys
import os
import argparse
from collections import defaultdict
try:
    import cPickle as pickle
except:
    import pickle

from scipy.io import loadmat
from scipy.stats import nbinom
from pylab import *
from numpy.ma import masked_less
from scipy.io import savemat

from analysis.polyfitnb import fit_nb 
from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, mark_stacks
from analysis.sample import sample_from_junctions
from analysis.psi import calc_psi, mean_psi, DirichletCalc, reads_given_psi, BINS_CENTER
from analysis.adjustdelta import adjustdelta
from grimoire.utils.utils import get_logger

def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname), bbox_inches='tight') 
        clf()
    else:
        show()

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

def _create_if_not_exists(my_dir):
    if not os.path.exists(my_dir):
        print "\nCreating directory %s..."%my_dir
        os.makedirs(my_dir)    

def __debug_stacks(args, all_junctions):
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
    return all_junctions


def __write_samples(sample, output, names, num):
    pickle.dump(sample, open("%s%s_exc%s_samples.pickle"%(output, names[num], num), 'w'))





def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', help='The matlab files pairs that we read. Can be used as a glob.')
    #required flags
    parser.add_argument('--tmp', required=True, help='Path to save the temporary files.')
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    #optional flags (everything below this)
    parser.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    parser.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')
    parser.add_argument('--trim', default=0, type=int, help='Trim the borders of the junctions because of poor mappability')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration. [Default: %(default)s]')
    parser.add_argument('--m', default=100, type=int, help='Number of bootstrapping samples. [Default: %(default)s]')  
    parser.add_argument('--trimborder', default=5, type=int, help='Trim the borders when sampling (keeping the ones with reads). [Default: %(default)s]')
    parser.add_argument('--names', nargs='+', required=True, help="The names that identify each of the experiments. [Default: %(default)s]")
    parser.add_argument('--n', default=1, type=int, help='Number of PSI samples per sample paired. [Default: %(default)s]') 
    parser.add_argument('--alpha', default=0.5, type=float, help='Alpha hyperparameter for the dirichlet distribution. [Default: %(default)s]') 
    parser.add_argument('--markstacks', default=0.001, type=float, help='Mark stack positions. Expects a p-value. Use a negative value in order to disable it. [Default: %(default)s]') 
    parser.add_argument('--nbdisp', default=0.1, type=int, help='Dispersion for the fallback NB function. [Default: %(default)s]')
    parser.add_argument('--psiparam', default=False, action='store_true', help='Instead of sampling, use a parametric form for the PSI calculation. [Default: %(default)s]')
    #parser.add_argument('--orfilter', default=False, action='store_true', help='When filtering, select sets of junctions where at least one passes the filter, instead of all passing the filter. [Default: %(default)s]')
    #EM flag
    parser.add_argument('--binsize', default=0.025, type=int, help='The bins for PSI values. With a --binsize of 0.025 (default), we have 40 bins')   
    parser.add_argument('--minreads', default=50, type=int, help='Minimum number of reads combining all positions in a junction to be considered. [Default: %(default)s]') 
    parser.add_argument('--minnonzero', default=10, type=int, help='Minimum number of positions for the best set.')
    parser.add_argument('--iter', default=10, type=int, help='Max number of iterations of the EM')
    parser.add_argument('--breakiter', default=0.01, type=float, help='If the log likelihood increases less that this flag, do not do another EM step')
    parser.add_argument('--V', default=0.1, type=float, help='Value of DeltaPSI used for initialization of the EM model [Default: %(default)s]')
    #Disable flags
    parser.add_argument('--nogc', dest="norm", action='store_true', default=True, help='Disable GC content normalization')
    parser.add_argument('--nodiscardb', dest="discardb", action='store_false',  default=True, help='Skip biscarding the b from the NB polynomial function, since we expect our fit to start from x=0, y=0')
    parser.add_argument('--nodiscardzeros', default=True, dest="discardzeros", action='store_false', help='Skip discarding zeroes')
    #Debug flags  
    parser.add_argument('--ONLYSTACKS', action='store_true', help="Useless flag for analysis. Used to test if stacks are worth masking.")
    parser.add_argument('--usetensor', action='store_true')
    parser.add_argument('--debug', type=int, default=0, help="Activate this flag for debugging purposes, activates logger and jumps some processing steps.")

    args = parser.parse_args()
    #create directories if they dont exist
    _create_if_not_exists(args.output)
    _create_if_not_exists(args.plotpath)
    #grab logger
    logger_path = args.logger
    if not logger_path:
        logger_path = args.output

    logger = get_logger("%smajiq.log"%logger_path, silent=args.silent, debug=args.debug)

    for path in args.files:
        logger.info("")
        logger.info("Processing %s..."%path)
        inc1, exc1, const1, inc2, exc2, const2 = load_data(path, args.names[0], args.names[1]) #loading the paired matrixes
        all_junctions = {"inc1": inc1, "exc1": exc1, "const1": const1, 
                         "inc2": inc2, "exc2": exc2, "const2": const2 }

        logger.info("GC content normalization...")
        for junc_set in all_junctions.keys():
            all_junctions[junc_set] = norm_junctions(all_junctions[junc_set]["junctions"], 
                                                     all_junctions[junc_set]["gc_content"])

        logger.info("Masking non unique...")
        for junc_set in all_junctions.keys():
            all_junctions[junc_set] = masked_less(all_junctions[junc_set], 0) 

        if args.debug:
            logger.debug("Skipping fitfunc because --debug!!")
            fitfunc1 = poly1d([1, 0])
            fitfunc2 = poly1d([1, 0])
        else:
            logger.info("Fitting NB function with constitutive events...")
            fitfunc1 = fit_nb(all_junctions["const1"], "%s_nbfit1"%args.output, args.plotpath, nbdisp=args.nbdisp)
            fitfunc2 = fit_nb(all_junctions["const2"], "%s_nbfit2"%args.output, args.plotpath, nbdisp=args.nbdisp)

        if args.markstacks >= 0:
            logger.info("Marking and masking stacks...")
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

            if args.ONLYSTACKS: #Just for analysis and debugging of stacks
                __debug_stacks(args, all_junctions)


        #TODO This should be happening for every pair of N/M experiments (for N in exp: for M in exp:)
        logger.info('Filtering to obtain "best set"...')
        best_set = defaultdict(array)
        best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardlow(args.minnonzero, True, logger, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"])
        best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads(args.minreads, True, logger, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
         
        logger.info("Bootstrapping for 'best set'...")
        mean_exc1, var_exc1, exc_best_samples1 = sample_from_junctions(best_set["exc1"], args.m, args.k, discardzeros=False, trimborder=args.trimborder, fitted_func=fitfunc1, debug=args.debug)
        mean_inc1, var_inc1, inc_best_samples1 = sample_from_junctions(best_set["inc1"], args.m, args.k, discardzeros=False, trimborder=args.trimborder, fitted_func=fitfunc1, debug=args.debug)      
        mean_exc2, var_exc2, exc_best_samples2 = sample_from_junctions(best_set["exc2"], args.m, args.k, discardzeros=False, trimborder=args.trimborder, fitted_func=fitfunc2, debug=args.debug)
        mean_inc2, var_inc2, inc_best_samples2 = sample_from_junctions(best_set["inc2"], args.m, args.k, discardzeros=False, trimborder=args.trimborder, fitted_func=fitfunc2, debug=args.debug)

        logger.info("\nCalculating PSI for 'best set' %s ..."%(args.names))
        best_psi1 = calc_psi(inc_best_samples1, exc_best_samples1, args.names[0], args.alpha, args.n, args.debug, args.psiparam)
        best_psi2 = calc_psi(inc_best_samples2, exc_best_samples2, args.names[1], args.alpha, args.n, args.debug, args.psiparam)
        
        logger.info("\nCalculating delta PSI for 'best set' %s ..."%(args.names))
        best_delta_psi = mean_psi(best_psi1) - mean_psi(best_psi2)

        logger.info("Obtaning prior matrix for 'best set'...")
        mixture_pdf = adjustdelta(best_delta_psi, args.output, plotpath=args.plotpath, title=" ".join(args.names), numiter=args.iter, breakiter=args.breakiter, V=args.V)

        logger.info("Calculating prior matrix...")
        numbins = int(round(len(mixture_pdf)/2)) #half the delta bins
        dircalc = DirichletCalc() 
        #Calculate prior matrix
        prior_matrix = []
        for i in xrange(numbins):
            prior_matrix.extend(mixture_pdf[numbins-i:(numbins*2)-i])

        prior_matrix = array(prior_matrix).reshape(numbins, -1)

        clf()
        title("Prior Matrix")
        imshow(prior_matrix)
        xlabel("PSI_i")
        ylabel("PSI_j")
        _save_or_show(args.plotpath, plotname="prior_matrix_no_jefferies")

        #Adjust prior matrix with Jefferies prior        
        jefferies = []
        psi_space = linspace(0, 1-args.binsize, num=numbins) + args.binsize/2
        for i in psi_space:
            jefferies.append([])
            for j in psi_space:
                jefferies[-1].append(dircalc.pdf([i, 1-i, j, 1-j], [args.alpha, args.alpha, args.alpha, args.alpha]))

        #jefferies = array([dircalc.pdf([x, 1-x], [0.5, 0.5]) for x in psi_space])
        jefferies = array(jefferies)
        jefferies /= sum(jefferies)

        clf()
        title("Jefferies Matrix")
        imshow(jefferies)
        xlabel("PSI_i")
        ylabel("PSI_j")
        _save_or_show(args.plotpath, plotname="jefferies_matrix")

        prior_matrix *= jefferies #Normalize PSI with jefferies
        prior_matrix /= sum(prior_matrix) #renormalize so it sums 1

        clf()
        title("Prior Matrix")
        imshow(prior_matrix)
        xlabel("PSI 1")
        ylabel("PSI 2")
        _save_or_show(args.plotpath, plotname="prior_matrix")

        logger.info("Saving prior matrix for %s..."%(args.names))
        pickle.dump(prior_matrix, open("%s%s_%s_priormatrix.pickle"%(args.output, args.names[0], args.names[1]), 'w'))

        logger.info("Bootstrapping for all samples...")
        mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(all_junctions["exc1"], args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc1, debug=args.debug)
        mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(all_junctions["inc1"], args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc1, debug=args.debug)      
        mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(all_junctions["exc2"], args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc2, debug=args.debug)
        mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(all_junctions["inc2"], args.m, args.k, discardzeros=args.discardzeros, trimborder=args.trimborder, fitted_func=fitfunc2, debug=args.debug)
  
        logger.info("Writing samples in disk...")
        __write_samples(exc_samples1, args.output, args.names, 0)
        __write_samples(inc_samples1, args.output, args.names, 1)
        __write_samples(exc_samples2, args.output, args.names, 0)
        __write_samples(inc_samples2, args.output, args.names, 1)

        logger.info("\nCalculating PSI %s ..."%(args.names))
        psi1 = calc_psi(inc_samples1, exc_samples1, args.names[0], args.alpha, args.n, args.debug, args.psiparam)
        psi2 = calc_psi(inc_samples2, exc_samples2, args.names[1], args.alpha, args.n, args.debug, args.psiparam)

        pickle.dump(psi1, open("%s%s_psi.pickle"%(args.output, args.names[0]), 'w'))
        pickle.dump(psi1, open("%s%s_psi.pickle"%(args.output, args.names[1]), 'w'))

        logger.info("Calculating P(Data | PSI_i, PSI_j)...")
        #P(Data | PSI_i, PSI_j) = P(vector_i | PSI_i) * P(vector_j | PSI_j)
        data_given_psi1 = reads_given_psi(inc_samples1, exc_samples1, psi_space)
        data_given_psi2 = reads_given_psi(inc_samples2, exc_samples2, psi_space)


        data_given_psi = []
        for sample in xrange(data_given_psi1.shape[0]):
            #TODO Tensor product is calculated with scipy.stats.kron. Probably faster, have to make sure I am using it correctly.
            data_given_psi.append(data_given_psi1[sample].reshape(-1, numbins) * data_given_psi2[sample].reshape(numbins, -1))
            print "Event %s"%sample
            print inc_samples1[sample], exc_samples1[sample]
            print
            print sum(data_given_psi[sample])
            clf()
            title("P(Data | PSI 1, PSI 2) Event %s (Inc1: %s, Exc1: %s Inc2: %s Exc2: %s)"%(sample, sum(inc_samples1[sample]), sum(exc_samples1[sample]), sum(inc_samples2[sample]), sum(exc_samples2[sample])))
            imshow(data_given_psi[sample])
            xlabel("PSI 1")
            ylabel("PSI 2")
            _save_or_show(args.plotpath, plotname="datagpsi_sample%s"%sample)


        #Finally, P(PSI_i, PSI_j | Data) proportional P(PSI_i, PSI_j)* P(Data | PSI_i, PSI_j) 
        logger.info("Calculate Posterior Delta Matrices...")
        posterior_matrix = []
        for sample in xrange(len(data_given_psi)):
            pm = (prior_matrix * data_given_psi[sample])
            posterior_matrix.append(pm / sum(pm))

            clf()
            title("Posterior Delta Event %s ()"%sample)
            imshow(posterior_matrix[-1])
            xlabel("PSI i")
            ylabel("PSI j")
            _save_or_show(args.plotpath, plotname="postdelta_sample%s"%sample)

        pickle.dump(posterior_matrix, open("%s%s_%s_deltamatrix.pickle"%(args.output, args.names[0], args.names[1]), 'w'))
        logger.info("Done!")

if __name__ == '__main__':
    main()



