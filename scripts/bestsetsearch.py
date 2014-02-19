import argparse
import pickle
from collections import defaultdict

from pylab import *
from numpy.ma import masked_less

from analysis.polyfitnb import fit_nb 
from analysis.filter import norm_junctions, discardlow, discardhigh, discardminreads, discardmaxreads, mark_stacks, discardminreads_and
from analysis.psi import calc_psi, mean_psi, simple_psi
from majiq_analysis import load_data
from analysis.sample import sample_from_junctions, mean_junction

"""
NOTE: This script should be copied one directory up to actually work. 
"""

def plot_densities(deltadata, my_title):
    deltadata = nan_to_num(deltadata) #substitute nan with zero, because histogram is a shitty function that cant take nans. Shame, shame on histogram. You should be a more manly function and take NaNs without crying, you are part of matplotlib.
    title(my_title)
    xlim(-1, 1)
    ylim(0, 2500)
    xlabel("Delta PSI")
    ylabel("Density")
    if len(deltadata[deltadata > 1]):
        print "DELTADATA BAD", deltadata[deltadata > 1]
        sys.exit(1)

    values, bins = histogram(deltadata, bins = 40)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    bar(center, values, align='center', width=width)
    print "VALUES", values
    print "EDGES", bins


def _save_or_show(plotpath, name):
    if plotpath:
        savefig("%s%s.png"%(plotpath, name), width=200, height=400, dpi=100)
        clf()
    else:
        show()  

def plot_save(deltas, plotpath, plot_title, plot_name):
    #check_deltas(deltas)
    print deltas.shape
    plot_densities(deltas, plot_title)
    _save_or_show(plotpath, plot_name)


def check_deltas(deltas):
    print "All deltas:", len(deltas)
    for i in arange(-0.9, 0, 0.1):
        print "Delta < %.1f:"%i, len(deltas[deltas < i])

    for i in arange(0.1, 1, 0.1):
        print "Delta > %.1f:"%i, len(deltas[deltas > i])

#best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"]
def doublezeroes(exc1, inc1, exc2, inc2):
    doubles1 = 0
    doubles2 = 0

    for i in xrange(len(inc1)):
        #print inc1[i].mean(), exc1[i].mean(), inc2[i].mean(), exc2[i].mean()
        if inc1[i].mean() == 0 and exc1[i].mean() == 0:
            doubles1 += 1

        if inc2[i].mean() == 0 and exc2[i].mean() == 0:
            doubles2 += 1
    print
    print "DOUBLE1", doubles1, "DOUBLE2", doubles2
    print 

def main():
    """
    Script to test what is happening with the empirical delta PSI that we intend to use as our prior
    """
    parser = argparse.ArgumentParser() 
    parser.add_argument('pair',  help='Path for samples of conditions 1 to N')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--debug', type=int,  default=0)
    args = parser.parse_args()

    #0 PREPROCESSING
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


    ######################################################
    """
    print "\n 1) Raw reads with raw PSI..."
    print all_junctions['inc1'].shape
    print all_junctions['inc2'].shape
    inc_reads1 = mean_junction(all_junctions['inc1'])
    inc_reads2 = mean_junction(all_junctions['inc2'])
    exc_reads1 = mean_junction(all_junctions['exc1'])
    exc_reads2 = mean_junction(all_junctions['exc2'])
    psi1 = simple_psi(inc_reads1, exc_reads1)
    psi2 = simple_psi(inc_reads2, exc_reads2)
    plot_save(array(psi1 - psi2), args.plotpath, "Raw reads, Raw PSI", "1_rawreads_rawpsi")

    #######################################################
    print '\n 2) Raw reads with dirichlet PSI ...'
    psi1 = calc_psi(all_junctions['inc1'], all_junctions['exc1'], name1, alpha=0.5, n=1, debug=args.debug, psinosample=True)
    psi2 = calc_psi(all_junctions['inc2'], all_junctions['exc2'], name2, alpha=0.5, n=1, debug=args.debug, psinosample=True) 
    plot_save(mean_psi(psi1) - mean_psi(psi2), args.plotpath, "Raw reads, Dirichlet PSI", "2_rawreads_psidirch")

    #######################################################
    print '\n 3) Bootstrapping reads with dirichlet PSI ...'
    mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(all_junctions["exc1"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)
    mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(all_junctions["inc1"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)      
    mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(all_junctions["exc2"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)
    mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(all_junctions["inc2"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)
    psi1 = calc_psi(inc_samples1, exc_samples1, name1, alpha=0.5, n=1, debug=args.debug, psinosample=True)
    psi2 = calc_psi(inc_samples2, exc_samples2, name2, alpha=0.5, n=1, debug=args.debug, psinosample=True)
    print "PSIBAD", psi1[psi1 > 1], psi1[psi1 < -1]
    plot_save(mean_psi(psi1) - mean_psi(psi2), args.plotpath, "Bootstrapping reads, Dirichlet PSI", "3_bootreads_psidirch")

    #######################################################
    print '\n 4)  ... filtering by minreads and minnonzero ...'
    best_set = defaultdict(array)
    best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardlow(10, True, False, all_junctions["exc1"], all_junctions["inc1"], all_junctions["exc2"], all_junctions["inc2"])
    best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads(50, True, False, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
    mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(best_set["exc1"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)
    mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(best_set["inc1"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)      
    mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(best_set["exc2"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)
    mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(best_set["inc2"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)
    psi1 = calc_psi(inc_samples1, exc_samples1, name1, alpha=0.5, n=1, debug=args.debug, psinosample=True)
    psi2 = calc_psi(inc_samples2, exc_samples2, name2, alpha=0.5, n=1, debug=args.debug, psinosample=True)
    print "PSIBAD", psi1[psi1 > 1], psi1[psi1 < -1]
    plot_save(mean_psi(psi1) - mean_psi(psi2), args.plotpath, "Bootstapping reads, Dirichlet PSI min10reads50", "4_bootreads_psidirch_min10reads50")    

    ########################################################
    print "\n 5) ... with AND filter ..."
    best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads(2, False, False, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
    mean_exc1, var_exc1, exc_samples1 = sample_from_junctions(best_set["exc1"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)
    mean_inc1, var_inc1, inc_samples1 = sample_from_junctions(best_set["inc1"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)      
    mean_exc2, var_exc2, exc_samples2 = sample_from_junctions(best_set["exc2"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)
    mean_inc2, var_inc2, inc_samples2 = sample_from_junctions(best_set["inc2"], m=10, k=5, discardzeros=True, trimborder=True, fitted_func=poly1d([1, 0]), debug=args.debug)    
    psi1 = calc_psi(inc_samples1, exc_samples1, name1, alpha=0.5, n=1, debug=args.debug, psinosample=True)
    psi2 = calc_psi(inc_samples2, exc_samples2, name2, alpha=0.5, n=1, debug=args.debug, psinosample=True)
    print "PSIBAD", psi1[psi1 > 1], psi1[psi1 < -1]
    plot_save(mean_psi(psi1) - mean_psi(psi2), args.plotpath, "Bootstapping reads, Dirichlet PSI min10reads50 AND some at least 2 reads", "5_bootreads_psidirch_min10reads50_AND") 


    ########################################################
    """
    for minandreads in xrange(1, 10):
        print "\n 6) Raw reads with raw PSI (filter ANDjunctions %s)..."%minandreads
        best_set = defaultdict(array)
        # The AND filter turned out to be more complex than I expected
        # minreads=0, orfilter=True, logger=False, returnempty=False

        best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads_and(incexcpairs=[[all_junctions["inc1"], all_junctions["exc1"]], [all_junctions["inc2"], all_junctions["exc2"]]], minreads=minandreads, logger=False)

        """
        mylogger = True
        orfilter = True
        returnempty = True
        best_set["exc1"], best_set["inc1"] = discardminreads(minandreads, orfilter, mylogger, returnempty, all_junctions["exc1"], all_junctions["inc1"]) 
        best_set["exc2"], best_set["inc2"] = discardminreads(minandreads, orfilter, mylogger, returnempty, all_junctions["exc2"], all_junctions["inc2"]) 
        orfilter = False
        returnempty = False
        best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads(-1, orfilter, mylogger, returnempty, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"]) 
        # end AND filter
        """
        inc_reads1 = mean_junction(best_set['inc1'])
        inc_reads2 = mean_junction(best_set['inc2'])
        exc_reads1 = mean_junction(best_set['exc1'])
        exc_reads2 = mean_junction(best_set['exc2'])

        doublezeroes(exc_reads1, inc_reads1, exc_reads2, inc_reads2)

        psi1 = simple_psi(inc_reads1, exc_reads1)
        psi2 = simple_psi(inc_reads2, exc_reads2)
        plot_save(array(psi1 - psi2), args.plotpath, "Raw reads, Raw PSI (filtering minANDreads %s)"%minandreads, "6_rawreads_rawpsi_filterand_%02d"%minandreads)
    

    best_ratio = 0
    minandreads = 1
    delta_significant = 0.2
    for minandreads in xrange(1, 5):
        for minpos in range(5, 30):
            for minorreads in range(5, 50):
                print "\n 7) Min pos %s Min OR reads %s"%(minpos, minorreads)
                best_set = defaultdict(array)

                #discardminreads_and()
                # The AND filter turned out to be more complex than I expected

                mylogger = False
                orfilter = True
                returnempty = True
                best_set["exc1"], best_set["inc1"] = discardminreads(minandreads, orfilter, mylogger, returnempty, all_junctions["exc1"], all_junctions["inc1"]) 
                best_set["exc2"], best_set["inc2"] = discardminreads(minandreads, orfilter, mylogger, returnempty, all_junctions["exc2"], all_junctions["inc2"]) 
                orfilter = False
                returnempty = False
                best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads(-1, orfilter, mylogger, returnempty, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"]) 
                # end AND filter
                orfilter = True
                best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardlow(minpos, orfilter, mylogger, best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])
                best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"] = discardminreads(minorreads, orfilter, mylogger, returnempty,best_set["exc1"], best_set["inc1"], best_set["exc2"], best_set["inc2"])

                inc_reads1 = mean_junction(best_set['inc1'])
                inc_reads2 = mean_junction(best_set['inc2'])
                exc_reads1 = mean_junction(best_set['exc1'])
                exc_reads2 = mean_junction(best_set['exc2'])

                doublezeroes(exc_reads1, inc_reads1, exc_reads2, inc_reads2)

                psi1 = simple_psi(inc_reads1, exc_reads1)
                psi2 = simple_psi(inc_reads2, exc_reads2)
                deltapsi = array(psi1 - psi2)
                print deltapsi
                ratio_significant = len(deltapsi[abs(deltapsi) > abs(delta_significant)]) / float(len(deltapsi)) 
                print "RATIO", ratio_significant

                #plot_save(array(psi1 - psi2), args.plotpath, "Raw reads, Raw PSI (filtering minANDreads %s minORreads %s minpos %s) "%(minandreads, minorreads, minpos), "7_rawreads_rawpsi_minor_%02d_minpos_%02d"%(minorreads, minpos))


if __name__ == '__main__':
    main()




