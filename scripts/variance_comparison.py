from matplotlib import use
use('Agg', warn=False)
import sys
import os
from analysis import polyfitnb
import argparse

from scipy.io import loadmat
from pylab import *
from numpy.ma import masked_less
from matplotlib import rcParams

# from junction_sample import sample_from_junctions, norm_junctions, calc_nonzeromeanvar, discardlow, discardhigh, DEBUG, TESTBREAK, EPSILON
import analysis.filter
import polyfit_ectf
import pipelines
import junction_sample

def _numzeros(junctions):
    "Obtain the number of zeros for every junction"
    return (junctions == 0).sum(axis=1)

def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname.replace(" ", "_")), bbox_inches='tight', width=1000, height=2000, dpi=300) #WNo spaces allowed, underscores!
        clf()
    else:
        show()

def calc_score_tmp(mean_sample1, var_sample1, mean_sample2, var_sample2):
    print abs(mean_sample1-mean_sample2)/((mean_sample1+mean_sample2)*0.5)
    return {'mean': abs(mean_sample1-mean_sample2)/((mean_sample1+mean_sample2)*0.5),
            'variance': abs(var_sample1-var_sample2)/((var_sample1+var_sample2)*0.5)}

def calc_score(mean_sample1, var_sample1, mean_sample2, var_sample2):
    return calc_score_tmp(mean_sample1, var_sample1, mean_sample2, var_sample2)
    # return abs(var_sample1-var_sample2)/((mean_sample1+mean_sample2)*0.5)

def plot_method1Vsmethod2(score1, score2, score1_name, score2_name, replica1_name, replica2_name, plotpath=None):
    """Compute 2 plots: one for the variance, one for the mean"""
    plotname="%sVs%s" % (score1_name, score2_name)
    total_junctions = float(len(score1))
    better_in_method1 = sum(score1 < score2)
    equal_in_both = sum(score1 == score2)
    better_in_method2 = sum(score1 > score2)

    print "Better in %s: %s (%.2f%%) Equal: %s Better in %s: %s (%.2f%%)"%(score1_name, better_in_method1, (better_in_method1/total_junctions)*100, equal_in_both, score2_name, better_in_method2, (better_in_method2/total_junctions)*100)
    title("%s vs %s\n%s and %s"%(score1_name, score2_name, replica1_name, replica2_name))
    max_value = max(max(score1), max(score2))
    xlabel(score1_name)
    ylabel(score2_name)
    xlim(0, max_value)
    ylim(0, max_value)
    plot([0, max_value], [0, max_value])
    plot(score1, score2, '.')
    _save_or_show(plotpath, plotname)


def compare_methods(m, k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, method1_name, method2_name, plotpath, scores_cached):

    if method1_name in scores_cached:
        score_method1 = scores_cached[method1_name]
    else:
        mean_method1_rep1, var_method1_rep1, samples_not_used = junction_sample.sample_from_junctions(replica1, m, k, fit_func=fit_func1, **methods[method1_name])
        mean_method1_rep2, var_method1_rep2, samples_not_used = junction_sample.sample_from_junctions(replica2, m, k, fit_func=fit_func2, **methods[method1_name])
        score_method1 = calc_score_tmp(mean_method1_rep1, var_method1_rep1, mean_method1_rep2, var_method1_rep2)
        scores_cached[method1_name] = score_method1

    if method2_name in scores_cached:
        score_method2 = scores_cached[method2_name]
    else:
        mean_method2_rep1, var_method2_rep1, samples_not_used = junction_sample.sample_from_junctions(replica1, m, k, fit_func=fit_func1, **methods[method2_name])
        mean_method2_rep2, var_method2_rep2, samples_not_used = junction_sample.sample_from_junctions(replica2, m, k, fit_func=fit_func2, **methods[method2_name])
        score_method2 = calc_score_tmp(mean_method2_rep1, var_method2_rep1, mean_method2_rep2, var_method2_rep2)
        scores_cached[method2_name] = score_method2

    plot_method1Vsmethod2(score_method1, score_method2, method1_name,  method2_name, rep1_name, rep2_name, plotpath)

    # # Debugging odd results...
    # if 'Majiq' == method1_name:
    #     score_to_test = score_method1
    # elif 'Majiq' == method2_name:
    #     score_to_test = score_method2
    # for idx, lsv in enumerate(score_to_test):
    #     if lsv > 15:
    #         print lsv_junc1[1][idx][1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('par1', help='Path for parameters of replica1')
    parser.add_argument('par2', help='Path for parameters of replica2')  
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration')
    parser.add_argument('--m', default=100, type=int, help='Number of samples')     
    parser.add_argument('--junctype', default='rand10k', help='The type of junction to analyze. (Inc, Exc or rand10k for now)')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--output', default=None, help="Path to save the results to.")
    parser.add_argument('--minnonzero', default=1, type=int, help='Minimum number of positive positions to consider the junction')
    parser.add_argument('--maxnonzero', default=1000000, type=int, help='Maximum number of positive positions to consider the junction')
    parser.add_argument('--minreads', default=1, type=int, help='Minimum number of reads combining all positions in a junction to be considered')
    parser.add_argument('--dispersion', default=0.1, type=float, help='Dispersion factor (used in junctions sampling).')

    args = parser.parse_args()

    # Parse LSV files
    lsv_junc1, const1 = pipelines.load_data_lsv(args.par1)
    lsv_junc2, const2 = pipelines.load_data_lsv(args.par2)

    # Filter non-quantifiable
    lsv_junc1 = analysis.filter.lsv_quantifiable( lsv_junc1, args.minnonzero, args.minreads, None )
    lsv_junc2 = analysis.filter.lsv_quantifiable( lsv_junc2, args.minnonzero, args.minreads, None )

    replica1, replica2 = junction_sample.check_junctions_in_replicates(lsv_junc1, lsv_junc2,discard_empty_junctions=True)

    #Get the experiment names
    rep1_name = os.path.basename(args.par1).split('.')[-2]
    rep2_name = os.path.basename(args.par2).split('.')[-2]

    fit_func1 = polyfitnb.fit_nb(const1, "%s_nbfit" % args.output, args.plotpath, nbdisp=args.dispersion, logger=None, discardb=True)
    fit_func2 = polyfitnb.fit_nb(const2, "%s_nbfit" % args.output, args.plotpath, nbdisp=args.dispersion, logger=None, discardb=True)

    methods = {
        'Poisson':                  {'discardzeros': False, 'trimborder': False,   'nb': False},
        'Naive_Boots':              {'discardzeros': False, 'trimborder': False,   'nb': False},
        'Naive_Boots_trim_borders': {'discardzeros': True,  'trimborder': True,    'nb': False},
        'Naive_Boots_no_zeros':     {'discardzeros': True,  'trimborder': False,   'nb': False},
        'Neg_Binomial':             {'discardzeros': False, 'trimborder': False,   'nb': True},
        'Majiq':                    {'discardzeros': True,  'trimborder': True,    'nb': True},
        'Majiq_with_zeros':         {'discardzeros': False, 'trimborder': True,    'nb': True},
        'Majiq_no_stacks':          {'discardzeros': True,  'trimborder': True,    'nb': True},
    }

    scores_cached = {}
    compare_methods(args.m, args.k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Poisson',       'Majiq', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',   'Neg_Binomial', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',   'Naive_Boots_trim_borders', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',   'Majiq', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots_no_zeros',  'Majiq', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_with_zeros',      'Majiq', args.plotpath, scores_cached)

    #Majiq with stacks
    # marks_stacks = [0.0001, 0.0005, 0.001]  # Note that this array has to be sorted in ascending order
    # for mark_stacks in marks_stacks:
    #     analysis.filter.lsv_mark_stacks(lsv_junc1[0], fit_func1, mark_stacks, args.dispersion)
    #     analysis.filter.lsv_mark_stacks(lsv_junc2[0], fit_func2, mark_stacks, args.dispersion)
    #
    #     # Check that all LSVs have the same num. of junctions in both replicates
    #     replica1, replica2 = junction_sample.check_junctions_in_replicates(lsv_junc1, lsv_junc2)
    #
    #     compare_methods(args.m, args.k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_no_stack', 'Majiq', args.plotpath, scores_cached)
    #
    #     print len(replica1), len(replica2)


    # replica1_gcnorm = polyfit_ectf.norm_junctions(replica1, gc_factors=replica1_gc, gcnorm=True)
    # replica2_gcnorm = polyfit_ectf.norm_junctions(replica2, gc_factors=replica2_gc, gcnorm=True)

    #Sampling methods    REFERENCE (junctions, m, k, discardzeros=False, nb=False, trimborder=False, parameters=None)
    #Naive bootstrapping
    # boots_mean1, boots_var1, samples = junction_sample.sample_from_junctions(replica1, args.m, args.k, discardzeros=False, nb=False, trimborder=False)
    # boots_mean2, boots_var2, samples = junction_sample.sample_from_junctions(replica2, args.m, args.k, discardzeros=False, nb=False, trimborder=False)
    # #Naive bootstrapping with GC normalization
    # # boots_norm_mean1, boots_norm_var1, samples = junction_sample.sample_from_junctions(replica1_gcnorm, args.m, args.k, parameters=args.par1)
    # # boots_norm_mean2, boots_norm_var2, samples = junction_sample.sample_from_junctions(replica2_gcnorm, args.m, args.k, parameters=args.par2)
    #
    # #Naive bootstrapping with trimming of borders
    # # boots_trim_mean1, boots_trim_var1, samples = junction_sample.sample_from_junctions(replica1, args.m, args.k, discardzeros=False, nb=False, trimborder=False, fit_func=fit_func1)
    # # boots_trim_mean2, boots_trim_var2, samples = junction_sample.sample_from_junctions(replica2, args.m, args.k, discardzeros=False, nb=False, trimborder=False, fit_func=fit_func2)
    #
    # # #Naive bootstrapping discarding zeros
    # boots_zeros_mean1, boots_zeros_var1, samples = junction_sample.sample_from_junctions(replica1, args.m, args.k, discardzeros=True, nb=False, trimborder=False, fit_func=fit_func1)
    # boots_zeros_mean2, boots_zeros_var2, samples = junction_sample.sample_from_junctions(replica2, args.m, args.k, discardzeros=True, nb=False, trimborder=False, fit_func=fit_func2)
    # #Bootstrapping sampling from a NB distribution
    # boots_nb_mean1, boots_nb_var1, samples = junction_sample.sample_from_junctions(replica1, args.m, args.k, discardzeros=False, nb=True, trimborder=False, fit_func=fit_func1)
    # boots_nb_mean2, boots_nb_var2, samples = junction_sample.sample_from_junctions(replica2, args.m, args.k, discardzeros=False, nb=True, trimborder=False, fit_func=fit_func2)
    # #Majiq
    # majiq_mean1, majiq_var1, samples = junction_sample.sample_from_junctions(replica1, args.m, args.k, discardzeros=False, nb=True, trimborder=True, fit_func=fit_func1)
    # majiq_mean2, majiq_var2, samples = junction_sample.sample_from_junctions(replica2, args.m, args.k, discardzeros=False, nb=True, trimborder=True, fit_func=fit_func2)
    # #Majiq with zeros
    # majiq_mean_zeros1, majiq_var_zeros1, samples = junction_sample.sample_from_junctions(replica1, args.m, args.k, discardzeros=True, nb=True, trimborder=True, fit_func=fit_func1)
    # majiq_mean_zeros2, majiq_var_zeros2, samples = junction_sample.sample_from_junctions(replica2, args.m, args.k, discardzeros=True, nb=True, trimborder=True, fit_func=fit_func2)
    #
    # #calculate scores for different methods
    # #scores_empirical    = calc_score(empirical_mean1, empirical_var1, empirical_mean2, empirical_var2)
    # scores_poisson      = calc_score(boots_mean1, boots_mean1, boots_mean2, boots_mean2)
    # scores_boots        = calc_score(boots_mean1, boots_var1, boots_mean2, boots_var2)
    # # scores_boots_norm   = calc_score(boots_norm_mean1, boots_norm_var1, boots_norm_mean2, boots_norm_var2)
    # # scores_boots_trim   = calc_score(boots_trim_mean1, boots_trim_var1, boots_trim_mean2, boots_trim_var2)
    # scores_boots_zero   = calc_score(boots_zeros_mean1, boots_zeros_var1, boots_zeros_mean2, boots_zeros_var2)
    # scores_boots_nb     = calc_score(boots_nb_mean1, boots_nb_var1, boots_nb_mean2, boots_nb_var2)
    # scores_majiq        = calc_score(majiq_mean1, majiq_var1, majiq_mean2, majiq_var2)
    # scores_majiq_zeros  = calc_score(majiq_mean_zeros1, majiq_var_zeros1, majiq_mean_zeros2, majiq_var_zeros2)


    # #Majiq with stacks
    # # marks_stacks = [0.0001, 0.0005, 0.001]
    # # for mark_stacks in marks_stacks:
    #
    # #
    # #     analysis.filter.lsv_mark_stacks(lsv_junc1[0], fit_func1, mark_stacks, args.dispersion)
    # #     analysis.filter.lsv_mark_stacks(lsv_junc2[0], fit_func2, mark_stacks, args.dispersion)
    # #     lsv_junc1 = analysis.filter.lsv_quantifiable( lsv_junc1, args.minnonzero, args.minreads, None )
    # #     lsv_junc2 = analysis.filter.lsv_quantifiable( lsv_junc2, args.minnonzero, args.minreads, None )
    # #
    # #     # Check that all LSVs have the same num. of junctions in both replicates
    # #     replica1, replica2 = junction_sample.check_junctions_in_replicates(lsv_junc1, lsv_junc2)
    # #
    # #     print len(replica1), len(replica2)
    # #
    # #     majiq_mean_stacks1, majiq_var_stacks1, samples = junction_sample.sample_from_junctions(replica1, args.m, args.k, discardzeros=True, nb=True, trimborder=True, fit_func=fit_func1)
    # #     majiq_mean_stacks2, majiq_var_stacks2, samples = junction_sample.sample_from_junctions(replica2, args.m, args.k, discardzeros=True, nb=True, trimborder=True, fit_func=fit_func2)
    # #
    # #     #calculate scores for different methods
    # #     scores_majiq_stacks = calc_score(majiq_mean_stacks1, majiq_var_stacks1, majiq_mean_stacks2, majiq_var_stacks2)
    # #
    # #     # Majiq vs Majiq Stacks
    # #     plot_varvsvar(scores_majiq_zeros, scores_majiq_stacks, "MAJIQ",  "MAJIQ (No Stacks)", name1, name2, args.plotpath, "MajiqVsMajiqStacks_%f" % mark_stacks)
    #
    #

    # print "#---- Summary -----#"
    # #Strawman
    # plot_varvsvar(scores_poisson, scores_majiq_zeros, "Poisson",  "MAJIQ", name1, name2, args.plotpath)
    # #Testing GC normalization
    # # plot_varvsvar(scores_boots, scores_boots_norm, "Naive Bootstrapping", "Naive Bootstrapping GC norm", name1, name2, args.plotpath, "GCnorm")
    # #testing Parametric NB sampling
    # plot_varvsvar(scores_boots, scores_boots_nb, "Naive Bootstrapping", "Negative Binomial",  name1, name2, args.plotpath)
    # #testing trimmming
    # # plot_varvsvar(scores_boots, scores_boots_trim, "Naive bootstrapping", "Naive bootstrapping (trimming borders)", name1, name2, args.plotpath, "Trimming")
    # #Naive Vs Majiq
    # plot_varvsvar(scores_boots, scores_majiq_zeros, "Naive Bootstrapping", "MAJIQ", name1, name2, args.plotpath)
    # #Naive (nozeros) Vs Majiq
    # plot_varvsvar(scores_boots_zero, scores_majiq_zeros, "Naive Bootstrapping (No Zeros)", "MAJIQ", name1, name2, args.plotpath)
    # # Naive Poisson vs Naive Boostraping
    # plot_varvsvar(scores_poisson, scores_boots, "Poisson",  "Naive Bootstraping", name1, name2, args.plotpath)
    # # Majiq vs Majiq (No Zeros)
    # plot_varvsvar(scores_majiq, scores_majiq_zeros, "MAJIQ",  "MAJIQ (No Zeros)", name1, name2, args.plotpath)
    #
    # # Debugging odd results...
    # for idx, lsv in enumerate(scores_majiq_zeros):
    #     if lsv > 15:
    #         print lsv_junc1[1][idx][1]



if __name__ == '__main__':
    main()
