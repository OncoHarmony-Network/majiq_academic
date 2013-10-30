import sys
import os
import argparse

from scipy.io import loadmat
from pylab import *
from numpy.ma import masked_less
from matplotlib import rcParams

from junction_sample import sample_from_junctions, norm_junctions, calc_nonzeromeanvar, discardlow, discardhigh, DEBUG, TESTBREAK, EPSILON


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


def calc_score(mean_sample1, var_sample1, mean_sample2, var_sample2):
    return abs(var_sample1-var_sample2)/((mean_sample1+mean_sample2)*0.5)



def plot_varvsvar(score1, score2, score1_name, score2_name, replica1_name, replica2_name, plotpath=None, plotname=None):
    total_junctions = float(len(score1))
    better_in_method1 = sum(score1 < score2)
    equal_in_both = sum(score1 == score2)
    better_in_method2 = sum(score1 > score2)
    """
    print score1_name
    print score1
    print
    print score2_name
    print score2
    print
    print better_in_method1
    print better_in_method2
    """
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('matpath', help='Path with matfile with replica1 and replica2')    
    parser.add_argument('par1', help='Path for parameters of replica1')
    parser.add_argument('par2', help='Path for parameters of replica2')  
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration')
    parser.add_argument('--m', default=100, type=int, help='Number of samples')     
    parser.add_argument('--junctype', default='rand10k', help='The type of junction to analyze. (Inc, Exc or rand10k for now)')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--minnonzero', default=0, type=int, help='Minimum number of positive positions to consider the junction')   
    parser.add_argument('--maxnonzero', default=0, type=int, help='Maximum number of positive positions to consider the junction') 
    args = parser.parse_args()

    my_mat = loadmat(args.matpath)
    name1, name2 = os.path.basename(args.matpath).split('.')[0].split("_") #Get the experiment names from the mat file
    main_title = "%s VS %s\n"%(name1, name2) #main title of the plot generation

    replica1 = my_mat[name1][args.junctype][0, 0][0, 0]['cov']
    replica2 = my_mat[name2][args.junctype][0, 0][0, 0]['cov']
    replica1_gc = my_mat[name1][args.junctype][0, 0][0, 0]['gc_val']
    replica2_gc = my_mat[name2][args.junctype][0, 0][0, 0]['gc_val']

    if args.maxnonzero:
        print "Before Maxnonzero %s. replica1: %s replica2: %s"%(args.maxnonzero, replica1.shape, replica2.shape)
        replica1, replica1_gc, replica2, replica2_gc = discardhigh(replica1, replica1_gc, replica2, replica2_gc, args.maxnonzero)
        print "After Maxnonzero %s. replica1: %s replica2: %s"%(args.maxnonzero, replica1.shape, replica2.shape)      

    if args.minnonzero:
        print "Before Minnonzero %s. replica1: %s replica2: %s"%(args.minnonzero, replica1.shape, replica2.shape)
        replica1, replica1_gc, replica2, replica2_gc = discardlow(replica1, replica1_gc, replica2, replica2_gc, args.minnonzero)
        print "After Minnonzero %s. replica1: %s replica2: %s"%(args.minnonzero, replica1.shape, replica2.shape)      

    #every method discards the -1, no one likes them :'(
    replica1 = masked_less(replica1, 0) 
    replica2 = masked_less(replica2, 0) 
    replica1_gcnorm = norm_junctions(replica1, gc_factors=replica1_gc, gcnorm=True)
    replica2_gcnorm = norm_junctions(replica2, gc_factors=replica2_gc, gcnorm=True)
    """
    empirical_mean1 = replica1.mean(axis=1)
    empirical_var1 = replica1.var(axis=1)
    empirical_mean2 = replica2.mean(axis=1)
    empirical_var2 = replica2.var(axis=1)
    empnozero_mean1, empnozero_var1 = calc_nonzeromeanvar(replica1)
    empnozero_mean2, empnozero_var2 = calc_nonzeromeanvar(replica2)
    if DEBUG:
        empirical_mean1 = empirical_mean1[:TESTBREAK]
        empirical_mean2 = empirical_mean2[:TESTBREAK]
        empnozero_mean1 = empnozero_mean1[:TESTBREAK]
        empnozero_mean2 = empnozero_mean2[:TESTBREAK]         
        empirical_var1 = empirical_var1[:TESTBREAK]
        empirical_var2 = empirical_var2[:TESTBREAK]
        empnozero_var1 = empnozero_var1[:TESTBREAK]
        empnozero_var2 = empnozero_var2[:TESTBREAK] 
    """
    replica1_gcnorm = norm_junctions(replica1, gc_factors=replica1_gc, gcnorm=True)
    replica2_gcnorm = norm_junctions(replica2, gc_factors=replica2_gc, gcnorm=True)

    #Sampling methods    REFERENCE (junctions, m, k, discardzeros=False, nb=False, trimborder=False, parameters=None)
    #Naive bootstrapping
    boots_mean1, boots_var1 = sample_from_junctions(replica1, args.m, args.k, parameters=args.par1)
    boots_mean2, boots_var2 = sample_from_junctions(replica2, args.m, args.k, parameters=args.par2)
    #Naive bootstrapping with GC normalization
    boots_norm_mean1, boots_norm_var1 = sample_from_junctions(replica1_gcnorm, args.m, args.k, parameters=args.par1)
    boots_norm_mean2, boots_norm_var2 = sample_from_junctions(replica2_gcnorm, args.m, args.k, parameters=args.par2)
    #Naive bootstrapping with trimming of borders
    boots_trim_mean1, boots_trim_var1 = sample_from_junctions(replica1_gcnorm, args.m, args.k, trimborder=True, parameters=args.par1)
    boots_trim_mean2, boots_trim_var2 = sample_from_junctions(replica2_gcnorm, args.m, args.k, trimborder=True, parameters=args.par2)
    #Naive bootstrapping discarding zeros
    boots_zeros_mean1, boots_zeros_var1 = sample_from_junctions(replica1_gcnorm, args.m, args.k, discardzeros=True, parameters=args.par1)
    boots_zeros_mean2, boots_zeros_var2 = sample_from_junctions(replica2_gcnorm, args.m, args.k, discardzeros=True, parameters=args.par2)
    #Bootstrapping sampling from a NB distribution
    boots_nb_mean1, boots_nb_var1  = sample_from_junctions(replica1_gcnorm, args.m, args.k, nb=True, parameters=args.par1)
    boots_nb_mean2, boots_nb_var2  = sample_from_junctions(replica2_gcnorm, args.m, args.k, nb=True, parameters=args.par2)
    #Majiq with zeros
    majiq_mean1, majiq_var1 = sample_from_junctions(replica1_gcnorm, args.m, args.k,  discardzeros=True, nb=True, trimborder=True, parameters=args.par1)
    majiq_mean2, majiq_var2  = sample_from_junctions(replica2_gcnorm, args.m, args.k, discardzeros=True, nb=True, trimborder=True, parameters=args.par2)

    #calculate scores for different methods
    #scores_empirical    = calc_score(empirical_mean1, empirical_var1, empirical_mean2, empirical_var2)    
    scores_poisson      = calc_score(boots_zeros_mean1, boots_zeros_var1, boots_zeros_mean2, boots_zeros_mean2)
    scores_boots        = calc_score(boots_mean1, boots_var1, boots_mean2, boots_var2)
    scores_boots_norm   = calc_score(boots_norm_mean1, boots_norm_var1, boots_norm_mean2, boots_norm_var2)
    scores_boots_trim   = calc_score(boots_trim_mean1, boots_trim_var1, boots_trim_mean2, boots_trim_var2)
    scores_boots_zero   = calc_score(boots_zeros_mean1, boots_zeros_var1, boots_zeros_mean2, boots_zeros_var2)
    scores_boots_nb     = calc_score(boots_nb_mean1, boots_nb_var1, boots_nb_mean2, boots_nb_var2)
    scores_majiq        = calc_score(majiq_mean1, majiq_var1, majiq_mean2, majiq_var2)

    #Strawman
    plot_varvsvar(scores_poisson, scores_majiq, "Poisson",  "MAJIQ", name1, name2, args.plotpath, "PoissonvsMAJIQ")
    #Testing GC normalization
    plot_varvsvar(scores_boots, scores_boots_norm, "Naive Bootstrapping", "Naive Bootstrapping GC norm", name1, name2, args.plotpath, "GCnorm")
    #testing Parametric NB sampling
    plot_varvsvar(scores_boots_norm, scores_boots_nb, "Naive bootstrapping", "Parametric Bootstrapping",  name1, name2, args.plotpath, "NBVSparametric")
    #testing trimmming 
    plot_varvsvar(scores_boots_norm, scores_boots_trim, "Naive bootstrapping", "Naive bootstrapping (trimming borders)", name1, name2, args.plotpath, "Trimming")
    #FINAL test
    plot_varvsvar(scores_boots_zero, scores_majiq, "Naive Bootstrapping", "MAJIQ", name1, name2, args.plotpath, "naiveNoZeroVSMAJIQ")




if __name__ == '__main__':
    main()
