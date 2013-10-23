import sys
import os
import argparse

from scipy.io import loadmat
from pylab import *
from numpy.ma import masked_less
from matplotlib import rcParams

from junction_sample import sample_from_junctions, norm_junctions, calc_nonzeromeanvar, DEBUG, TESTBREAK, EPSILON


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


def plot_varvsvar(score1, score2, my_title="", score1_name="", score2_name="", plotpath=None):
    title(my_title)
    max_value = max(max(score1), max(score2))
    xlabel(score2_name)
    ylabel(score1_name)
    xlim(0, max_value)
    ylim(0, max_value)
    plot([0, max_value], [0, max_value])
    plot(score1, score2, '.')
    _save_or_show(plotpath)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('matpath', help='Path with matfile with replica1 and replica2')    
    parser.add_argument('par1', help='Path for parameters of replica1')
    parser.add_argument('par2', help='Path for parameters of replica2')  
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration')
    parser.add_argument('--m', default=100, type=int, help='Number of samples') 
    parser.add_argument('--junctype', default='rand10k', help='The type of junction to analyze. (Inc, Exc or rand10k for now)')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    args = parser.parse_args()

    my_mat = loadmat(args.matpath)
    name1, name2 = os.path.basename(args.matpath).split('.')[0].split("_") #Get the experiment names from the mat file
    main_title = "%s VS %s\n"%(name1, name2) #main title of the plot generation

    replica1 = my_mat[name1][args.junctype][0, 0][0, 0]['cov']
    replica2 = my_mat[name2][args.junctype][0, 0][0, 0]['cov']
    replica1_gc = my_mat[name1][args.junctype][0, 0][0, 0]['gc_val']
    replica2_gc = my_mat[name2][args.junctype][0, 0][0, 0]['gc_val']

    #every method discards the -1, no one likes them :'(
    #replica1 = replica1[numpy > -EPSILON]  
    #replica2 = replica2[replica2 > -EPSILON] 
    replica1 = masked_less(replica1, 0) 
    replica2 = masked_less(replica2, 0) 

    print replica1

    empirical_mean1 = replica1.mean(axis=1)
    empirical_var1 = replica1.var(axis=1)
    empirical_mean2 = replica2.mean(axis=1)
    empirical_var2 = replica2.var(axis=1)
    empnozero_mean1, empnozero_var1 = calc_nonzeromeanvar(replica1)
    empnozero_mean2, empnozero_var2 = calc_nonzeromeanvar(replica2)

    if DEBUG:
        empirical_var1 = empirical_var1[:TESTBREAK]
        empirical_var2 = empirical_var2[:TESTBREAK]
        empnozero_var1 = empnozero_var1[:TESTBREAK]
        empnozero_var2 = empnozero_var2[:TESTBREAK] 

    replica1_gcnorm = norm_junctions(replica1, gc_factors=replica1_gc, gcnorm=True)
    replica2_gcnorm = norm_junctions(replica2, gc_factors=replica2_gc, gcnorm=True)

    #Bootstrapping method
    bootstrapping_zeros_mean1, bootstrapping_zeros_var1 = sample_from_junctions(replica1, args.m, args.k, discardzeros=False, nb=False, parameters=args.par1)
    bootstrapping_zeros_mean2, bootstrapping_zeros_var2 = sample_from_junctions(replica2, args.m, args.k, discardzeros=False, nb=False, parameters=args.par2)
    #Majiq without normalizing
    majiq_nonorm_mean1, majiq_nonorm_var1 = sample_from_junctions(replica1, args.m, args.k, discardzeros=False, nb=True, parameters=args.par1)
    majiq_nonorm_mean2, majiq_nonorm_var2  = sample_from_junctions(replica2, args.m, args.k, discardzeros=False, nb=True, parameters=args.par2)
    #Majiq with zeros
    majiq_zero_mean1, majiq_zero_var1 = sample_from_junctions(replica1_gcnorm, args.m, args.k, discardzeros=False, nb=True, parameters=args.par1)
    majiq_zero_mean2, majiq_zero_var2  = sample_from_junctions(replica1_gcnorm, args.m, args.k, discardzeros=False, nb=True, parameters=args.par2)

    minnonzero = 5
    num_before1 = replica1_gcnorm.shape[0]
    num_before2 = replica2_gcnorm.shape[0]
    numzeros1 = _numzeros(replica1_gcnorm) 
    numzeros2 = _numzeros(replica2_gcnorm)   
    pass_threshold1 = (numzeros1 < replica1_gcnorm.shape[1]-minnonzero)
    pass_threshold2 = (numzeros2 < replica2_gcnorm.shape[1]-minnonzero)
    pass_threshold = all([pass_threshold1, pass_threshold2], axis=0)
    print sum(pass_threshold)
    print sum(pass_threshold1)
    print sum(pass_threshold2)
    #junctions = junctions[pass_threshold]
    replica1 = replica1[pass_threshold]
    replica2 = replica2[pass_threshold]  
    replica1_gcnorm = replica1_gcnorm[pass_threshold]
    replica2_gcnorm = replica2_gcnorm[pass_threshold] 

    #Bootstrapping method without zeroes
    bootstrapping_nozeros_mean1, bootstrapping_nozeros_var1 = sample_from_junctions(replica1, args.m, args.k, discardzeros=True, nb=False, parameters=args.par1)
    bootstrapping_nozeros_mean2, bootstrapping_nozeros_var2 = sample_from_junctions(replica2, args.m, args.k, discardzeros=True, nb=False, parameters=args.par2)
    #Majiq without zeros
    majiq_nozero_mean1, majiq_nozero_var1 = sample_from_junctions(replica1_gcnorm, args.m, args.k, discardzeros=True, nb=True, parameters=args.par1)
    majiq_nozero_mean2, majiq_nozero_var2 = sample_from_junctions(replica2_gcnorm, args.m, args.k, discardzeros=True, nb=True, parameters=args.par2)
    #calculate scores for different methods
    scores_poisson       = abs(bootstrapping_zeros_mean1-empirical_var2)
    scores_zero_boots = abs(bootstrapping_zeros_var1-empirical_var2)
    scores_nozero_boots = abs(bootstrapping_nozeros_var1-empnozero_var2)
    scores_notnorm_majiq = abs(majiq_nonorm_var1-empirical_var2)
    scores_zero_majiq    = abs(majiq_zero_var1-empirical_var2)
    scores_nonzero_majiq = abs(majiq_nozero_var1-empnozero_var2)

    print "Better junctions in MAJIQ than Poisson: %s (%.2f%%)"%(sum(scores_poisson > scores_zero_majiq), (sum(scores_poisson > scores_zero_majiq)/float(len(scores_zero_majiq)))*100)
    print "Better junctions in MAJIQ than Bootstrapping (with zeros): %s (%.2f%%)"%(sum(scores_zero_boots > scores_zero_majiq), (sum(scores_zero_boots > scores_zero_majiq)/float(len(scores_zero_majiq)))*100)
    print "Better junctions in MAJIQ than Bootstrapping (without zeros): %s (%.2f%%)"%(sum(scores_nozero_boots > scores_nonzero_majiq), (sum(scores_nozero_boots > scores_nonzero_majiq)/float(len(scores_nonzero_majiq)))*100)
    print "Better junctions in MAJIQ without 0s than MAJIQ with 0s: %s (%.2f%%)"%(sum(scores_zero_majiq > scores_nonzero_majiq), (sum(scores_zero_majiq > scores_nonzero_majiq)/float(len(scores_zero_majiq)))*100)

    #MAJIQ vs Poisson
    plot_varvsvar(scores_poisson, scores_zero_majiq, "Poisson sampling variance(=mean) vs MAJIQ variance\n(with 0s)", 
                                                     "Poisson %s - Empirical %s"%(name1, name2), 
                                                     "MAJIQ %s - Empirical %s"%(name1, name2))
    plot_varvsvar(scores_zero_boots, scores_zero_majiq, "Bootstrapping sampling variance vs MAJIQ variance\n(with 0s)", 
                                                     "Bootstrapping %s - Empirical %s"%(name1, name2), 
                                                     "MAJIQ %s - Empirical %s"%(name1, name2))

    plot_varvsvar(scores_zero_majiq, scores_nonzero_majiq, "MAJIQ (with 0s) variance vs MAJIQ (without 0s)\n(with 0s)", 
                                                           "MAJIQ (with 0) %s - Empirical (with 0) %s"%(name1, name2), 
                                                           "MAJIQ (without 0) %s - Empirical (without 0) %s"%(name1, name2))


    plot_varvsvar(scores_zero_majiq, scores_nonzero_majiq, "MAJIQ (with 0s) variance vs MAJIQ (without 0s)\n(with 0s)", 
                                                           "MAJIQ (with 0) %s - Empirical (with 0) %s"%(name1, name2), 
                                                           "MAJIQ (without 0) %s - Empirical (without 0) %s"%(name1, name2))


if __name__ == '__main__':
    main()
