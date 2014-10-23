from collections import defaultdict
from pdb import set_trace
from matplotlib import use
from numpy.ma import masked_less
from scripts.utils import _save_or_show

use('Agg', warn=False)

from scipy.stats import pearsonr
import argparse
from pylab import *
import analysis.filter
import junction_sample
from scipy.stats import nbinom
import os.path
from analysis import filter
from analysis import polyfitnb


def _numzeros(junctions):
    "Obtain the number of zeros for every junction"
    return (junctions == 0).sum(axis=1)


def calc_score_tmp(mean_sample1, var_sample1, mean_sample2, var_sample2):
    return {'mean': abs(mean_sample1-mean_sample2)/((mean_sample1+mean_sample2)*0.5),
            'variance': abs(var_sample1-var_sample2)/((var_sample1+var_sample2)*0.5)}

def calc_score(mean_sample1, var_sample1, mean_sample2, var_sample2):
    return calc_score_tmp(mean_sample1, var_sample1, mean_sample2, var_sample2)
    # return abs(var_sample1-var_sample2)/((mean_sample1+mean_sample2)*0.5)

def plot_method1Vsmethod2(scores1, scores2, score1_name, score2_name, replica1_name, replica2_name, plotpath=None):
    """Compute 2 plots: one for the variance, one for the mean"""
    plotname="%sVs%s" % (score1_name, score2_name)

    total_junctions = float(len(scores1['mean']))  # The same in all mean/variance scores1/scores2

    stats_names = ['mean', 'variance']

    f, axarr = subplots(2, 1, sharex=True)
    axarr[0].set_ylabel(score2_name)
    axarr[-1].set_xlabel(score1_name)

    for name_i, stat_name in enumerate(stats_names):

        axarr[name_i].set_title(stat_name, fontsize=10)

        score1, score2 = scores1[stat_name], scores2[stat_name]
        better_in_method1 = sum(score1 < score2)
        equal_in_both = sum(score1 == score2)
        better_in_method2 = sum(score1 > score2)

        max_value = max(max(score1), max(score2))

        xlim(0, max_value)
        ylim(0, max_value)
        axarr[name_i].plot([0, max_value], [0, max_value])
        axarr[name_i].plot(score1, score2, '.')

        pear, pvalue = pearsonr(score1, score2)
        r_squared = pear**2

        # axarr[name_i].text(abs(max_value)*0.1, max_value-abs(max_value)*0.2, r'$R^2$: %.2f (p-value: %.2E)'%(r_squared, pvalue), fontsize=12, bbox={'facecolor':'yellow', 'alpha':0.3, 'pad':10})
        axarr[name_i].text(abs(max_value)*0.1, max_value-abs(max_value)*0.2, '%.2f%%' % ((better_in_method1/total_junctions)*100), fontsize=16)
        axarr[name_i].text(abs(max_value)*0.8, max_value-abs(max_value)*0.8, '%.2f%%' % ((better_in_method2/total_junctions)*100), fontsize=16)

        print "[For %8s]:: Better in %s: %s (%.2f%%) Equal: %s Better in %s: %s (%.2f%%)"%(stat_name, score1_name, better_in_method1, (better_in_method1/total_junctions)*100, equal_in_both, score2_name, better_in_method2, (better_in_method2/total_junctions)*100)

    suptitle("%s vs %s\n%s and %s"%(score1_name, score2_name, replica1_name, replica2_name))

    _save_or_show(plotpath, plotname)


def plot_method1Vsmethod2_by_coverage(scores1_list, scores2_list, score1_name, score2_name, replica1_name, replica2_name, plotpath=None):
    """Compute 6 plots, one for each coverage subset, variance and mean"""

    # TODO: THIS IS NOT WORKING PROPERLY YET. Review why scores1 has zeros
    plotname="%sVs%s-Coverage" % (score1_name, score2_name)

    stats_names = ['mean', 'variance']

    f, axarr = subplots(len(stats_names), len(scores1_list), sharex=True, sharey=True)
    axarr[0][0].set_ylabel(score2_name)
    axarr[0][-1].set_xlabel(score1_name)


    for score_idx in range(len(scores1_list)):
        total_junctions = float(len(scores1_list[score_idx]['mean']))  # The same in all mean/variance scores1/scores2

        for name_i, stat_name in enumerate(stats_names):
            axarr[name_i][0].set_title(stat_name, fontsize=10)
            score1, score2 = scores1_list[score_idx][stat_name], scores2_list[score_idx][stat_name]

            better_in_method1 = sum(score1 < score2)
            equal_in_both = sum(score1 == score2)
            better_in_method2 = sum(score1 > score2)

            max_value = max(max(score1), max(score2))

            xlim(0, max_value)
            ylim(0, max_value)
            axarr[name_i][score_idx].plot([0, max_value], [0, max_value])
            axarr[name_i][score_idx].plot(score1, score2, '.')

            pear, pvalue = pearsonr(score1, score2)
            r_squared = pear**2

            # axarr[name_i][score_idx].text(abs(max_value)*0.1, max_value-abs(max_value)*0.2, r'$R^2$: %.2f (p-value: %.2E)'%(r_squared, pvalue), fontsize=12, bbox={'facecolor':'yellow', 'alpha':0.3, 'pad':10})
            axarr[name_i][score_idx].text(abs(max_value)*0.2, max_value-abs(max_value)*0.2, "%.2f%%" % ((better_in_method1/total_junctions)*100), fontsize=12)
            axarr[name_i][score_idx].text(abs(max_value)*0.7, max_value-abs(max_value)*0.8, "%.2f%%" % ((better_in_method2/total_junctions)*100), fontsize=12)

            print "[For %8s]:: Better in %s: %s (%.2f%%) Equal: %s Better in %s: %s (%.2f%%)"%(stat_name, score1_name, better_in_method1, (better_in_method1/total_junctions)*100, equal_in_both, score2_name, better_in_method2, (better_in_method2/total_junctions)*100)

    suptitle("%s vs %s\n%s and %s"%(score1_name, score2_name, replica1_name, replica2_name))

    _save_or_show(plotpath, plotname)


def plot_sigmaVsMu_met1_met2(mean_var_by_replicate, scores1, scores2, score1_name, score2_name, plotpath=None):
    """Compute 4 plots: for the variance and the mean, for method 1 and method 2."""
    stats_names = ('mean', 'variance')
    plotname="%sVs%s-%sVs%s" % (stats_names[0], stats_names[1], score1_name, score2_name)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)    # The big subplot
    axs = [
        fig.add_subplot(2, 2, 1),
        fig.add_subplot(2, 2, 2),
        fig.add_subplot(2, 2, 3),
        fig.add_subplot(2, 2, 4)
    ]

    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    # method1 == 'poisson'
    method1_var_wins = scores1['variance'] < scores2['variance']
    method1_var_lose = scores1['variance'] > scores2['variance']

    max_value = 50

    for idx in range(4):
        if idx % 2:
            score_name = score2_name
            axs[idx].set_ylabel('variance', fontsize=8)
            method_index = '2'

        else:
            score_name = score1_name
            if idx:
                tail = "win"
            else:
                tail = "lose"
            axs[idx].set_ylabel("%s_%s\nvariance" % (score1_name, tail), fontsize=8)
            method_index = '1'

        if idx < 2:
            axs[idx].set_title(score_name, fontsize=12)
            method = method1_var_lose
        else:
            method = method1_var_wins

        mean_index_str = 'mean' + method_index
        var_index_str = 'var' + method_index

        axs[idx].plot(mean_var_by_replicate['rep1'][mean_index_str]*method, mean_var_by_replicate['rep1'][var_index_str]*method, '.')
        axs[idx].plot(mean_var_by_replicate['rep2'][mean_index_str]*method, mean_var_by_replicate['rep2'][var_index_str]*method, '.')
        axs[idx].plot([0, max_value], [0, max_value])

        axs[idx].set_xlabel('mean', fontsize=8)

        axs[idx].set_xlim(0, 50)
        if idx == 3:
            axs[idx].set_ylim(0, 500)

        axs[idx].grid(True)
    _save_or_show(plotpath, plotname)


def compare_methods(m, k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, method1_name, method2_name, plotpath, scores_cached, coverage=False, rep1_gc=None, rep2_gc=None):

    if coverage:
        rep1_pool, rep2_pool = junction_sample.split_junction_pool(replica1, replica2)
        score_method1_list = []
        score_method2_list = []

        for rep_index in range(0, len(rep1_pool)):
            mean_method1_rep1, var_method1_rep1, samples_not_used = junction_sample.sample_from_junctions(rep1_pool[rep_index], m, k, fit_func=fit_func1, poisson=method1_name=='Poisson', **methods[method1_name])
            mean_method1_rep2, var_method1_rep2, samples_not_used = junction_sample.sample_from_junctions(rep2_pool[rep_index], m, k, fit_func=fit_func2, poisson=method1_name=='Poisson', **methods[method1_name])
            score_method1_list.append(calc_score_tmp(mean_method1_rep1, var_method1_rep1, mean_method1_rep2, var_method1_rep2))
            mean_method2_rep1, var_method2_rep1, samples_not_used = junction_sample.sample_from_junctions(rep1_pool[rep_index], m, k, fit_func=fit_func1, poisson=method2_name=='Poisson', **methods[method2_name])
            mean_method2_rep2, var_method2_rep2, samples_not_used = junction_sample.sample_from_junctions(rep2_pool[rep_index], m, k, fit_func=fit_func2, poisson=method2_name=='Poisson', **methods[method2_name])
            score_method2_list.append(calc_score_tmp(mean_method2_rep1, var_method2_rep1, mean_method2_rep2, var_method2_rep2))

        plot_method1Vsmethod2_by_coverage(score_method1_list, score_method2_list, method1_name,  method2_name, rep1_name, rep2_name, plotpath)

    elif rep1_gc:  # GC normalitation

        if 'gc' in method1_name:
            replica_data1 = rep1_gc
            replica_data2 = rep2_gc
        else:
            replica_data1 = replica1
            replica_data2 = replica2

        mean_method1_rep1, var_method1_rep1, samples_not_used = junction_sample.sample_from_junctions(replica_data1, m, k, fit_func=fit_func1, poisson=method1_name=='Poisson', **methods[method1_name])
        mean_method1_rep2, var_method1_rep2, samples_not_used = junction_sample.sample_from_junctions(replica_data2, m, k, fit_func=fit_func2, poisson=method1_name=='Poisson', **methods[method1_name])
        score_method1 = calc_score_tmp(mean_method1_rep1, var_method1_rep1, mean_method1_rep2, var_method1_rep2)

        if 'gc' in method2_name:
            replica_data1 = rep1_gc
            replica_data2 = rep2_gc
        else:
            replica_data1 = replica1
            replica_data2 = replica2

        mean_method2_rep1, var_method2_rep1, samples_not_used = junction_sample.sample_from_junctions(replica_data1, m, k, fit_func=fit_func1, poisson=method2_name=='Poisson', **methods[method2_name])
        mean_method2_rep2, var_method2_rep2, samples_not_used = junction_sample.sample_from_junctions(replica_data2, m, k, fit_func=fit_func2, poisson=method2_name=='Poisson', **methods[method2_name])
        score_method2 = calc_score_tmp(mean_method2_rep1, var_method2_rep1, mean_method2_rep2, var_method2_rep2)

        plot_method1Vsmethod2(score_method1, score_method2, method1_name,  method2_name, rep1_name, rep2_name, plotpath)

    else:
        mean_method1_rep1, var_method1_rep1, samples_not_used = junction_sample.sample_from_junctions(replica1, m, k, fit_func=fit_func1, poisson=method1_name=='Poisson', **methods[method1_name])
        mean_method1_rep2, var_method1_rep2, samples_not_used = junction_sample.sample_from_junctions(replica2, m, k, fit_func=fit_func2, poisson=method1_name=='Poisson', **methods[method1_name])
        score_method1 = calc_score_tmp(mean_method1_rep1, var_method1_rep1, mean_method1_rep2, var_method1_rep2)

        mean_method2_rep1, var_method2_rep1, samples_not_used = junction_sample.sample_from_junctions(replica1, m, k, fit_func=fit_func1, poisson=method2_name=='Poisson', **methods[method2_name])
        mean_method2_rep2, var_method2_rep2, samples_not_used = junction_sample.sample_from_junctions(replica2, m, k, fit_func=fit_func2, poisson=method2_name=='Poisson', **methods[method2_name])
        score_method2 = calc_score_tmp(mean_method2_rep1, var_method2_rep1, mean_method2_rep2, var_method2_rep2)

        plot_method1Vsmethod2(score_method1, score_method2, method1_name,  method2_name, rep1_name, rep2_name, plotpath)

    # mean_var_by_replicate = {
    #         'rep1': {
    #             'name'  : rep1_name,
    #             'mean1' : mean_method1_rep1,
    #             'var1'  : var_method1_rep1,
    #             'mean2' : mean_method2_rep1,
    #             'var2'  : var_method2_rep1
    #         },
    #         'rep2': {
    #             'name'  : rep2_name,
    #             'mean1' : mean_method1_rep2,
    #             'var1'  : var_method1_rep2,
    #             'mean2' : mean_method2_rep2,
    #             'var2'  : var_method2_rep2
    #         }
    #     }
    # plot_sigmaVsMu_met1_met2(mean_var_by_replicate, score_method1, score_method2, method1_name,  method2_name, plotpath)


def lsv_mark_stacks(lsv_list, fitfunc_r, pvalue_limit, dispersion, logger=None):

    minstack = sys.maxint
     #the minimum value marked as stack
    numstacks = 0
    filtered = [False]*len(lsv_list[0])
    stack_junc_idxs = defaultdict(int)
    for lidx, junctions in enumerate(lsv_list[0]):

        for i, junction in enumerate(junctions):
            if np.count_nonzero(junction) == 0:
                continue
            for j, value in enumerate(junction):
                if value > 0:
                    #TODO Use masker, and marking stacks will probably be faster.
                    copy_junc = list(junction)
                    copy_junc.pop(j)
                    copy_junc = np.array(copy_junc)
                    copy_junc = copy_junc[copy_junc > 0]

                    #FINISH TODO
                    mean_rest = np.mean(copy_junc)
                    pval = polyfitnb.get_negbinom_pval(fitfunc_r, mean_rest, value)
                    if pval < pvalue_limit:
                        lsv_list[0][lidx][i, j] = -2
                        minstack = min(minstack, value)
                        numstacks += 1
                        filtered[lidx] = True
                        stack_junc_idxs[lsv_list[1][lidx][1]] = i

        masked_less(lsv_list[0][lidx], 0)
    filtered = np.array(filtered)
    print "LSVs marked with stacks: %d" % np.count_nonzero(filtered)
    return [lsv_list[0][filtered], lsv_list[1][filtered]], stack_junc_idxs


def mark_stacks(junctions, fitted_1_r, pvalue_limit, dispersion, logger=False):
    """Mark stacks with zeros and return a logic index vector for modified junctions

    NOTE: Use ONLY for testing the variances (otherwise use analysis.filter in majiq)"""

    minstack = sys.maxint #the minimum value marked as stack
    numstacks = 0
    filtered = [False]*len(junctions)
    for i, junction in enumerate(junctions):
        for j, value in enumerate(junction):
            if value > 0:
                #TODO Use masker, and marking stacks will probably be faster.
                copy_junc = list(junction)
                copy_junc.pop(j)
                copy_junc = array(copy_junc)
                copy_junc = copy_junc[copy_junc > 0]
                #FINISH TODO
                mean_rest = np.mean(copy_junc)
                pval = polyfitnb.get_negbinom_pval(fitted_1_r, mean_rest, value)
                if pval < pvalue_limit:
                    junctions[i, j] = -2
                    minstack = min(minstack, value)
                    numstacks += 1
                    filtered[i] = True
        masked_less(junction, 0)
    if logger: logger.info("Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"%(junctions.size, numstacks, pvalue_limit, (float(numstacks)/junctions.size)*100))

    return filtered


def plot_stacks_method1Vsmethod2(stacks_data, score1_name, score2_name, replica1_name, replica2_name, plotpath=None):
    """Compute 2 plots: one for the variance, one for the mean"""

    def autolabel(rects, direction=1):
        # attach some text labels
        pos_number=1.02
        for rect in rects:
            if direction<0:
                pos_number=1.06
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., pos_number*height*direction, '%.2f' % height, ha='center', va='bottom', size='x-small')

    stats_names = ['mean', 'variance']
    for name_i, stat_name in enumerate(stats_names):
        better = []
        worse = []
        plotname="Majiq_stack_analysis_%s_%sVs%s" % (stat_name, score1_name, score2_name)

        for pval in sorted(stacks_data.keys(), reverse=True):
            total_junc = len(stacks_data[pval][score1_name][stat_name])
            if total_junc == 0:
                print "[%8s for pval=%.1e]:: Total junctions: %d; Better in %s: %s (%.2f%%); Better in %s: %s (%.2f%%)" % \
                      (stat_name, pval, total_junc, score1_name, 0, 0, score2_name, 0, 0)
                continue

            score1 = stacks_data[pval][score1_name][stat_name]
            score2 = stacks_data[pval][score2_name][stat_name]

            better_in_method1 = sum(score1 < score2)
            better_in_method2 = sum(score1 > score2)

            better.append((float(better_in_method1)/total_junc)*100)
            worse.append((float(better_in_method2)/total_junc)*100)

            print "[%8s for pval=%.1e]:: Total junctions: %d; Better in %s: %s (%.2f%%); Better in %s: %s (%.2f%%)" % \
                  (stat_name, pval, total_junc, score1_name, better_in_method1, (float(better_in_method1)/total_junc)*100, score2_name, better_in_method2, (float(better_in_method2)/total_junc)*100)

        ind = np.arange(len(better))    # the x locations for the groups
        width = 0.35                    # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, better, width, color='g')
        rects2 = ax.bar(ind, [-val for val in worse], width, color='r')
        rects3 = ax.bar(ind+width, add(better, [-val for val in worse]), width, color='b')

        # add some
        ax.set_ylabel('Percentage of junctions')
        ax.set_title('Threshold p-value')
        ax.set_xticks(ind+width)
        ax.set_xticklabels(sorted(stacks_data.keys(), reverse=True), size='x-small')

        ax.legend((rects1[0], rects2[0], rects3[0]), ('improved', 'worsen', 'difference'),  prop={'size': 6})
        autolabel(rects1)
        autolabel(rects2, direction=-1)
        autolabel(rects3)

        suptitle(plotname)

        _save_or_show(plotpath, plotname)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('par1', help='Path for parameters of replica1')
    parser.add_argument('par2', help='Path for parameters of replica2')  
    parser.add_argument('--k', default=50, type=int, help='Number of positions to sample per iteration')
    parser.add_argument('--m', default=100, type=int, help='Number of samples')     
    parser.add_argument('--junctype', default='rand10k', help='The type of junction to analyze. (Inc, Exc or rand10k for now)')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    parser.add_argument('--output', default=None, help="Path to save the results to.")
    parser.add_argument('--minnonzero', default=5, type=int, help='Minimum number of positive positions to consider the junction')
    parser.add_argument('--maxnonzero', default=1000000, type=int, help='Maximum number of positive positions to consider the junction')
    parser.add_argument('--minreads', default=10, type=int, help='Minimum number of reads combining all positions in a junction to be considered')
    parser.add_argument('--dispersion', default=0.1, type=float, help='Dispersion factor (used in junctions sampling).')

    args = parser.parse_args()

    #Get the experiment names
    rep1_name = os.path.basename(args.par1).split('.')[-2]
    rep2_name = os.path.basename(args.par2).split('.')[-2]

    replica1, replica2, fit_func1, fit_func2, rep1_gc, rep2_gc, const_rep1, const_rep2 = junction_sample.load_junctions(args.par1, args.par2, args, fromlsv=True)

    # methods = {
    #     'Poisson':                  {'discardzeros': 0, 'trimborder': False,   'nb': False},
    #     'Naive_Boots':              {'discardzeros': 0, 'trimborder': False,   'nb': False},
    #     'Naive_Boots_trim_borders': {'discardzeros': 1, 'trimborder': 5,    'nb': False},
    #     'Naive_Boots_no_zeros':     {'discardzeros': 1, 'trimborder': False,   'nb': False},
    #     'Neg_Binomial':             {'discardzeros': 0, 'trimborder': False,   'nb': True},
    #     'Majiq':                    {'discardzeros': 1, 'trimborder': 5,    'nb': True},
    #     'Majiq_with_zeros':         {'discardzeros': 0, 'trimborder': 5,    'nb': True},
    #     'Majiq_no_stacks':          {'discardzeros': 1, 'trimborder': 5,    'nb': True}
    #     'Majiq_padding_5':          {'discardzeros': 5, 'trimborder': 5,    'nb': True},
    #     'Majiq_padding_10':         {'discardzeros': 10,'trimborder': 5,    'nb': True},
    #     'Majiq_gc_norm':            {'discardzeros': 1, 'trimborder': 5,    'nb': True},
    #
    # }
    #
    # scores_cached = {}
    #
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_10','Majiq', args.plotpath, scores_cached)
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_10','Majiq_with_zeros', args.plotpath, scores_cached)
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_5', 'Majiq', args.plotpath, scores_cached)
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_5', 'Majiq_with_zeros', args.plotpath, scores_cached)
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_10','Majiq_padding_5', args.plotpath, scores_cached)
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Poisson',         'Majiq', args.plotpath, scores_cached, coverage=True )
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Poisson',         'Majiq_padding_5', args.plotpath, scores_cached )
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Poisson',         'Majiq_padding_10', args.plotpath, scores_cached )
    #
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',         'Majiq', args.plotpath, scores_cached)
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots_no_zeros','Majiq', args.plotpath, scores_cached)
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_with_zeros',    'Majiq', args.plotpath, scores_cached)
    #
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',     'Neg_Binomial', args.plotpath, scores_cached)
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',     'Naive_Boots_trim_borders', args.plotpath, scores_cached)
    #
    # compare_methods(args.m, args.k, np.array(replica1), np.array(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_gc_norm',      'Majiq', args.plotpath, scores_cached, rep1_gc=rep1_gc, rep2_gc=rep2_gc)

    methods = {
        'Majiq':                    {'discardzeros': 1, 'trimborder': 5,    'nb': True},
        'Majiq_no_stacks':          {'discardzeros': 1, 'trimborder': 5,    'nb': True}
    }

    # const_rep1, const_rep2 = intersect_const_juncs(const_rep1, const_rep2)
    #
    # stacks_data = defaultdict(lambda: defaultdict())
    # pvals_stacks = [1.0/10**power for power in (3, 5, 7, 9)]
    #
    # for pval_stack in pvals_stacks:
    #     rep1_fresh = np.array(const_rep1)
    #     rep2_fresh = np.array(const_rep2)
    #     stacks_filtered_rep1 = mark_stacks(rep1_fresh, fit_func1, pval_stack, args.dispersion)
    #     stacks_filtered_rep2 = mark_stacks(rep2_fresh, fit_func2, pval_stack, args.dispersion)
    #     filtered_all = map(lambda a, b: a or b, stacks_filtered_rep1, stacks_filtered_rep2)
    #
    #     rep1_filtered = rep1_fresh[np.array(filtered_all)]
    #     rep2_filtered = rep2_fresh[np.array(filtered_all)]
    #
    #     # compare_methods(args.m, args.k, rep1_filtered, rep2_filtered, rep1_name, rep2_name, fit_func1, fit_func2, methods, method_stacks_name, 'Majiq', args.plotpath, scores_cached)  NOT USED!!
    #
    #     mean_method1_rep1, var_method1_rep1, samples_not_used = junction_sample.sample_from_junctions(rep1_filtered, args.m, args.k, fit_func=fit_func1, poisson=False, **methods['Majiq_no_stacks'])
    #     mean_method1_rep2, var_method1_rep2, samples_not_used = junction_sample.sample_from_junctions(rep2_filtered, args.m, args.k, fit_func=fit_func2, poisson=False, **methods['Majiq_no_stacks'])
    #     score_method1 = calc_score_tmp(mean_method1_rep1, var_method1_rep1, mean_method1_rep2, var_method1_rep2)
    #     mean_method2_rep1, var_method2_rep1, samples_not_used = junction_sample.sample_from_junctions(const_rep1[np.array(filtered_all)], args.m, args.k, fit_func=fit_func1, poisson=False, **methods['Majiq'])
    #     mean_method2_rep2, var_method2_rep2, samples_not_used = junction_sample.sample_from_junctions(const_rep2[np.array(filtered_all)], args.m, args.k, fit_func=fit_func2, poisson=False, **methods['Majiq'])
    #     score_method2 = calc_score_tmp(mean_method2_rep1, var_method2_rep1, mean_method2_rep2, var_method2_rep2)
    #
    #     stacks_data[pval_stack]['Majiq_no_stacks'] = score_method1
    #     stacks_data[pval_stack]['Majiq'] = score_method2
    #
    #     plot_method1Vsmethod2(score_method1, score_method2, 'Majiq_no_stacks_%.10f' % pval_stack,  'Majiq_%.10f' % pval_stack, rep1_name, rep2_name, args.plotpath)
    #
    # plot_stacks_method1Vsmethod2(stacks_data, 'Majiq_no_stacks',  'Majiq', rep1_name, rep2_name, args.plotpath)

    replica_quan1, info_quan1 = filter.quantifiable_in_group([replica1], args.minnonzero, args.minreads, None)
    replica_quan2, info_quan2 = filter.quantifiable_in_group([replica2], args.minnonzero, args.minreads, None)
    lreps_quan, linfos_quan = filter.lsv_intersection([replica_quan1, info_quan1], [replica_quan2, info_quan2])
    
    # Majiq with stacks
    stacks_data = defaultdict(lambda: defaultdict())

    pvals_stacks = [1.0/10**power for power in (3, 5, 7, 9)]

    for pval_stack in pvals_stacks:
        print "Computing Majiq comparison for stack removal, p-value=%.10f" %pval_stack

        stacks_filtered_rep1, junc_filt_rep1 = lsv_mark_stacks(np.array(replica1), fit_func1, pval_stack, .1, logger=None)
        stacks_filtered_rep2, junc_filt_rep2 = lsv_mark_stacks(np.array(replica2), fit_func2, pval_stack, .1, logger=None)

        filtered_lsv1, info_filt1 = filter.quantifiable_in_group([stacks_filtered_rep1], args.minnonzero, args.minreads, None)
        filtered_lsv2, info_filt2 = filter.quantifiable_in_group([stacks_filtered_rep2], args.minnonzero, args.minreads, None)

        lreps_filtered, linfos_filtered = filter.lsv_intersection([filtered_lsv1, info_filt1], [filtered_lsv2, info_filt2])

        lreps, linfos = majiq_intersec([lreps_filtered, linfos_filtered], [lreps_quan, linfos_quan])

        # compare_methods(args.m, args.k, rep1_filtered, rep2_filtered, rep1_name, rep2_name, fit_func1, fit_func2, methods, method_stacks_name, 'Majiq', args.plotpath, scores_cached)  NOT USED!!

        mean_method1_rep1, var_method1_rep1, samples_not_used = junction_sample.sample_from_junctions(only_juncs_stacked(lreps[0][0], linfos, junc_filt_rep1, junc_filt_rep2), args.m, args.k, fit_func=fit_func1, poisson=False, **methods['Majiq_no_stacks'])
        mean_method1_rep2, var_method1_rep2, samples_not_used = junction_sample.sample_from_junctions(only_juncs_stacked(lreps[0][1], linfos, junc_filt_rep1, junc_filt_rep2), args.m, args.k, fit_func=fit_func2, poisson=False, **methods['Majiq_no_stacks'])
        # mean_method1_rep1, var_method1_rep1, samples_not_used = junction_sample.sample_from_junctions([e for j in lreps[0][0] for e in j], args.m, args.k, fit_func=fit_func1, poisson=False, **methods['Majiq_no_stacks'])
        # mean_method1_rep2, var_method1_rep2, samples_not_used = junction_sample.sample_from_junctions([e for j in lreps[0][1] for e in j], args.m, args.k, fit_func=fit_func2, poisson=False, **methods['Majiq_no_stacks'])
        score_method1 = calc_score_tmp(mean_method1_rep1, var_method1_rep1, mean_method1_rep2, var_method1_rep2)

        mean_method2_rep1, var_method2_rep1, samples_not_used = junction_sample.sample_from_junctions(only_juncs_stacked(lreps[1][0], linfos, junc_filt_rep1, junc_filt_rep2), args.m, args.k, fit_func=fit_func1, poisson=False, **methods['Majiq_no_stacks'])
        mean_method2_rep2, var_method2_rep2, samples_not_used = junction_sample.sample_from_junctions(only_juncs_stacked(lreps[1][1], linfos, junc_filt_rep1, junc_filt_rep2), args.m, args.k, fit_func=fit_func2, poisson=False, **methods['Majiq_no_stacks'])

        # mean_method2_rep1, var_method2_rep1, samples_not_used = junction_sample.sample_from_junctions([e for j in lreps[1][0] for e in j], args.m, args.k, fit_func=fit_func1, poisson=False, **methods['Majiq'])
        # mean_method2_rep2, var_method2_rep2, samples_not_used = junction_sample.sample_from_junctions([e for j in lreps[1][1] for e in j], args.m, args.k, fit_func=fit_func2, poisson=False, **methods['Majiq'])
        score_method2 = calc_score_tmp(mean_method2_rep1, var_method2_rep1, mean_method2_rep2, var_method2_rep2)

        stacks_data[pval_stack]['Majiq_no_stacks'] = score_method1
        stacks_data[pval_stack]['Majiq'] = score_method2

    plot_stacks_method1Vsmethod2(stacks_data, 'Majiq_no_stacks',  'Majiq', rep1_name, rep2_name, args.plotpath)


def only_juncs_stacked(lreps, linfos, junc_filt_rep1, junc_filt_rep2):
    juns = []
    for ii, jun in enumerate(lreps):
        for jj, e in enumerate(jun):
            if junc_filt_rep1[linfos[ii][1]] == jj or junc_filt_rep2[linfos[ii][1]] == jj:
                juns.append(e)
    return juns

def majiq_intersec(lsv_list1, lsv_list2):

    lsv_match = [[[], []], [[], []]]
    match_info = []

    ids1 = set([xx[1] for xx in lsv_list1[1]])
    ids2 = set([xx[1] for xx in lsv_list2[1]])
    matched_names = ids1.intersection(ids2)

    for ii in matched_names:
        for idx, nm in enumerate(lsv_list1[1]):
            if nm[1] == ii:
                lsv_match[0][0].append(lsv_list1[0][0][idx][0])
                lsv_match[0][1].append(lsv_list1[0][1][idx][0])
                match_info.append(nm)
                break
        for idx, nm in enumerate(lsv_list2[1]):
            if nm[1] == ii:
                lsv_match[1][0].append(lsv_list2[0][0][idx][0])
                lsv_match[1][1].append(lsv_list2[0][1][idx][0])
                break

    return lsv_match, match_info


def intersect_const_juncs(const_rep1, const_rep2):
    lsv_match = [[], []]

    ids1 = set(const_rep1[1])
    ids2 = set(const_rep2[1])
    matched_names = ids1.intersection(ids2)

    for ii in matched_names:
        for idx, nm in enumerate(const_rep1[1]):
            if nm == ii:
                lsv_match[0].append(const_rep1[0][idx])
                break
        for idx, nm in enumerate(const_rep2[1]):
            if nm == ii:
                lsv_match[1].append(const_rep2[0][idx])
                break

    return np.array(lsv_match[0]), np.array(lsv_match[1])



if __name__ == '__main__':
    main()
