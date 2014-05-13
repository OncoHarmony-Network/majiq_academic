from matplotlib import use
use('Agg', warn=False)

from scipy.stats import pearsonr
import argparse
from pylab import *
import analysis.filter
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


def plot_method1Vsmethod2_by_coverage(scores1, scores2, score1_name, score2_name, replica1_name, replica2_name, plotpath=None):
    """Compute 6 plots, one for each coverage subset, variance and mean"""

    # TODO: THIS IS NOT WORKING PROPERLY YET. Review why scores1 has zeros
    plotname="%sVs%s-Coverage" % (score1_name, score2_name)

    print scores1

    stats_names = ['mean', 'variance']

    f, axarr = subplots(len(stats_names), len(scores1), sharex=True, sharey=True)
    axarr[0][0].set_ylabel(score2_name)
    axarr[0][-1].set_xlabel(score1_name)


    for score_idx in range(len(scores1)):
        total_junctions = float(len(scores1[score_idx]['mean']))  # The same in all mean/variance scores1/scores2

        for name_i, stat_name in enumerate(stats_names):
            axarr[name_i][0].set_title(stat_name, fontsize=10)
            score1, score2 = scores1[score_idx][stat_name], scores2[score_idx][stat_name]

            print score1

            print "\n\n"
            print score2[0]

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
            axarr[name_i][score_idx].text(abs(max_value)*0.1, max_value-abs(max_value)*0.2, "(%.2f%%)" % (better_in_method1/total_junctions)*100, fontsize=12, bbox={'facecolor':'yellow', 'alpha':0.3, 'pad':10})
            axarr[name_i][score_idx].text(abs(max_value)*0.8, max_value-abs(max_value)*0.8, "(%.2f%%)" % (better_in_method2/total_junctions)*100, fontsize=12, bbox={'facecolor':'yellow', 'alpha':0.3, 'pad':10})

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


def compare_methods(m, k, replica1, replica2, rep1_name, rep2_name, fit_func1, fit_func2, methods, method1_name, method2_name, plotpath, scores_cached, coverage=False):

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


    mean_method1_rep1, var_method1_rep1, samples_not_used = junction_sample.sample_from_junctions(replica1, m, k, fit_func=fit_func1, poisson=method1_name=='Poisson', **methods[method1_name])
    mean_method1_rep2, var_method1_rep2, samples_not_used = junction_sample.sample_from_junctions(replica2, m, k, fit_func=fit_func2, poisson=method1_name=='Poisson', **methods[method1_name])

    if method1_name in scores_cached:
        score_method1 = scores_cached[method1_name]
    else:
        score_method1 = calc_score_tmp(mean_method1_rep1, var_method1_rep1, mean_method1_rep2, var_method1_rep2)
        scores_cached[method1_name] = score_method1

    mean_method2_rep1, var_method2_rep1, samples_not_used = junction_sample.sample_from_junctions(replica1, m, k, fit_func=fit_func1, poisson=method2_name=='Poisson', **methods[method2_name])
    mean_method2_rep2, var_method2_rep2, samples_not_used = junction_sample.sample_from_junctions(replica2, m, k, fit_func=fit_func2, poisson=method2_name=='Poisson', **methods[method2_name])

    if method2_name in scores_cached:
        score_method2 = scores_cached[method2_name]
    else:
        score_method2 = calc_score_tmp(mean_method2_rep1, var_method2_rep1, mean_method2_rep2, var_method2_rep2)
        scores_cached[method2_name] = score_method2

    plot_method1Vsmethod2(score_method1, score_method2, method1_name,  method2_name, rep1_name, rep2_name, plotpath)

    mean_var_by_replicate = {
            'rep1': {
                'name'  : rep1_name,
                'mean1' : mean_method1_rep1,
                'var1'  : var_method1_rep1,
                'mean2' : mean_method2_rep1,
                'var2'  : var_method2_rep1
            },
            'rep2': {
                'name'  : rep2_name,
                'mean1' : mean_method1_rep2,
                'var1'  : var_method1_rep2,
                'mean2' : mean_method2_rep2,
                'var2'  : var_method2_rep2
            }
        }

    plot_sigmaVsMu_met1_met2(mean_var_by_replicate, score_method1, score_method2, method1_name,  method2_name, plotpath)


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

    #Get the experiment names
    rep1_name = os.path.basename(args.par1).split('.')[-2]
    rep2_name = os.path.basename(args.par2).split('.')[-2]

    replica1, replica2, fit_func1, fit_func2 = junction_sample.load_junctions(args.par1, args.par2, args, fromlsv=False)

    methods = {
        'Poisson':                  {'discardzeros': 0, 'trimborder': False,   'nb': False},
        'Naive_Boots':              {'discardzeros': 0, 'trimborder': False,   'nb': False},
        'Naive_Boots_trim_borders': {'discardzeros': 1, 'trimborder': True,    'nb': False},
        'Naive_Boots_no_zeros':     {'discardzeros': 1, 'trimborder': False,   'nb': False},
        'Neg_Binomial':             {'discardzeros': 0, 'trimborder': False,   'nb': True},
        'Majiq':                    {'discardzeros': 1, 'trimborder': True,    'nb': True},
        'Majiq_with_zeros':         {'discardzeros': 0, 'trimborder': True,    'nb': True},
        'Majiq_no_stacks_0001':     {'discardzeros': 1, 'trimborder': True,    'nb': True},
        'Majiq_no_stacks_0005':     {'discardzeros': 1, 'trimborder': True,    'nb': True},
        'Majiq_no_stacks_001':      {'discardzeros': 1, 'trimborder': True,    'nb': True},
        'Majiq_padding_5':          {'discardzeros': 5, 'trimborder': True,    'nb': True},
        'Majiq_padding_10':         {'discardzeros': 10,'trimborder': True,    'nb': True},
    }

    scores_cached = {}

    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_10',      'Majiq', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_10',      'Majiq_with_zeros', args.plotpath, scores_cached)

    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_5',      'Majiq', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_5',      'Majiq_with_zeros', args.plotpath, scores_cached)

    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_padding_10',      'Majiq_padding_5', args.plotpath, scores_cached)

    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Poisson',       'Majiq', args.plotpath, scores_cached )
    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Poisson',       'Majiq_padding_5', args.plotpath, scores_cached )
    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Poisson',       'Majiq_padding_10', args.plotpath, scores_cached )

    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',   'Majiq', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots_no_zeros',  'Majiq', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_with_zeros',      'Majiq', args.plotpath, scores_cached)

    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',   'Neg_Binomial', args.plotpath, scores_cached)
    compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Naive_Boots',   'Naive_Boots_trim_borders', args.plotpath, scores_cached)


    # Majiq with stacks
    marks_stacks = [0.0001, 0.0005, 0.001]  # Note that this array has to be sorted in ascending order
    for mark_stacks in marks_stacks:
        analysis.filter.mark_stacks(replica1, fit_func1, mark_stacks, args.dispersion)
        analysis.filter.mark_stacks(replica2, fit_func2, mark_stacks, args.dispersion)

        compare_methods(args.m, args.k, list(replica1), list(replica2), rep1_name, rep2_name, fit_func1, fit_func2, methods, 'Majiq_no_stacks_'+str(mark_stacks).split('.')[1], 'Majiq', args.plotpath, scores_cached)


if __name__ == '__main__':
    main()
