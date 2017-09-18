from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
from majiq.src.config import Config

import os
import numpy as np

__author__ = 'jordi'


def save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if a flag was given"""
    if plotname:
        plt.title(plotname)

    if plotpath:
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)
        plt.savefig("%s/%s.png" % (plotpath, plotname.replace(" ", "_")), bbox_inches='tight')
        # WNo spaces allowed, underscores!
        plt.clf()
    else:
        plt.show()


def plot_mappability_zeros(junctions, plotpath, numzeros, plotname):
    """
    Plot the percentage of zero positions in the junction sites against the total sum of the positions in the junction.
    This should be proportional. Dots in the right mean that there are some zero positions in the junctions with high
    coverage, meaning that there is probably some mappability bias.
    """
    junction_sum = junctions.sum(axis=1)
    # Cumulative positional value per junction
    plt.plot(junction_sum, numzeros / float(junctions.shape[1]), '*')
    save_or_show(plotpath, plotname)


def plot_negbinomial_fit(mean_junc, std_junc, fit_function, plotpath, plotname):
    # plot the fit of the line
    plt.xlabel("Mean")
    plt.ylabel("Std")
    plt.xlim(0, max(mean_junc) * 1.1)  # adjust the x axis of the plot
    plt.ylim(0, max(std_junc) * 1.1)  # adjust the y axis of the plot
    plt.plot(mean_junc, std_junc, '*')  # the mean and std for e very junction that passes the filter
    plt.plot(mean_junc, fit_function(mean_junc), '-r')
    save_or_show(plotpath, plotname)


def plot_fitting(ecdf, plotpath, extra=[], title='', title_extra=[], plotname=None):
    if plotpath:
        plt.xlabel("P-value")
        plt.ylabel("non_corrected ECDF")
        plt.plot(np.linspace(0, 1, num=len(ecdf)), ecdf, label='Fitted NB')
        for ex_idx, extraval in enumerate(extra):
            plt.plot(np.linspace(0, 1, num=len(extraval)), extraval, label=title_extra[ex_idx])
        plt.plot([0, 1], 'k')
        plt.legend(loc='upper left')
        if plotname is None:
            fname = title
        else:
            fname = plotname
        save_or_show(plotpath, fname)




def plot_gc_content():
    majiq_config = Config()
    idx = 0
    for tissue, list_idx in majiq_config.tissue_repl.items():
        plt.figure(idx)
        for exp_n in list_idx:
#            f = interpolate.interp1d(mglobals.gc_means[exp_n], mglobals.gc_bins_vaL[exp_n])
#            print mglobals.gc_means[exp_n]
            mn = majiq_config.gc_means[exp_n].min()
            mx = majiq_config.gc_means[exp_n].max()
            xx = np.arange(mn, mx, 0.001)
            yy = majiq_config.gc_factor[exp_n](xx)
            # print "XX ",exp_n, xx
            # print "Yy", exp_n, yy
            plt.plot(xx, yy, label=majiq_config.exp_list[exp_n])
            plt.axis((0.3, 0.7, 0.5, 1.5))
            plt.title("Gc factor")
            plt.grid()
            plt.legend(loc='upper left')
#        pyplot.show()
        plt.savefig('%s/gcontent_%s.png' % (majiq_config.outDir, tissue))
        idx += 1


def plot_pearsoncorr(var1, var2, my_title, my_xlabel, my_ylabel, plotpath=None, max_value=None):
    var1 = np.array(var1)
    var2 = np.array(var2)
    plt.xlabel(my_xlabel)
    plt.ylabel(my_ylabel)

    if not max_value:
        max_value = max(max(var1), max(var2))

    plt.xlim(0, max_value)
    plt.ylim(0, max_value)

    # plot([0, max_value], [0, max_value])
    pear, pvalue = pearsonr(var1, var2)
    r_squared = pear ** 2
    a, b = np.polyfit(var1, var2, 1)
    fit_func = np.poly1d([a, b])
    plt.plot(var1, fit_func(var1), '-r')
    # percentage_under = sum(var1 < var2)/float(len(var1))
    plt.text(abs(max_value) * 0.1, max_value - abs(max_value) * 0.2,
             r'$R^2$: %.2f (p-value: %.2E)' % (r_squared, pvalue),
             fontsize=18, bbox={'facecolor': 'yellow', 'alpha': 0.3, 'pad': 10})
    plt.title(my_title)
    print(r"%s R^2: %.2f (p-value: %.2E)" % (my_title, r_squared, pvalue))
    plt.plot(var1, var2, '.')
    if plotpath:
        save_or_show(plotpath, my_title)


def plot_matrix(matrix, my_title, plotname, plotpath):
    plt.clf()
    ax = plt.subplot(1, 1, 1)
    plt.title(my_title)
    plt.imshow(matrix)
    plt.xlabel(u"PSI i")
    plt.ylabel(u"PSI j")
    ax.set_xticklabels([0, 0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0, 0, 0.25, 0.5, 0.75, 1])

    _save_or_show(plotpath, plotname=plotname)


def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        plt.savefig("%s%s.png" % (plotpath, plotname), bbox_inches='tight')
        plt.clf()
    else:
        plt.show()