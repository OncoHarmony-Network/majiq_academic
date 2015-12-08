import matplotlib
matplotlib.use('Agg')
import argparse
import numpy as np
import matplotlib.pyplot as plt
import cPickle as pkl
import sys
import colorbrewer as cb
import math
from scripts.utils import save_or_show
from voila.utils import utils_voila as uvoila


def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb


def get_n1(n1_files, rep_name):
    for n1_file in n1_files:
        if rep_name in n1_file:
            return pkl.load(open(n1_file))[0]


def add_method(method_name, method_num, targetdir, group_folder):
    print "Processing %s ..." % method_name
    ratios = uvoila.list_files_or_dir([targetdir], prefix='ratios', containing=method_name)
    n1s = uvoila.list_files_or_dir([targetdir], prefix='n1', containing=method_name)
    rep_names = [fie.split('/')[-2] for fie in ratios]

    bins_area=200
    bins_all=[]
    bins_n1=[]
    n_list=[]

    # plot the lines of all ratios (pairs of ranks)
    for i, ratio_path in enumerate(ratios):
        alpha_plt=0.4
        ratio = np.array(pkl.load(open(ratio_path)))
        numevents = ratio.shape[0]
        if group_folder in rep_names[i]:
            #alpha_plt=.9
            #plot the diagonal if we are in the first step
            diagonal = np.linspace(0, 1, num=numevents+1)
            plt.plot(diagonal, diagonal, '--', color="#cccccc")
            plt.plot(np.linspace(0, 1, len(ratio)), ratio, '-', label="%s Groups - %.0f%%" % (method_name.upper(), math.ceil(ratio[-1]*100)) , linewidth=4, color=rgb_to_hex(cb.Paired[10][2*method_num + 1]), alpha=1)
            continue

        n_list.append(len(ratio))
        aux=[]
        for aa in xrange(0,bins_area,1):
            aux.append(ratio[int((aa/(bins_area*1.0))*len(ratio))])
        bins_all.append(aux)

        n1 = get_n1(n1s, rep_names[i])
        bins_n1.append(np.arange(0, len(ratio)*1./n1, 1./n1)[:len(ratio)])

    print "[%s] Avg. N=%.2f" % (method_name, np.mean(n_list))

    ndbins=np.array(bins_all, ndmin=2)
    nbins_n1=np.array(bins_n1)

    plt.plot(np.arange(0+1./(2*bins_area), 1, 1./bins_area), np.mean(ndbins, axis=0), label="%s Pairs mean & STD - %.0f%%" % (method_name.upper(), math.ceil(np.mean(ndbins, axis=0)[-1]*100)), color=rgb_to_hex(cb.Paired[10][2*method_num + 1]), lw=2)
   # plt.plot(np.arange(0, 1, 1./nbins_n1.shape[1]), np.mean(nbins_n1, axis=0), '--', label="%s Pairs Null Model & STD - %.0f%%" % (method_name.upper(), math.ceil(np.mean(nbins_n1, axis=0)[-1]*100)), color=rgb_to_hex(cb.Paired[10][2*method_num + 1]), lw=2)

    plt.fill_between(np.arange(0+1./(2*bins_area),1,1./bins_area),
        np.mean(ndbins, axis=0)-np.std(ndbins, axis=0),
        np.mean(ndbins, axis=0)+np.std(ndbins, axis=0), facecolor=rgb_to_hex(cb.Paired[10][2*method_num]), lw=0, alpha=alpha_plt)

    return np.mean(n_list)


def main():
    """
    Delta PSI reproducibility plot. Comparison of how reproducible are the rankings of the most confident
    changing/non-changing events between two pairs of conditions.
    """
    methods=[
        'majiq',
        'miso',
        'mats',
        'naive',
	'nondenovo'
    ]

    parser = argparse.ArgumentParser(description="")
    parser.add_argument('targetdir', type=str, help='Folder containing the ratios (Pickle files) for each method.')
    parser.add_argument('--plotpath', required=True, help='Output file path.')
    parser.add_argument('--plotname', default='figure2b', help='Filename of the plot.')
    parser.add_argument('--plottitle', default='dPSI Reproducibility - H124L123 Vs H56L45', help='Title of the plot.')
    parser.add_argument('--groupfolder', default='H124L123H56L45', help='Group subfolder name with the ratios for delta psi groups. [Only if there are group comparisons].')
    parser.add_argument('--extension', choices=('pdf', 'png'), default='pdf', help='Format of the plot generated.')
    parser.add_argument('--methods', dest='methods_set', nargs='*', default=methods, help='Methods to be compared. It could be any combination of: %s.' % ', '.join(methods))
    args = parser.parse_args()

    if not len(set(methods).intersection(set(args.methods_set))):
        print 'Please choose between the list of methods currently available to compare: [%s].' % ', '.join(methods)
        sys.exit(1)

    font = {'size': 16} #here also 'weight' and 'family'
    matplotlib.rc('font', **font)

    avg_ns=[]
    plt.figure(figsize=[10, 10]) # In inches
    for met_num, met in enumerate(args.methods_set):
        avg_ns.append(add_method(met, met_num, args.targetdir, args.groupfolder))

    plt.xlabel("Fraction of selected events", fontsize=20)
    plt.ylabel("Fraction of events reproduced", fontsize=20)
    plt.xlim(0, 1)
    plt.ylim(0, 1) #a ratio

    plt.title("%s (N=%d)" % (args.plottitle, int(np.mean(avg_ns))), fontsize=16)
    plt.legend(loc=2, fontsize=11)
    save_or_show(args.plotpath, args.plotname, exten=args.extension)

if __name__ == '__main__':
    main()
