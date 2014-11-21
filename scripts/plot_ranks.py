import matplotlib
matplotlib.use('Agg')
import argparse
import pickle
import os
from pylab import *


def _save_or_show(plotpath, name):
    if plotpath:
        if os.path.isdir(plotpath):
            plot_base_path, plot_name = plotpath, name
        else:
            plot_base_path, plot_name = os.path.split(plotpath)
            if not os.path.exists(plot_base_path):
                os.makedirs(plot_base_path)
            if not plot_name:
                plot_name = name
        savefig("%s/%s.png"%(plot_base_path, plot_name), width=300, height=300, dpi=100)
        print "Saved in:\n%s/%s" % (plot_base_path, plot_name)

        clf()
    else:
        show()  

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('ratios', nargs='+', help='The ratios (pairs of ranks)')
    parser.add_argument('--plotpath', help='')
    parser.add_argument('--labels', nargs='+', help='The labels for the plot lines of the ratios.')
    parser.add_argument('--title', help='The title of the plot')
    parser.add_argument('--fdr', nargs='+', type=int, help="Determine which plots are FDR lines (1) and which are not (0), and paint them as a dotted line [Example: --fdr 0 1 0 0 1]")
    parser.add_argument('--colors', nargs='*',  help="Steps of best events to take")
    parser.add_argument('--plotname', default='rankcomp', help='Plot name')
    parser.add_argument('--grouppairs', action='store_true', default='False', help='Flag for group vs pairs comparison')
    args = parser.parse_args()

    fig = figure(figsize=[10, 10]) # In inches
    #figure out how many groups of events exist

    font = {'size': 16} #here also 'weight' and 'family'
    matplotlib.rc('font', **font)

    first = True
    title = args.title
    lsizes = []
    #plot the lines of all ratios (pairs of ranks)
    for i, ratio_path in enumerate(args.ratios):

        ratio = array(pickle.load(open(ratio_path))) 
        numevents = ratio.shape[0]
        if first:
            #plot the diagonal if we are in the first step
            diagonal = linspace(0, 1, num=numevents+1)
            plot(diagonal, diagonal, '--', color="#cccccc") 
            first = False


        xlabel("Fraction of selected events", fontsize=20)
        ylabel("fraction of events reproduced", fontsize=20)
        xlim(0, 1)
        ylim(0, 1) #a ratio
        linetype = '-'
        if args.fdr:
            if args.fdr[i]: 
                linetype = '--'

        #label is file path if not specified
        if args.labels:
            my_label = "%s (N=%d)" % (args.labels[i], numevents)
        else:
            my_label = ratio_path.split(".pickle")[0].split("/")[-1]

        x_space = linspace(0, 1, len(ratio))
        if args.colors:
            plot(x_space, ratio, linetype, label=my_label, linewidth=2, color=args.colors[i])          
        else: 
            plot(x_space, ratio, linetype, label=my_label, linewidth=2)

        if args.grouppairs:
            if i == len(args.ratios)-1:
                title += "\nGroup N=%d; Pairs Avg. N=%.2f" % (len(args.ratios), np.mean(lsizes))
            else:
                lsizes.append(numevents)

    title("%s" % (title) , fontsize=16)

    if not args.grouppairs:
        legend(loc=2)
    _save_or_show(plotpath=args.plotpath, name=args.plotname)


if __name__ == '__main__':
    main()
