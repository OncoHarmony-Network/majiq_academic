import argparse
import pickle

from pylab import *

def _save_or_show(plotpath, name):
    if plotpath:

        savefig("%s_%s.png"%(plotpath, name), width=300, height=300, dpi=100)
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
    parser.add_argument('--max',  default=1000, type = int, help="Max number of events to analyze")
    parser.add_argument('--step', default=1, type = int, help="Steps of best events to take")
    parser.add_argument('--colors', nargs='+',  default = ["blue", "green", "red"], help="Steps of best events to take")
    args = parser.parse_args()

    fig = figure(figsize=[7, 7])
    #figure out how many groups of events exist
    numbins = int(round(args.max/args.step))
    #plot the diagonal
    diagonal = linspace(0, 1, num=numbins+1)
    plot(range(0, numbins+1), diagonal, '--', color="#cccccc") 


    #plot the lines of all ratios (pairs of ranks)
    for i, ratio in enumerate(args.ratios):
        ratio = array(pickle.load(open(ratio))) 
        xlabel("Events (ranked)", fontsize=20)
        ylabel("Ratio (total %s events)"%args.max, fontsize=20)
        xlim(0, numbins)
        ylim(0, 1) #a ratio
        linetype = '-'
        if args.fdr:
            if args.fdr[i]: 
                linetype = '--'


        plot(range(0, numbins+1), ratio, linetype, label=args.labels[i], linewidth=2, color=args.colors[i])


    title(args.title, fontsize=16)    
    legend(loc=2, fontsize=12)
    _save_or_show(args.plotpath, "rankcomp")


if __name__ == '__main__':
    main()
