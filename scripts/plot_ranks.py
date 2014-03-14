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
    parser.add_argument('--colors', nargs='*',  default = ["blue", "green", "red"], help="Steps of best events to take")
    args = parser.parse_args()

    fig = figure(figsize=[7, 7])
    #figure out how many groups of events exist

    first = True
    #plot the lines of all ratios (pairs of ranks)
    for i, ratio_path in enumerate(args.ratios):

        ratio = array(pickle.load(open(ratio_path))) 
        numevents = ratio.shape[0]
        if first:
            #plot the diagonal if we are in the first step
            diagonal = linspace(0, 1, num=numevents+1)
            plot(range(0, numevents+1), diagonal, '--', color="#cccccc") 
            first = False

        xlabel("Events (ranked)", fontsize=20)
        ylabel("Ratio (total %s events)"%numevents, fontsize=20)
        xlim(0, numevents)
        ylim(0, 1) #a ratio
        linetype = '-'
        if args.fdr:
            if args.fdr[i]: 
                linetype = '--'

        #label is file path if not specified
        if args.labels:
            my_label = args.labels[i]
        else:
            my_label = ratio_path.split(".pickle")[0].split("/")[-1] #.replace('_', ' V=')

        if args.colors:
            plot(range(numevents), ratio, linetype, label=my_label, linewidth=2, color=args.colors[i])          
        else: 
            plot(range(numevents), ratio, linetype, label=my_label, linewidth=2)




    title(args.title, fontsize=16)    
    legend(loc=2, fontsize=12)
    _save_or_show(args.plotpath, "rankcomp")


if __name__ == '__main__':
    main()
