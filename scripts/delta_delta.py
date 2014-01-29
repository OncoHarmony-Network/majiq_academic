import argparse
import pickle

from pylab import *


def _save_or_show(plotpath):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig(plotpath, bbox_inches='tight', width=500, height=500, dpi=100) 
        clf()
    else:
        show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('delta1', help='Path for delta1')    
    parser.add_argument('delta2', help='Path for delta2')
    parser.add_argument('--filterindex2', help="Index of the junctions that were filtered on delta2")
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--name1')
    parser.add_argument('--name2')   
    args = parser.parse_args()

    delta1 = pickle.load(open(args.delta1))
    delta2 = pickle.load(open(args.delta2))

    if args.filterindex2:
        index2 = pickle.load(open(args.filterindex2))
        #apply filter
        filter_delta1 = []
        for i in range(delta1.shape[0]):
            if i in index2:
                filter_delta1.append(delta1[i])

        delta1 = array(filter_delta1)

    total_delta = float(len(delta1))
    min_distance = 0.01 
    better_in_delta1 = sum(abs(delta1)+min_distance < abs(delta2))
    better_in_delta2 = sum(abs(delta1) > abs(delta2)+min_distance)
    deltadelta_score = sum(abs(delta1) - abs(delta2))/total_delta
    print "\nBetter in %s: %s (%.2f%%) Better in %s: %s (%.2f%%)"%(args.name1, better_in_delta1, (better_in_delta1/total_delta)*100, args.name2, better_in_delta2, (better_in_delta2/total_delta)*100)
    print "\nScore = %.5f"%(float(better_in_delta1)/better_in_delta2)
    title("%s vs %s"%(args.name1, args.name2))
    max_value = max(max(delta1), max(delta2))*0.7
    xlabel(args.name1)
    ylabel(args.name2)
    xlim(0, max_value)
    ylim(0, max_value)
    plot([0, max_value], [0, max_value])
    plot(abs(delta1), abs(delta2), '.r')
    _save_or_show(args.plotpath)


if __name__ == '__main__':
    main()







