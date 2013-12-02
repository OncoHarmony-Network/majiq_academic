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
    parser.add_argument('delta1', help='Path with matfile with replica1 and replica2')    
    parser.add_argument('delta2', help='Path for parameters of replica1')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--name1')
    parser.add_argument('--name2')   
    args = parser.parse_args()

    delta1 = pickle.load(open(args.delta1))
    delta2 = pickle.load(open(args.delta2))
    title("%s vs %s"%(args.name1, args.name2))
    max_value = max(max(delta1), max(delta2))*0.7
    xlabel(args.name1)
    ylabel(args.name2)
    xlim(0, max_value)
    ylim(0, max_value)
    plot([0, max_value], [0, max_value])
    plot(delta1, delta2, '.r')
    _save_or_show(args.plotpath)



if __name__ == '__main__':
    main()
