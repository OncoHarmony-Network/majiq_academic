import argparse
import pickle

from pylab import *

def _save_or_show(plotpath):
    """Generic function that either shows in a popup or saves the figure, depending if a flag was given"""
    if plotpath:
        savefig("%s.png"%(plotpath), bbox_inches='tight')
        clf()
    else:
        show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('delta1', help='Delta PSI values for comparison 1.')
    parser.add_argument('delta2', help='Delta PSI values for comparison 2.')
    parser.add_argument('--name1', help='Name of method 1.')
    parser.add_argument('--name2', help='Name of method 2.')
    parser.add_argument('--compname', default="", help='Name of comparison made (Ex: tissue1 VS tissue2).')        
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    args = parser.parse_args()

    delta1 = pickle.load(open(args.delta1))
    delta2 = pickle.load(open(args.delta2))
    print args.name1
    print delta1
    print
    print args.name2
    print delta2
    print
    max_value = 1
    title(args.compname)
    xlabel("%s %s"%(args.name1, args.compname))
    ylabel("%s %s"%(args.name2, args.compname))
    xlim(0, max_value)
    ylim(0, max_value)
    plot([0, max_value], [0, max_value])
    plot(delta1, delta2, '.r')
    _save_or_show(args.plotpath)

if __name__ == '__main__':
    main()



