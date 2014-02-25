import argparse
import pickle

from pylab import *
import matplotlib



def _save_or_show(plotpath):
    if plotpath:
        savefig(plotpath, width=200, height=400, dpi=100)
        clf()
    else:
        show()  

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('v', help='')
    parser.add_argument('--plotpath', help="the plotpath")
    parser.add_argument('--title', help="Plot title")
    args = parser.parse_args()

    font = {'size': 22} #here also 'weight' and 'family'b6 
    matplotlib.rc('font', **font)

    v_values = pickle.load(open(args.v))
    title(args.title)
    xlabel("Ranked Events")
    ylabel("P(Delta PSI > V)")
    plot(v_values)
    _save_or_show(args.plotpath)

if __name__ == '__main__':
    main()
