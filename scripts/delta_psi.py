import sys
import argparse
import pickle


from pylab import *
from scipy.io import savemat
   
BINS = linspace(0, 1, num=39)

def mean_psi(psi_events):
    "Calculate the mean for every junction"
    ret = []
    for psi_dist in psi_events:    
        ret.append(sum(psi_dist*BINS))

    return array(ret)

#deprecated
def sample_psi(psi_events):
    from numpy.random import choice
    "Get a random point estimate NOTE: REQUIRES python 2.7 for numpy new version"
    ret = []
    for psi_dist in psi_events:
        ret.append(choice(BINS, p=psi_dist))

    return array(ret)


def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname), bbox_inches='tight') 
        clf()
    else:
        show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('psi1', help='Path for pickle with the psi values')
    parser.add_argument('psi2', help='Path for pickle with the psi values')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--title', default=None, help='') 
    parser.add_argument('--output', required=True, help='Path to save the pickle output to.')
    args = parser.parse_args()
    psi_values1 = pickle.load(open(args.psi1))
    psi_values2 = pickle.load(open(args.psi2))
    print "Calculate mean PSI"
    psi1 = mean_psi(psi_values1[:,0])
    psi2 = mean_psi(psi_values2[:,0])

    delta_psi = psi1 - psi2
    print "PSI1 min:%s max:%s\nPSI2 min:%s max:%s\n Delta PSI min:%s max:%s "%(max(psi1), min(psi1), max(psi2), min(psi2),max(delta_psi), min(delta_psi))
    #delta_dist, bin_edges = histogram(delta_psi, range=[0,1], bins=len(delta_psi)/10)
    #plot(bin_edges, delta_dist)
    xlim(-1, 1)
    savemat("%s.mat"%args.plotpath, {"DeltaPSI:" : list(delta_psi)})
    pickle.dump(delta_psi, open("%s_deltapsi.pickle"%(args.output), 'w'))
    hist(delta_psi, bins = 60, histtype='step')

    if args.title:
        title(args.title)

    xlabel("Delta PSI")
    _save_or_show(args.plotpath, "deltapsi")


if __name__ == '__main__':
    main()

