import argparse
import pickle


from pylab import *
#from numpy import histogram2d


BSIZE = 0.025 #TODO To parameters
BINS = arange(0, 1, BSIZE) # The bins for PSI values. With a BSIZE of 0.025, we have 40 BINS
BINS_CENTER = arange(0+BSIZE/2, 1, BSIZE) #The center of the previous BINS. This is used to calculate the mean value of each bin.


def _save_or_show(plotpath, plotname=None):
    """Generic function that either shows in a popup or saves the figure, depending if the plotpath flag"""
    if plotpath:
        savefig("%s%s.png"%(plotpath, plotname), bbox_inches='tight') 
        clf()
    else:
        show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('prior', help='Pickle for the prior matrix')
    parser.add_argument('posterior', help='Pickle for the posterior matrices')
    parser.add_argument('--plotpath', default=None, help='Path to save the plot to, if not provided will show on a matplotlib popup window') 
    parser.add_argument('--numex', default=0, type=int, help='Event ID to show from the pickle')   
    args = parser.parse_args()

    prior_matrix = pickle.load(open(args.prior))
    #prior_matrix = prior_matrix.reshape(-1)
    #heatmap, xedges, yedges = histogram2d(prior_matrix, BINS_CENTER, bins=prior_matrix.shape[0])

    imshow(prior_matrix)
    xlabel("PSI_i")
    ylabel("PSI_j")
    legend()
    #extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    #imshow(heatmap, extent=extent)
    _save_or_show(args.plotpath, "prior_posterior")



if __name__ == '__main__':
    main()