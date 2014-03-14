import argparse

from pylab import *
from scipy.stats import norm


def delta_info(newdist, norm_space, delta_limit = 0.2):
    """
    Print out some information on the delta PSI
    """

    for i in xrange(len(norm_space)):
        if norm_space[i] > -delta_limit:
            pos_index_neg = i
            break

    for i in xrange(len(norm_space)-1, 0, -1):
        if norm_space[i] < delta_limit:
            pos_index_pos = i
            break

    #return probability of delta PSI > V AND probability of delta PSI < -V
    return newdist[0:pos_index_neg].sum(), newdist[pos_index_pos+1:].sum()


def _save_or_show(plotpath, name):
    if plotpath:
        savefig("%s_%s.png"%(plotpath, name), width=200, height=400, dpi=100)
        clf()
    else:
        show()  

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--plotpath', default=None, help='')
    args = parser.parse_args()

    binsize = 0.025
    numbins = 80

    for std in arange(0.1, 0.25, 0.05):
        for uniform in arange(0., 6., 1.):
            mydist = norm(loc=0, scale=std)
            norm_space = linspace(-1, 1-binsize, num=numbins*2) + binsize/2
            pdfnorm = mydist.pdf(norm_space)
            label_uniform = uniform
            uniform = uniform/numbins
            newdist = (pdfnorm+uniform)/(pdfnorm+uniform).sum()
            v_pos, v_neg = delta_info(newdist, norm_space)
            title("Center Std. %s, Uniform %s. \nP(Delta PSI) > 0.2: %.3f P(Delta PSI) < -0.2: %.3f"%(std, label_uniform, v_pos, v_neg))
            plot(linspace(-1, 1, num=len(list(newdist))), newdist)
            _save_or_show(args.plotpath, "synth_std%s_uni%s"%(std, label_uniform))

if __name__ == '__main__':
    main()
