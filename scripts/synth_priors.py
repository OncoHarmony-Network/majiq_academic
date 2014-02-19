import argparse

from pylab import *
from scipy.stats import norm
from scipy.stats import beta


def _save_or_show(plotpath, name):
    if plotpath:
        savefig("%s_%s.png"%(plotpath, name), width=200, height=400, dpi=100)
        clf()
    else:
        show()  

def main():
    parser = argparse.ArgumentParser(description="Testing synthetic priors")
    parser.add_argument('pickles', nargs='+', help='')
    args = parser.parse_args()
    
    numbins = 80
    binsize = 2./numbins
    mean = 0
    std = 0.1
    delta_limit = 0.2

    norm_space = linspace(-1, 1-binsize, num=numbins) + binsize/2

    mydist = norm(loc=mean, scale=std)
    pdfnorm = mydist.pdf(norm_space)
    pdfnorm /= sum(pdfnorm)

    #some stats about delta PSI
    for i in xrange(len(norm_space)):
        if norm_space[i] > -delta_limit:
            pos_index_neg = i
            break

    for i in xrange(len(norm_space)-1, 0, -1):
        if norm_space[i] < delta_limit:
            pos_index_pos = i
            break

    uniform = 0.2/numbins
    newdist = (pdfnorm+uniform)/(pdfnorm+uniform).sum()
    print "Delta PSI > %s (%.4f)"%(delta_limit, newdist[0:pos_index_neg].sum())
    print "Delta PSI < -%s (%.4f)"%(delta_limit, newdist[pos_index_pos+1:].sum())
    title("norm(mean=%s, std=%s) + uniform Delta PSI > 0.2 = %.2f Delta PSI < -0.2 = %.2f"%(mean, std, newdist[pos_index_pos+1:].sum(), newdist[0:pos_index_neg].sum()), fontsize=20)
    xlabel("Delta PSI", fontsize=15)
    ylabel("Density", fontsize=15)
    plot(norm_space, newdist)
    show()


if __name__ == '__main__':
    main()
