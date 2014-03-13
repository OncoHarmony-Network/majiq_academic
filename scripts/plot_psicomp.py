import argparse
import pickle

from pylab import *

def _save_or_show(plotpath):
    if plotpath:
        savefig("%s.png"%(plotpath, name), width=200, height=400, dpi=100)
        clf()
    else:
        show()  

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('file1', help='')
    parser.add_argument('file2', help='')
    parser.add_argument('--plotpath', default=None)
    args = parser.parse_args()

    majiqrank = pickle.load(open(args.file1))
    misorank = pickle.load(open(args.file2))

    print len(majiqrank)
    best_miso = 0.
    best_majiq = 0.
    for i, majiq in enumerate(majiqrank):

        if majiq-misorank[i] > 0:
            best_majiq += 1
        else: 
            best_miso += 1

    total = best_miso+best_majiq
    print "%s (%.2f %%) better in MAJIQ (lower triangle). Better in MISO (upper triangle) %s (%.2f %%) total=%s"%(best_majiq, (best_majiq/total)*100, best_miso, (best_miso/total)*100, total) 

    title("L1 distance of PSI between replicas (Liver1 vs Liver4)")
    max_value = max(max(majiqrank), max(misorank))
    xlabel("MISO")
    ylabel("MAJIQ")
    xlim(0, 1)
    ylim(0, 1)
    plot([0, max_value], [0, max_value], linewidth=2)
    font = {'size': 22} #here also 'weight' and 'family'
    
    matplotlib.rc('font', **font)
    plot(majiqrank, misorank, '.')
    _save_or_show(args.plotpath)

if __name__ == '__main__':
    main()
