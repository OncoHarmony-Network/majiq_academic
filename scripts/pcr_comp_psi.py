import argparse
import pickle
from collections import defaultdict

from pylab import *

def psi_e(psi_dist):
    "Calculates the expected value of PSI"
    delta_space = list(linspace(0, 1, num=len(psi_dist))) #the values that correspond to the different delta psi [-1, 1] 
    e = 0
    for i, value in enumerate(psi_dist):
        e += value*delta_space[i]

    return e

def main():
    parser = argparse.ArgumentParser(description="rnaseq is a suite of tools for the analysis of Alternative Splicing Events and Alternative Splicing Quantification.")
    parser.add_argument("pcr")
    parser.add_argument("rnaseq")
    parser.add_argument("names")
    parser.add_argument('-V', default=0.1, type=float)
    parser.add_argument('-Y', default=10, type=float, help="RT-PCR scores")
    parser.add_argument('--stim', action='store_true', default=False)
    parser.add_argument('--plotpath')
    args = parser.parse_args()
    
    pcr_results = defaultdict(float)
    print "Loading PCR data..."
    for line in open(args.pcr):
        if line:
            sline = line.split()
            symbol = sline[0]
            name = sline[1]
            if args.stim:
                psi_pcrs = sline[2:6]
            else:
                psi_pcrs = sline[6:]

            found = 0
            for psi_pcr in psi_pcrs:
                try:
                    psi_pcr = float(psi_pcr)
                    pcr_results[name] += psi_pcr
                    found += 1
                except ValueError:
                    pass

            if found: 
                pcr_results[name] /= found


    #print pcr_results
    print "Loading names..."
    names = [line for line in pickle.load(open(args.names))]

    print "Checking rnaseq..."
    found = 0
    hit = 0
    changing_hit = 0
    changing_pcr = 0
    e_scores = []
    pcr_scores = []
    for i, matrix in enumerate(pickle.load(open(args.rnaseq))):
        if names[i] in pcr_results: 
            found += 1 
            pcr_score = pcr_results[names[i]]
            pcr_scores.append(pcr_score)
            e_score = psi_e(matrix) #minus because its the opposite comparison to MAJIQ (rest-stim VS stim-rest)
            e_scores.append(e_score)

    print "FOUND", found
    #plot configuration
    font = {'size': 16} #here also 'weight' and 'family'
    matplotlib.rc('font', **font)
    figure(figsize=[10, 10])
    xlabel("PCR score")
    ylabel("E(PSI)")
    pcr_lim = 100
    rnaseq_lim = 1
    xlim(0, pcr_lim)
    ylim(0, rnaseq_lim)
    #regression line
    fit = polyfit(pcr_scores, e_scores, 1)
    fit_fn = poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    plot(pcr_scores, e_scores, 'yo', pcr_scores, fit_fn(pcr_scores), '--k')
    #scatter plot
    #plot(pcr_scores, e_scores, '.')
    plot([-pcr_lim, pcr_lim], [-rnaseq_lim, rnaseq_lim], '--')
    savefig(args.plotpath, width=200, height=200, dpi=100)




if __name__ == '__main__':
    main()