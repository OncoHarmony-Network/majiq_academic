import argparse
import pickle
from collections import defaultdict

from pylab import *

def _find_delta_border(V, numbins):
    "Finds the border index to which a V corresponds in its delta_space given the number of bins the matrix will have"
    delta_space = list(linspace(-1, 1, num=numbins+1))
    delta_space.pop(0) #first border to the left is -1, and we are not interested in it
    #get the index position that corresponds to the V threshold
    for i, value in enumerate(delta_space):
        if value > V:
            return i  
    #if nothing hit, V = 1
    return numbins


def collapse_matrix(matrix):
    "Collapse the diagonals probabilities in 1-D and return them"
    collapse = []
    #FOR TEST matrix = array([[0, 1, 2, 3, 4, 500], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], [100, 1, 2, 3, 4, 5], ])
    matrix_corner = matrix.shape[0]+1
    for i in xrange(-matrix_corner, matrix_corner):
        collapse.append(diagonal(matrix, offset=i).sum())   

    return array(collapse)

def matrix_area(matrix, V=0.2, absolute=True):
    """Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute"""
    collapse = collapse_matrix(matrix)
    #get the delta psi histogram borders based on the size of 'collapse'
    border = _find_delta_border(V, collapse.shape[0])
    #grab the values inside the area of interest
    area = []
    if V < 0: 
        area.append(collapse[0:border+1])
        if absolute: #if absolute V, pick the other side of the array too
            area.append(collapse[-border-1:])
    else:
        area.append(collapse[border:])
        if absolute: #if absolute V, pick the other side of the array too
            area.append(collapse[0:len(collapse)-border])

    return sum(area)

def matrix_e(matrix):
    "Calculates the expected value of delta PSI E()"
    collapse = collapse_matrix(matrix) #one dimensional discrete distribution of psi
    delta_space = list(linspace(-1, 1, num=len(collapse))) #the values that correspond to the different delta psi [-1, 1] 
    e = 0
    for i, value in enumerate(collapse):
        e += value*delta_space[i]

    return e

def main():
    parser = argparse.ArgumentParser(description="rnaseq is a suite of tools for the analysis of Alternative Splicing Events and Alternative Splicing Quantification.")
    parser.add_argument("pcr")
    parser.add_argument("rnaseq")
    parser.add_argument("names")
    parser.add_argument('-V', default=0.1, type=float)
    parser.add_argument('-Y', default=10, type=float, help="RT-PCR scores")
    parser.add_argument('--plotpath')
    args = parser.parse_args()
    
    pcr_results = defaultdict(float)
    print "Loading PCR data..."
    for line in open(args.pcr):
        symbol, name, pcr_score = line.split() 
        pcr_results[name] = float(pcr_score)

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
            positive_score = matrix_area(matrix, V=args.V, absolute=False)
            negative_score = matrix_area(matrix, V=-args.V, absolute=False)
            neutral_score = 1 - positive_score - negative_score
            e_score = -matrix_e(matrix) #minus because its the opposite comparison to MAJIQ (rest-stim VS stim-rest)
            e_scores.append(e_score)
            if abs(pcr_score) > args.Y:
                changing_pcr += 1

            print "%s \nPCR score: %.3f"%(names[i], pcr_results[names[i]])
            positive_condition = (pcr_score > args.Y and positive_score > negative_score and positive_score > neutral_score)
            negative_condition = (pcr_score < -args.Y and negative_score > positive_score and negative_score > neutral_score)
            neutral_condition = ((pcr_score > -args.Y and pcr_score < args.Y) and neutral_score > positive_score and neutral_score > negative_score)

            if positive_condition:
                hit += 1
                changing_hit += 1
                print "P(V > %s): %.3f"%(args.V, positive_score)
                print "Positive hit"
            elif negative_condition:
                hit += 1
                changing_hit += 1
                print "P(V < -%s): %.3f"%(args.V, negative_score)
                print "Negative hit"
            elif neutral_condition:
                hit += 1
                print "P(%s > V > -%s): %.3f"%(args.V, args.V, neutral_score)
                print "Neutral hit"
            else:
                print "P(V > %s): %.3f"%(args.V, positive_score)
                print "P(V < -%s): %.3f"%(args.V, negative_score)
                print "P(%s > V > -%s): %.3f"%(args.V, args.V, neutral_score)
                print "E(P) %.3f"%(e_score)
                print "MISS"

            print 

    #plot configuration
    font = {'size': 16} #here also 'weight' and 'family'
    matplotlib.rc('font', **font)
    figure(figsize=[10, 10])
    xlabel("PCR score")
    ylabel("E(Delta PSI)")
    pcr_lim = 50
    rnaseq_lim = 0.5
    xlim(-pcr_lim, pcr_lim)
    ylim(-rnaseq_lim, rnaseq_lim)
    #regression line
    fit = polyfit(pcr_scores, e_scores, 1)
    fit_fn = poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    plot(pcr_scores, e_scores, 'yo', pcr_scores, fit_fn(pcr_scores), '--k')
    #scatter plot
    #plot(pcr_scores, e_scores, '.')
    plot([-pcr_lim, pcr_lim], [-rnaseq_lim, rnaseq_lim], '--')
    savefig(args.plotpath, width=200, height=200, dpi=100)

    print "RESULTS"
    print "Quantifiable events: %s out of %s"%(found, len(pcr_results))
    print "Hits: %s out of %s (%.3f%%)"%(hit, found, (float(hit)/found)*100)
    print "Changing hits: %s out of %s (%.3f%%)"%(changing_hit, changing_pcr, (float(changing_hit)/changing_pcr)*100)
    print ""


if __name__ == '__main__':
    main()