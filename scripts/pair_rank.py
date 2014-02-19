"""

Rank MAJIQ, MISO or MATS events to test delta PSI reproducibility


"""
import sys
import argparse
import pickle
from collections import defaultdict 

from pylab import *


def print_matrix(matrix):
    "Print MAJIQ delta PSI matrix in screen"
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            print "%.4f"%matrix[i][j],
        print

    print ret
    print


def collapse_matrix(matrix):
    "Collapse the diagonals probabilities in 1-D and return them"
    collapse = []
    #FOR TEST matrix = array([[0, 1, 2, 3, 4, 500], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], [100, 1, 2, 3, 4, 5], ])

    matrix_corner = matrix.shape[0]+1
    for i in xrange(-matrix_corner, matrix_corner):
        collapse.append(diagonal(matrix, offset=i).sum())   

    return array(collapse)

def mean_matrix(matrix):
    "Get the best point mean"
    #collapse in 1 dimension
    collapse = collapse_matrix(matrix)
    delta_values = linspace(-1, 1, num=len(collapse))
    print delta_values+collapse
    print collapse
    plot(collapse)
    show()
    #UNFINISHED
    
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

def matrix_area(matrix, V=0.2, absolute=True):
    """Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute"""
    collapse = collapse_matrix(matrix)
    #get the delta psi histogram borders based on the size of 'collapse'
    border = _find_delta_border(V, collapse.shape[0])
    #grab the values inside the area of interest
    area = []
    if V < 0: 
        area.append(collapse[0:border+1])
        if absolute: #if absolute V, pick the other side of the array
            area.append(collapse[-border-1:])
    else:
        area.append(collapse[border:])
        if absolute: #if absolute V, pick the other side of the array
            area.append(collapse[0:len(collapse)-border])

    return sum(area)

def v_sum(matrix):
    """
    Calculate sum_v v*P(Delta PSI > V)
    """
    absolute = True
    ret = 0.
    for v in arange(0, 1, 0.1):
        ret += matrix_area(matrix, V=v, absolute=absolute)*v

    return ret

def rank_majiq(path, names, V=0.2, absolute=True, dofilter=True, E=False):
    rank = []
    names = pickle.load(open(names))
    for i, dmatrix in enumerate(pickle.load(open(path))):
        if E:
            rank.append([names[i].split(":")[0], v_sum(dmatrix)])
        else:
            area = matrix_area(dmatrix, V, absolute)
            if area > 0.2 or not dofilter:
                rank.append([names[i].split(":")[0], matrix_area(dmatrix, V, absolute)])

    rank.sort(key=lambda x: -x[1])
    return rank


def miso_reader(path, dofilter=True):
    ret = []
    for line in open(path):
        sline = line.split('\t')
        if sline[0] != "event_name":
            #event_name  sample1_posterior_mean  sample1_ci_low  sample1_ci_high sample2_posterior_mean  
            #sample2_ci_low  sample2_ci_high diff    bayes_factor    isoforms    sample1_counts  
            #sample1_assigned_counts sample2_counts  sample2_assigned_counts chrom   strand  mRNA_starts mRNA_ends
            transcripts = sline[9].split(",")
            delta_psi = sline[7].split(",")
            bayes_factor = sline[8]
            if len(transcripts) == 2: #only interested in 2 transcripts events for now
                if float(bayes_factor) > 2 or not dofilter: # below 2 means no change according to tables
                    exons1 = transcripts[0].split('_')
                    exons2 = transcripts[1].split('_')                
                    trans_name = exons1[0].split('.')[0].strip("'")
                    ret.append([trans_name, float(delta_psi[0]), float(bayes_factor)])
    return ret


def rank_miso(path, dofilter=True):
    rank = miso_reader(path, dofilter)
    rank.sort(key=lambda x: (-abs(x[1]), -x[2])) #sort first by delta PSI, then by bayes factor
    return rank


def rank_mats(path, dofilter=True):
    """
    ID      GeneID      geneSymbol  chr strand  exonStart_0base exonEnd upstreamES  upstreamEE  downstreamES    downstreamEE    ID  IC_SAMPLE_1 SC_SAMPLE_1 IC_SAMPLE_2 SC_SAMPLE_2 IncFormLen  SkipFormLen PValue  FDR IncLevel1   IncLevel2   IncLevelDifference
    10357   "AK076509"  NA  chr16   +   23109119    23109259    23108718    23108851    23109984    23110153    10357   2471,4440,3021,3043,4596    6,35,10,8,11    321,601,628,702,730,572 7,17,17,13,15,22    187 61  0.0 0.0 0.993,0.976,0.99,0.992,0.993    0.937,0.92,0.923,0.946,0.941,0.895  0.062
    10650   "BC013699"  NA  chr11   +   83364794    83364836    83350978    83351171    83365646    83365793    10650   157,260,162,157,196 14,29,18,14,22  6,4,12,
    """
    rank = []
    for line in open(path):
        sline = line.split()
        if sline[0] != "ID":
            geneID = sline[1]
            pvalue = float(sline[-4])
            #fdr = float(sline[-3])
            delta_psi =  float(sline[-1])
            if pvalue < 0.05 or not dofilter:
                rank.append([geneID, delta_psi, pvalue])

    rank.sort(key=lambda x: (-abs(x[1]), x[2]))
    return rank


def _save_or_show(plotpath, name):
    if plotpath:
        figure(figsize=[50, 10])
        savefig("%s_%s.png"%(plotpath, name), width=200, height=200, dpi=100)
        clf()
    else:
        show()  

def _is_in_chunk(event1, chunk):
    for event2 in chunk: 
        if event1[0] == event2[0]: #event[0] is the name of the event
            return 1   

    return 0

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('pair', nargs='+', help='')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--miso', default=False, action = 'store_true', help='Analyze MISO input')
    parser.add_argument('--mats', default=False, action= 'store_true', help='Analyze MATS input')
    parser.add_argument('--max',  default=1000, type = int, help="Max number of events to analyze")
    parser.add_argument('--V', default=0.2, type = float, help="Steps of best events to take")
    parser.add_argument('--E', default=False, action="store_true", help="For MAJIQ, calculate sum_v P(deltaPSI > V)")
    parser.add_argument('--proximity', default=0, type = int, help="How close the 2 events have to be in the ranking in order to ")
    parser.add_argument('--evnames', default=None, nargs="*", help="Event names for both pairs in MAJIQ")
    parser.add_argument('--noabsolute', dest="absolute", default=True, action="store_false", help="Determine if V is absolute")
    parser.add_argument('--nofilter', dest="filter", default=True, action="store_false", help="Skip filtering by BF, p-value, etc")
    parser.add_argument('--fullrank', default=False, action='store_true', help="Benchmark searching for events in full ranking on rank2")
    parser.add_argument('--fdr', default=None, help="In addition to the rank, calculate the False Discovery Rate. Only works with --fullrank")
    
    args = parser.parse_args()

    print "Calculating ranks..."
    if args.miso:
        rank1 = array(rank_miso(args.pair[0], args.filter))
        rank2 = array(rank_miso(args.pair[1], args.filter))
    elif args.mats:
        rank1 = array(rank_mats(args.pair[0], args.filter))
        rank2 = array(rank_mats(args.pair[1], args.filter))
    else:
        rank1 = rank_majiq(args.pair[0], args.evnames[0], args.V, args.absolute, args.filter, args.E)
        rank2 = rank_majiq(args.pair[1], args.evnames[1], args.V, args.absolute, args.filter, args.E)

    print "Calculating the ratios..."
    #calculate the ratios
    ratios = []


    if args.proximity or args.fullrank:
        #Using proximity or full rank window
        if args.proximity: print "Using proximity window of %s..."%args.proximity
        else: print "Using full rank2 for all events %s..."%args.max
        found = 0
        fdr = [0] #zero to avoid using "first" flag for first element
        v_values = []
        for i in xrange(args.max+1):
            if args.proximity:
                min_chunk = max(0, i-args.proximity/2)
                max_chunk = min_chunk+args.proximity

            elif args.fullrank: #check in the whole set instead of subsets
                min_chunk = 0
                max_chunk = args.max

            if i % 20 == 0:
                print "Event rank1 n=%s. Window rank2: %s-%s"%(i, min_chunk, max_chunk)  
            
            #check if event1 is inside the window of rank2
            found += _is_in_chunk(rank1[i], list(rank2[min_chunk:max_chunk]))
            if args.fdr:
                v_values.append(rank1[i][1])
                fdr.append(fdr[-1]+v_values[-1])

            ratios.append(float(found))

        fdr.pop(0) #remove now useless first item
        #normalize ratios
        ratios = array(ratios)
        ratios /= ratios.shape[0]
        if args.fdr: #normalize fdr if we are calculating it
            fdr = array(fdr)
            fdr /= fdr.shape[0]

    else: #"equalrank" chunks of same n size in both ranks
        for i in xrange(args.max+1):
            chunk1 = list(rank1[0:i]) 
            chunk2 = list(rank2[0:i])
            #check if event1 is into chunk2
            found = 0
            for event1 in chunk1:
                found += _is_in_chunk(event1, chunk2)
                        
            ratios.append(float(found) / args.max)
            if i % 20 == 0:
                print "%s..."%i,
                sys.stdout.flush()

        ratios = array(ratios)
    
    print "RESULT:", ratios[0:10], "...", ratios[-10:], "length", ratios.shape
    print "Saving..."
    pickle.dump(ratios, open(args.output, 'w'))

    if args.fdr:
        print "FDR:", fdr[0:10], "...", fdr[-10:], "length", fdr.shape
        pickle.dump(fdr, open("%s.pickle"%args.fdr, 'w'))
        pickle.dump(v_values, open("%s_v.pickle"%args.fdr, 'w'))

    print "Done!"

if __name__ == '__main__':
    main()
