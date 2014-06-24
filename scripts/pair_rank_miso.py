"""

Rank MAJIQ, MISO or MATS events to test delta PSI reproducibility


"""
import sys
import argparse
import pickle
from collections import defaultdict 

from pylab import *



def miso_reader(paths, dofilter=True):
    delta_psis = defaultdict(list)
    for i, path in enumerate(paths):
        print path
        for line in open(path):
            sline = line.split('\t')
            if sline[0] != "event_name":
                #event_name  sample1_posterior_mean  sample1_ci_low  sample1_ci_high sample2_posterior_mean sample2_ci_low  sample2_ci_high diff    bayes_factor    isoforms    sample1_counts  sample1_assigned_counts sample2_counts  sample2_assigned_counts chrom   strand  mRNA_starts mRNA_ends
                event_name = sline[0]
                delta_psi = sline[7].split(",")[0]
                bayes_factor = sline[8]
                delta_psis[event_name].append(float(delta_psi)) #deltapsi

    ret = []
    for name, deltapsi in delta_psis.items():
        ret.append([name, random.uniform(-1, 1)])

    return ret


def rank_miso(paths, dofilter=True, ranknochange=False):

    rank = miso_reader(paths, dofilter)

    if ranknochange: 
        rank.sort(key=lambda x: (abs(x[1]))) #sort first by smallest delta PSI, then by bayes factor
    else:
        rank.sort(key=lambda x: (-abs(x[1]))) #sort first by biggest delta PSI, then by biggest bayes factor
    
    #print rank[0:50]
    return rank


def _is_in_chunk(event1, chunk):
    for event2 in chunk: 
        if event1[0] == event2[0]: #event[0] is the name of the event
            return 1   

    return 0

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--pair1', nargs='+', help='')
    parser.add_argument('--pair2', nargs='+', help='')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--max',  default=1000, type = int, help="Max number of events to analyze")
    parser.add_argument('--proximity', default=0, type = int, help="How close the 2 events have to be in the ranking in order to ")
    parser.add_argument('--noabsolute', dest="absolute", default=True, action="store_false", help="Determine if V is absolute")
    parser.add_argument('--nofilter', dest="filter", default=True, action="store_false", help="Skip filtering by BF, p-value, etc")
    parser.add_argument('--fullrank', default=False, action='store_true', help="Benchmark searching for events in full ranking on rank2")
    parser.add_argument('--fdr', default=None, help="In addition to the rank, calculate the False Discovery Rate. Only works with --fullrank")
    parser.add_argument('--ranknochange', default=False, action='store_true', help="Calculate P(deltaPSI < V) instead of P(deltaPSI > V) to rank first the ones with low delta PSI")
    args = parser.parse_args()

    print "Calculating ranks..."
    rank1 = array(rank_miso(args.pair1, args.filter, args.ranknochange))
    rank2 = array(rank_miso(args.pair2, args.filter, args.ranknochange))

    print "Num events", len(rank1), len(rank2)
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
        #print "FDR:", fdr[0:10], "...", fdr[-10:], "length", fdr.shape
        pickle.dump(fdr, open("%s.pickle"%args.fdr, 'w'))
        pickle.dump(v_values, open("%s_v.pickle"%args.fdr, 'w'))


    
    print "Done!"

if __name__ == '__main__':
    main()
