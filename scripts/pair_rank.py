"""

Rank MAJIQ, MISO or MATS events to test delta PSI reproducibility


"""
import argparse
import pickle
import os
from pylab import *


def print_matrix(matrix):
    "Print MAJIQ delta PSI matrix in screen"
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            print "%.4f"%matrix[i][j],
        print

    print
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

def rank_majiq(bins_list, names, V=0.2, absolute=True, dofilter=True, E=False, ranknochange=False):
    MINTHRESHOLD = 0.2 # minimum threshold in order to consider a prob. significant enough to be included in the ranking
    rank = []
    # lsv_types_dict = {
    #     's|1e1.1|1e2.1':'SE',
    #     't|1e1.1|1e2.1':'SE',
    #     # 's|1e1.1|1e1.2':'A3SS',
    #     # 't|1e1.1|2e1.1':'A3SS',
    #     # 't|1e1.1|1e1.2':'A5SS',
    #     # 's|1e1.1|2e1.1':'A5SS'
    # }

    print len(names), len(bins_list)
    for i, lsv_bins in enumerate(bins_list):
        # if names[i][2] not in lsv_types_dict.keys():
        #     continue
        if len(lsv_bins)>2:
            continue
        # for dmatrix in lsv_bins:
        dmatrix = lsv_bins[0]
        if E:
            v_prob = v_sum(dmatrix)
            rank.append([names[i], round(v_prob, 2)])
        else:
            area = matrix_area(dmatrix, V, absolute)
            if ranknochange: #P(Delta PSI < V) = 1 - P(Delta PSI > V)
                area = 1 - area

            if area > MINTHRESHOLD or not dofilter:
                rank.append([names[i], area])
    rank.sort(key=lambda x: x[1], reverse=True)
    return rank


def miso_reader(path, dofilter=True):
    ret = []
    for line in open(path):
        sline = line.split('\t')
        if sline[0] != "event_name":
            #event_name  sample1_posterior_mean  sample1_ci_low  sample1_ci_high sample2_posterior_mean  
            #sample2_ci_low  sample2_ci_high diff    bayes_factor    isoforms    sample1_counts  
            #sample1_assigned_counts sample2_counts  sample2_assigned_counts chrom   strand  mRNA_starts mRNA_ends
            try:
                transcripts = sline[9].split(",")
            except:
                print line
                raise
            delta_psi = sline[7].split(",")
            bayes_factor = sline[8]
            trans_name = sline[0]
            if len(transcripts) == 2: #only interested in 2 transcripts events for now
                ret.append([trans_name, float(delta_psi[0]), float(bayes_factor)])
    return ret


def rank_miso(path, dofilter=True, ranknochange=False):
    rank = miso_reader(path, dofilter)
    if ranknochange: 
        rank.sort(key=lambda x: (abs(x[1]), x[2])) #sort first by smallest delta PSI, then by bayes factor
    else:
        rank.sort(key=lambda x: (abs(x[1]), x[2]), reverse=True) #sort first by biggest delta PSI, then by inverse bayes factor
    return rank


def rank_mats(path, dofilter=True, ranknochange=False):
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

    if ranknochange:
        rank.sort(key=lambda x: (abs(x[1]), x[2])) #biggest delta PSI first, small p-value
    else:
        rank.sort(key=lambda x: (-abs(x[1]), x[2])) #biggest delta PSI first, small p-value

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


def event_names_set_mats(rank):
    pass


def event_names_set_miso(rank):
    return set([ranked_pair[0] for ranked_pair in rank])


def event_names_set_majiq(rank):
    return set([ranked_pair[0][1] for ranked_pair in rank])


def skim_rank(rank, common_names, method):
    # print "Before filter: %d" % len(rank)
    if method == 'majiq':
        names = [ranked_pair[0][1] for ranked_pair in rank]
    if method == 'miso':
        names = [ranked_pair[0] for ranked_pair in rank]

    common_index = np.array([name in common_names for name in names])
    # print "After filter: %d" % np.array(rank)[common_index].size
    return np.array(rank)[common_index]

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--majiq-files', required=True, dest='majiq_files', nargs=2, help='MAJIQ files with paired events to analyze')
    parser.add_argument('--miso-files', dest='miso_files', nargs=2, help='MISO files with paired events to analyze')
    parser.add_argument('--mats-files', dest='mats_files', nargs=2,help='MATS files with paired events to analyze')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--max',  default=1000, type = int, help="Max number of events to analyze")
    parser.add_argument('--V', default=0.2, type = float, help="Steps of best events to take")
    parser.add_argument('--E', default=False, action="store_true", help="For MAJIQ, calculate sum_v P(deltaPSI > V)")
    parser.add_argument('--proximity', default=0, type = int, help="How close the 2 events have to be in the ranking in order to ")
    parser.add_argument('--evnames', default=None, nargs="*", help="Event names for both pairs in MAJIQ")
    parser.add_argument('--noabsolute', dest="absolute", default=True, action="store_false", help="Determine if V is absolute")
    parser.add_argument('--nofilter', dest="filter", default=True, action="store_false", help="Skip filtering by BF, p-value, etc")
    parser.add_argument('--fullrank', default=False, action='store_true', help="Benchmark searching for events in full ranking on rank2")
    parser.add_argument('--fdr', default=None, help="In addition to the rank, calculate the False Discovery Rate. Only works with --fullrank")
    parser.add_argument('--ranknochange', default=False, action='store_true', help="Calculate P(deltaPSI < V) instead of P(deltaPSI > V) to rank first the ones with low delta PSI")
    parser.add_argument('--intersect-events', dest='intersect_events', default=False, action='store_true', help="Intersect the events among all the pairs")
    args = parser.parse_args()

    print args

    print "Calculating ranks..."
    ranks = {'majiq': [],
             'miso': [],
             'mats': []}

    if args.miso_files:
        for file in args.miso_files:
            ranks['miso'].append(array(rank_miso(file, args.filter, args.ranknochange)))
    if args.mats_files:
        for file in args.mats_files:
            ranks['mats'].append(array(rank_mats(file, args.filter, args.ranknochange)))
    if args.majiq_files:
        for file in args.majiq_files:
            majiq_data = pickle.load(open(file, 'r'))
            ranks['majiq'].append(rank_majiq(majiq_data[0], majiq_data[1], args.V, args.absolute, args.filter, args.E, args.ranknochange))

    if args.intersect_events:
        print "Computing intersection of events..."
        common_names = event_names_set_majiq(ranks['majiq'][0]).intersection(event_names_set_majiq(ranks['majiq'][1]))

        print "Initial set (from MAJIQ): %d LSVs" % len(common_names)
        if args.mats_files:
            mats_set = event_names_set_mats(ranks['mats'][0]).intersection(event_names_set_mats(ranks['mats'][1]))
            common_names = common_names.intersection(mats_set)
            print "After intersected with MATS: %d LSVs" % len(common_names)
        else:
            del ranks['mats']
        if args.miso_files:
            miso_set = event_names_set_miso(ranks['miso'][0]).intersection(event_names_set_miso(ranks['miso'][1]))
            common_names = common_names.intersection(miso_set)
            print "After intersected with MISO: %d LSVs" % len(common_names)
        else:
            del ranks['miso']


    for method_name, ranks_pair in ranks.items():
        print "Ranking %s...." % method_name
        if args.intersect_events:
            print "Skimming rankings, discarding events (LSVs) not common for %s..." % method_name
            for ii, rank in enumerate(ranks_pair):
                ranks_pair[ii] = skim_rank(rank, common_names, method_name)

        print "Final lengths: %d, %d"  % (len(ranks_pair[0]), len(ranks_pair[1]))
        rank1, rank2 = ranks_pair
        print "Num events", len(rank1), len(rank2)
        print "Calculating the ratios..."
        #calculate the ratios
        ratios = []

        max_events = min(args.max +1, min(len(rank1), len(rank2)))

        fdr = []
        if args.proximity or args.fullrank:
            #Using proximity or full rank window
            if args.proximity: print "Using proximity window of %s..."%args.proximity
            else: print "Using full rank2 for all events %s..."%args.max
            found = 0
            fdr = [0] #zero to avoid using "first" flag for first element
            v_values = []
            for i in xrange(max_events):
                if args.proximity:
                    min_chunk = max(0, i-args.proximity/2)
                    max_chunk = min_chunk+args.proximity

                elif args.fullrank: #check in the whole set instead of subsets
                    min_chunk = 0
                    max_chunk = max_events

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
            for i in xrange(max_events):
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
        print "Saving... in %s" % args.output

        if not os.path.exists(args.output):
            os.makedirs(args.output)

        pickle.dump(ratios, open(args.output+"/ratios.%s.pickle" % method_name, 'w'))

        if args.fdr:
            #print "FDR:", fdr[0:10], "...", fdr[-10:], "length", fdr.shape
            pickle.dump(fdr, open("%s.%s.pickle" % (args.fdr, method_name), 'w'))
            # pickle.dump(v_values, open("%s_v.pickle"%args.fdr, 'w'))

    print "Done!"

if __name__ == '__main__':
    main()
