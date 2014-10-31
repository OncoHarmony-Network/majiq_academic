"""

Rank MAJIQ, MISO or MATS events to test delta PSI reproducibility


"""

import matplotlib as mplot
mplot.use('Agg')
from analysis import psi


import scripts.utils
import analysis.psi
# import prettyplotlib as ppl


from collections import defaultdict
import argparse
import pickle
import os
from pylab import *
from pdb import set_trace as st

RANK_TYPES = ['all', 'intersected', 'only_exp1', 'exp1_and_exp2']


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
    try:
        matrix_corner = matrix.shape[0]+1
    except:
        import pdb
        pdb.set_trace()
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

def expected_dpsi(matrix):
    """
    Calculate sum_dpsi=Prob(dpsi)*dpsi == sum_v = v*P(Delta PSI)
    """
    absolute = True
    ret = 0.
    collapsed = collapse_matrix(matrix)
    for i, v in enumerate(linspace(-1, 1, num=collapsed.shape[0])):
        ret += collapsed[i]*abs(v)

    return ret


def rank_majiq(bins_list, names, V=0.2, absolute=True, dofilter=True, E=False, ranknochange=False, complex_lsvs=False, prior=None):
    MINTHRESHOLD = 0. # minimum threshold in order to consider a prob. significant enough to be included in the ranking
    rank = []
    # lsv_types_dict = {
    #     's|1e1.1|1e2.1':'SE',
    #     't|1e1.1|1e2.1':'SE',
    #     # 's|1e1.1|1e1.2':'A3SS',
    #     # 't|1e1.1|2e1.1':'A3SS',
    #     # 't|1e1.1|1e1.2':'A5SS',
    #     # 's|1e1.1|2e1.1':'A5SS'
    # }
    olderr = np.seterr(divide='ignore')

    for i, lsv_bins in enumerate(bins_list):
        # if names[i][2] not in lsv_types_dict.keys():
        #     continue
        if not complex_lsvs and len(lsv_bins) > 2:
            continue
        if ranknochange:
            dmatrix = lsv_bins[0]
            dmatrix /= sum(dmatrix)
        else:
            dmatrix = np.array(lsv_bins[0])

        if E:
            # v_prob = v_sum(dmatrix)
            v_prob = expected_dpsi(dmatrix)
            rank.append([names[i], v_prob])
        else:
            area = matrix_area(dmatrix, V, absolute)
            if ranknochange:
                #P(Delta PSI < V) = 1 - P(Delta PSI > V)
                area = 1.0 - area

            if area > MINTHRESHOLD or not dofilter:
                rank.append([names[i], area])
    if ranknochange:
        rank.sort(key=lambda x: x[1])
    else:
        rank.sort(key=lambda x: x[1], reverse=True)

    for idx, v in enumerate(rank):
        if v[1]<0.95:
            print "FDR=%d" % idx
            break
    return rank


def rank_miso(path, dofilter=True, ranknochange=False, complex_lsvs=False):
    rank = scripts.utils.miso_delta_reader(path, dofilter=dofilter, complex_lsvs=complex_lsvs)
    if ranknochange: 
        rank.sort(key=lambda x: (abs(x[1]), x[2]))
        #sort first by smallest delta PSI, then by bayes factor
    else:
        rank.sort(key=lambda x: (abs(x[1]), x[2]), reverse=True)
        #sort first by biggest delta PSI, then by inverse bayes factor
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
            # if pvalue < 0.05 or not dofilter:
            #     rank.append([geneID, delta_psi, pvalue])
            rank.append([geneID.replace('"', ''), delta_psi, pvalue])

    if ranknochange:
        rank.sort(key=lambda x: (abs(x[1]), x[2])) #biggest delta PSI first, small p-value
    else:
        rank.sort(key=lambda x: (-abs(x[1]), x[2])) #biggest delta PSI first, small p-value

    return rank


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


def create_restrict_plot(ratios_list):
    from scipy.integrate import simps, trapz

    method_name = 'majiq'
    for ratios in ratios_list:
        area_simp = simps(ratios, dx=1)
        area_trap = trapz(ratios, dx=1)
        print "Method: %s. Last point\t\t: %.3f " % (method_name, ratios[-1])
        print "Method: %s. Area Under Curve (Simpson)\t: %.3f " % (method_name, area_simp)
        print "Method: %s. Area Under Curve (Trapezoid)\t: %.3f " % (method_name, area_trap)

def plot_fdr(output, method_name, fdr):

    fig = figure(figsize=[10, 10]) # In inches
    #figure out how many groups of events exist

    font = {'size': 16} #here also 'weight' and 'family'
    matplotlib.rc('font', **font)

    plot(fdr, label='FDR %s' % method_name)
    legend(loc=2)
    scripts.utils._save_or_show(output, "fdr.%s" % method_name)



def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--majiq-files', dest='majiq_files', nargs='*', help='MAJIQ files with paired events to analyze')
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
    parser.add_argument('--nofdr', action="store_false", default=True, help="Do not calculate the False Discovery Rate. Only applicable with --fullrank")
    parser.add_argument('--ranknochange', default=False, action='store_true', help="Calculate P(deltaPSI < V) instead of P(deltaPSI > V) to rank first the ones with low delta PSI")
    parser.add_argument('--type-rank', dest='type_rank', default='all', choices=RANK_TYPES, help='Configure which events are chosen for the ranking.')
    parser.add_argument('--create_restrict_plot', dest='create_restrict_plot', default=False, action='store_true', help="Create plot for only_ex1 ranks in different restrictive conditions. Only works with --type-rank only_exp1")
    parser.add_argument('--complex-lsvs', dest="complex_lsvs", default=False, action="store_true", help="Include complex LSVs")
    args = parser.parse_args()

    ranks = []
    import os

    if args.majiq_files:
        count_pairs = 0
        for file_nr, file in enumerate(args.majiq_files):
            majiq_data = pickle.load(open(file, 'r'))
            count_pairs += 1
            majiq_file1_names = majiq_data[1]
            ranks.append(rank_majiq(majiq_data[0], majiq_data[1], args.V, args.absolute, args.filter, args.E, args.ranknochange, args.complex_lsvs))

        plot_fdr(args.output, os.path.split(file)[1].split('.')[0], [r[1] for r in ranks[-1][:500]])

if __name__ == '__main__':
    main()
