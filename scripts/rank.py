"""

Rank MAJIQ, MISO or MATS events to test delta PSI reproducibility


"""
import matplotlib as mplot
mplot.use('Agg')
try:
    import cPickle as pkl
except ImportError:
    import Pickle as pkl

import scripts.utils
import argparse
from pylab import *

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
    MINTHRESHOLD = 0.95
    if E:
        MINTHRESHOLD = 0.20
    rank = []
    
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
            if np.isnan(v_prob): continue
            rank.append([names[i], v_prob])
        else:
            area = matrix_area(dmatrix, V, absolute)
            if np.isnan(area): continue
            if ranknochange:
                #P(Delta PSI < V) = 1 - P(Delta PSI > V)
                area = 1.0 - area
            rank.append([names[i], area])
    if ranknochange:
        rank.sort(key=lambda x: x[1])
    else:
        rank.sort(key=lambda x: x[1], reverse=True)

    num_conf_events = len(rank)
    for idx, v in enumerate(rank):
        if v[1]<MINTHRESHOLD:
            print "FDR=%d" % idx
            num_conf_events = idx
            break
    return rank, num_conf_events


def rank_miso(path, dofilter=True, ranknochange=False, complex_lsvs=False):
    rank = scripts.utils.miso_delta_reader(path, dofilter=dofilter, complex_lsvs=complex_lsvs)
    if ranknochange: 
        rank.sort(key=lambda x: (abs(x[1]), x[2]))
        #sort first by smallest delta PSI, then by bayes factor
    else:
        rank.sort(key=lambda x: (abs(x[1]), x[2]), reverse=True)
        #sort first by biggest delta PSI, then by inverse bayes factor
    return rank


def _is_in_chunk(event1, chunk):
    for event2 in chunk: 
        if event1[0] == event2[0]: #event[0] is the name of the event
            return 1
    return 0

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

def plot_fdr(output, method_name, rank, num_conf_events):

    fig = figure(figsize=[10, 10]) # In inches
    #figure out how many groups of events exist

    font = {'size': 16} #here also 'weight' and 'family'
    matplotlib.rc('font', **font)

    plot(rank, label='FDR %s' % method_name)
    plot((num_conf_events, num_conf_events), (0,1), 'k-')
    legend(loc=2)
    scripts.utils.save_or_show(output, "fdr.%s" % method_name)



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
            majiq_data = pkl.load(open(file, 'r'))
            count_pairs += 1
            rank, num_conf_events = rank_majiq(majiq_data[0], majiq_data[1], V=args.V, absolute=args.absolute, dofilter=args.filter, E=args.E, ranknochange=args.ranknochange, complex_lsvs=args.complex_lsvs)
            ranks.append(rank)

            plot_fdr(args.output, os.path.split(file)[1].split('.')[0], [r[1] for r in ranks[-1][:500]], num_conf_events)

if __name__ == '__main__':
    main()
