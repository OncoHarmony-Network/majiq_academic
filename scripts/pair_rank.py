"""

Rank MAJIQ, MISO or MATS events to test delta PSI reproducibility


"""
import matplotlib as mplot
mplot.use('Agg')
from voila.vlsv import collapse_matrix

import scripts.utils
# import prettyplotlib as ppl


from collections import defaultdict
import argparse
import os
from pylab import *
try:
    import cPickle as pickle
except:
    import pickle

RANK_TYPES = ['all', 'intersected', 'only_exp1', 'exp1_and_exp2']


def _find_delta_border(V, numbins):
    """Finds the border index to which a V corresponds in its delta_space given the number of bins the matrix will have"""
    delta_space = list(linspace(-1, 1, num=numbins + 1))
    delta_space.pop(0)  # first border to the left is -1, and we are not interested in it
    # get the index position that corresponds to the V threshold
    for i, value in enumerate(delta_space):
        if value > V:
            return i
    # if nothing hit, V = 1
    return numbins


def matrix_area(matrix, V=0.2, absolute=True, collapsed_mat=False):
    """Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute"""
    collapse = matrix
    if not collapsed_mat:
        collapse = collapse_matrix(matrix)
    # get the delta psi histogram borders based on the size of 'collapse'
    border = _find_delta_border(V, collapse.shape[0])
    # grab the values inside the area of interest
    area = []
    if V < 0:
        area.append(collapse[0:border + 1])
        if absolute:  # if absolute V, pick the other side of the array
            area.append(collapse[-border - 1:])
    else:
        area.append(collapse[border:])
        if absolute:  # if absolute V, pick the other side of the array
            area.append(collapse[0:len(collapse) - border])

    return sum(area)


def v_sum(matrix):
    """
    Calculate sum_v v*P(Delta PSI > V)
    """
    absolute = True
    ret = 0.
    for v in arange(0, 1, 0.1):
        ret += matrix_area(matrix, V=v, absolute=absolute) * v

    return ret


def expected_dpsi(matrix):
    """
    Calculate sum_dpsi=Prob(dpsi)*dpsi == sum_v = v*P(Delta PSI)
    """
    absolute = True
    ret = 0.
    collapsed = collapse_matrix(matrix)
    for i, v in enumerate(linspace(-1, 1, num=collapsed.shape[0])):
        ret += collapsed[i] * abs(v)

    return ret


def rank_majiq(vlsv_list, V=0.2, absolute=True, dofilter=True, E=False, ranknochange=False, complex_lsvs=False,
               prior=None, shrink=True, junc_selection=None, majiq_n=None):
    MINTHRESHOLD = 0.95
    if E:
        MINTHRESHOLD = 0.20
    rank = []

    print "Num of LSVs in majiq: %d" % len(vlsv_list)
    for i, vlsv in enumerate(vlsv_list):
        lsv_bins = vlsv.get_bins()
        junc_n = 0
        if 'i' in vlsv.get_type(): continue
        if len(lsv_bins) > 2:
            if junc_selection:
                bins_selected = lsv_bins[junc_selection[vlsv.get_id()]]
                junc_n = junc_selection[vlsv.get_id()]
            else:
                most_change = 0
                for jj, junc_bins in enumerate(lsv_bins):
                    if abs(expected_dpsi(np.array(junc_bins))) > most_change:
                        bins_selected = junc_bins
                        junc_n = jj
        else:
            bins_selected = lsv_bins[0]

        if ranknochange:
            dmatrix = np.exp(np.log(bins_selected))
            dmatrix /= sum(dmatrix)
        else:
            dmatrix = np.array(bins_selected)

        v_expected = expected_dpsi(dmatrix)
        area = 1.0 - matrix_area(dmatrix, V, absolute)  # P(Delta PSI < V) = 1 - P(Delta PSI > V)
        rank.append(["%s#%d" % (vlsv.get_id(), junc_n), v_expected, area, int(abs(v_expected)>=V and area<=0.05)])

    expected_mask = np.array([abs(r[1]) >= V for r in rank])
    fdr_mask = np.array([r[2] <= 0.05 for r in rank])
    expected_fdr_mask = np.logical_and(expected_mask, fdr_mask)
    # rank = np.array(rank)[expected_fdr_mask].tolist()
    if majiq_n:
        majiq_n[0] = np.count_nonzero(expected_fdr_mask)
    # rank = np.array(rank)[np.logical_and(expected_mask, fdr_mask)].tolist() #TODO: Remove this!!

    print "#FDR < 0.05: %d" % np.count_nonzero(fdr_mask)
    print "#E(Delta(PSI))>%.2f: %d" % (V, np.count_nonzero(expected_mask))
    print "#E(Delta(PSI))>%.2f and FDR<0.05: %d" % (V, np.count_nonzero(expected_fdr_mask))

    rank.sort(key=lambda x: (x[1], x[2]), reverse=True)

    return np.array(rank)


def rank_naive(bins_list, names, V=0.2, absolute=True, E=False, ranknochange=False, complex_lsvs=False):
    """Similar to MAJIQ files with the difference that the matrix is already collapsed"""
    rank = []

    print "Num of LSVs in naive_bootstrapping: %d" % len(bins_list)
    for i, lsv_bins in enumerate(bins_list):

        junc_n = -1
        most_change = 0
        dmatrix = lsv_bins

        # Expected dpsi
        v_prob = 0.
        for ii, v in enumerate(linspace(-1, 1, num=dmatrix.shape[0])):
            v_prob += dmatrix[ii] * abs(v)

        area = 1. - matrix_area(dmatrix, V, absolute, collapsed_mat=True)
        rank.append([names[i][1], v_prob, area, 1])

    expected_mask = np.array([abs(r[1]) >= V for r in rank])
    fdr_mask = np.array([r[2] <= 0.05 for r in rank])
    expected_fdr_mask = np.logical_and(expected_mask, fdr_mask)

    print "#FDR < 0.05: %d" % np.count_nonzero(fdr_mask)
    print "#E(Delta(PSI))>%.2f: %d" % (V, np.count_nonzero(expected_mask))
    print "#E(Delta(PSI))>%.2f and FDR<0.05: %d" % (V, np.count_nonzero(expected_fdr_mask))

    rank.sort(key=lambda x: (x[1], x[2]), reverse=True)
    return rank


def rank_miso(path, dofilter=True, ranknochange=False, complex_lsvs=False):
    rank = scripts.utils.miso_delta_reader(path, dofilter=dofilter, complex_lsvs=complex_lsvs)
    if ranknochange:
        rank.sort(key=lambda x: (abs(x[1]), x[2]))
        # sort first by smallest delta PSI, then by bayes factor
    else:
        rank.sort(key=lambda x: (abs(x[1]), x[2]), reverse=True)
        # sort first by biggest delta PSI, then by inverse bayes factor
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
            # pvalue = float(sline[-4])
            fdr = float(sline[-4])
            delta_psi = float(sline[-1])
            # if pvalue < 0.05 or not dofilter:
            # rank.append([geneID, delta_psi, pvalue])
            rank.append([geneID.replace('"', ''), delta_psi, fdr])

    if ranknochange:
        rank.sort(key=lambda x: (abs(x[1]), x[2]))  # biggest delta PSI first, small p-value
    else:
        rank.sort(key=lambda x: (-abs(x[1]), x[2]))  # biggest delta PSI first, small p-value

    return rank


def rank_mats_original(mats_file, dofilter=True, ranknochange=False, majiq_n=None):
    """Rank Splicing Events detected as differentially expressed by MATS. Compute the FDR for downstream analysis

    :param mats_file:
    :param dofilter:
    :param ranknochange:
    :return: list of events. Events are represented as a list where the 1st pos is the ID and the 2nd the deltapsi.
             NOTE: Since MATS generates IDs that are not unique among runs, we create it as the concatenation of all
             exonic coordinates (3 exons, except special case of MXE with 4).
    """

    # ID	GeneID	geneSymbol	chr	strand	longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE	ID	IC_SAMPLE_1	SC_SAMPLE_1	IC_SAMPLE_2	SC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
    # 1292	"ENSMUSG00000032366"	"Tpm1"	9	-	67022592	67024565	67022592	67023441	67031028	67031098	1292	1484,2687,1928	7,12,10	17,11,18	24,23,39	1171	61	0.0	0.0	0.917,0.921,0.909	0.036,0.024,0.023	0.888
    # 194	"ENSMUSG00000025199"	"Chuk"	19	-	44078888	44079266	44078888	44078985	44081945	44081995	194	59,87,60	3,8,2	34,51,61	57,70,38	310	43	0.0	0.0	0.732,0.601,0.806	0.076,0.092,0.182	0.596

    first_exon_field = 5
    last_exon_field = 11
    rank = []
    mats_nn = 0
    for line in open(mats_file):
        sline = line.split()
        if sline[0] == "ID":
            if sline[5] == "riExonStart_0base":
                last_exon_field = 13
            else:
                is_mxi = last_exon_field = 11
        else:
            geneID = "".join(sline[first_exon_field - 2:last_exon_field])
            pvalue = float(sline[-5])
            fdr = float(sline[-4])
            delta_psi = float(sline[-1])
            pass_thres =int(abs(delta_psi)>=0.2 and fdr<=0.05) 
            mats_nn += pass_thres 
            rank.append([geneID, delta_psi, pvalue, fdr, pass_thres])

    expected_mask = np.array([abs(r[1]) >= 0.2 for r in rank])
    fdr_cutoff = 0.05
    if majiq_n[0]:
        while np.count_nonzero(np.logical_and(expected_mask, np.array([r[3] <= fdr_cutoff for r in rank]))) < majiq_n[
            0]:
            fdr_cutoff += 0.05
    else:
        majiq_n[0] = np.count_nonzero(np.array([r[3] <= 0.05 for r in rank]))
    fdr_mask = np.array([r[3] <= fdr_cutoff for r in rank])

    print "MATS:"
    print "#FDR < 0.05: %d" % np.count_nonzero(np.array([r[3] <= 0.05 for r in rank]))
    #rank = np.array(rank)[np.logical_and(expected_mask, fdr_mask)].tolist() # TODO: Remove this!!
    print "#FDR < %.2f: %d" % (fdr_cutoff, np.count_nonzero(fdr_mask))
    print "#E(Delta(PSI))>0.20: %d" % np.count_nonzero(expected_mask)
    print "#E(Delta(PSI))>0.20 and FDR<0.05: %d" % np.count_nonzero(np.logical_and(expected_mask, fdr_mask))

    if ranknochange:
        rank.sort(key=lambda x: (abs(float(x[1])), float(x[2])))  # biggest delta PSI first, small p-value
    else:
        rank.sort(key=lambda x: (-abs(float(x[1])), float(x[2])))  # biggest delta PSI first, small p-value

    return mats_nn, rank


def _is_in_chunk(event1, chunk, report_rank2_expec=False):
    for event2 in chunk:
        if event1[0] == event2[0]:  # event[0] is the name of the event
            if report_rank2_expec:
                return [1, event2[1]]
            return 1
    if report_rank2_expec:
        return [0, '']
    return 0


def event_names_set_mats(rank):
    pass


def event_names_set_miso(rank):
    return set([ranked_pair[0] for ranked_pair in rank])


def event_names_set_majiq(rank):
    return set([ranked_pair[0][1] for ranked_pair in rank])


def event_names_set(rank):
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


# python ~/Projects/majiq/scripts/pair_rank.py --majiq-files
# ~/workspace/majiq/data/deltapsi/genomewise/Hippo1_Liver1_deltamatrix.pickle
# ~/workspace/majiq/data/deltapsi/genomewise/results_nofilter/Hippo1_Liver1_deltamatrix.pickle
# ~/workspace/majiq/data/deltapsi/genomewise/Hippo2_Liver2_deltamatrix.pickle
# ~/workspace/majiq/data/deltapsi/genomewise/results_nofilter/Hippo2_Liver2_deltamatrix.pickle
# --output output/repro/pair_rank_all/nooutlier/exp1_union_exp2 --nofilter --type-rank exp1_and_exp2 --max 300 $conf_ranks

# python ~/Projects/majiq/scripts/pair_rank.py
# --majiq-files /data/MGP/majiq/dpsi/Hippo_Liver/pos5reads10/Hippo1_Liver1.deltapsi.pickle /data/MGP/majiq/dpsi/Hippo_Liver/pos2reads2/Hippo5_Liver5.deltapsi.pickle
# --miso-files miso/Hip1Liv1/Hip1_vs_Liv1/bayes-factors/Hip1_vs_Liv1.miso_bf miso/Hip5Liv5/Hip5_vs_Liv5/bayes-factors/Hip5_vs_Liv5.miso_bf
# --mats-files mats/Hip1_Liv1/MATS_output/all.txt mats/Hip5_Liv5/MATS_output/all.txt
# --output repro/pair_rank/H1L1H5L5/full.E/ --fullrank --E --type-rank only_exp1 --max 154


def plot_fdr(output, method_name, fdr):
    diagonaly = np.linspace(0, 1, len(fdr))
    diagonalx = np.linspace(0, len(fdr), len(fdr))

    fig = figure(figsize=[10, 10])  # In inches
    # figure out how many groups of events exist

    font = {'size': 16}  # here also 'weight' and 'family'
    matplotlib.rc('font', **font)

    plot(diagonalx, diagonaly, '--', color="#cccccc")
    plot(fdr, label='FDR %s' % method_name)
    legend(loc=2)
    scripts.utils.save_or_show(output, "fdr.%s" % method_name)


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--majiq-files', dest='majiq_files', nargs='+',
                        help='MAJIQ files with paired events to analyze')
    parser.add_argument('--miso-files', dest='miso_files', nargs=2, help='MISO files with paired events to analyze')
    parser.add_argument('--mats-files', dest='mats_files', nargs=2, help='MATS files with paired events to analyze')
    parser.add_argument('--naive-files', dest='naive_files', nargs=2,
                        help='Naive Bootstrapping files with paired events to analyze')
    parser.add_argument('-o', '--output', required=True, help='Output file_str path')
    parser.add_argument('--max', default=1000, type=int, help="Max number of events to analyze")
    parser.add_argument('--V', default=0.2, type=float, help="Steps of best events to take")
    parser.add_argument('--E', default=False, action="store_true", help="For MAJIQ, calculate sum_v P(deltaPSI > V)")
    parser.add_argument('--proximity', default=0, type=int,
                        help="How close the 2 events have to be in the ranking in order to ")
    parser.add_argument('--evnames', default=None, nargs="*", help="Event names for both pairs in MAJIQ")

    parser.add_argument('--noabsolute', dest="absolute", default=True, action="store_false",
                        help="Determine if V is absolute")
    parser.add_argument('--nofilter', dest="filter", default=True, action="store_false",
                        help="Skip filtering by BF, p-value, etc")
    parser.add_argument('--fullrank', default=False, action='store_true',
                        help="Benchmark searching for events in full ranking on rank2")
    parser.add_argument('--fdr', action="store_true", default=None,
                        help="In addition to the rank, calculate the False Discovery Rate. Only works with --fullrank")
    parser.add_argument('--ranknochange', default=False, action='store_true',
                        help="Calculate P(deltaPSI < V) instead of P(deltaPSI > V) to rank first the ones with low delta PSI")
    parser.add_argument('--intersect-events', dest='intersect_events', default=False, action='store_true',
                        help="Intersect the events among all the pairs")
    parser.add_argument('--type-rank', dest='type_rank', default=RANK_TYPES[2], choices=RANK_TYPES,
                        help='Configure which events are chosen for the ranking.')
    parser.add_argument('--create_restrict_plot', dest='create_restrict_plot', default=False, action='store_true',
                        help="Create plot for only_ex1 ranks in different restrictive conditions. Only works with --type-rank only_exp1")
    parser.add_argument('--complex-lsvs', dest="complex_lsvs", default=False, action="store_true",
                        help="Include complex LSVs")
    parser.add_argument('--noshrink', dest='shrink', default=True, action='store_false',
                        help="Shrink ranks with the FDR number.")
    parser.add_argument('--mats_n', default=False, action='store_true',
                        help="Use MATS number of confident changing events (N of FDR<0.05, |E(Delta(PSI))|>.2)")
    parser.add_argument('--events-only', dest='only_events', default=False, action='store_true',
                        help="Create files with ONLY the events files")
    args = parser.parse_args()
    args = parser.parse_args()

    print args

    print "Calculating ranks..."
    ranks = defaultdict(list)
    n1 = defaultdict(list)

    majiq_N = [None]
    if args.majiq_files:

        count_pairs = 0
        for file_nr, file_str in enumerate(args.majiq_files):
            majiq_data = pickle.load(open(file_str, 'r'))
            # prior = pickle.load(open(str(file_str).replace('deltamatrix', 'priormatrix_jun_0')))
            if file_nr % 2 == 0:
                count_pairs += 1
                # majiq_file1_names = [vlsv.get_id() for vlsv in majiq_data.get_lsvs()]
                ranks['majiq_' + str(count_pairs)].append(
                    rank_majiq(majiq_data.get_lsvs(), args.V, args.absolute, args.filter, args.E, args.ranknochange,
                               args.complex_lsvs, shrink=args.shrink, majiq_n=majiq_N))
                majiq_file1_names = [a[0].split('#')[0] for a in ranks['majiq_1'][0]]
                continue

            if args.type_rank == 'only_exp1':
                # Select events from experiment 1
                exp1_index = np.array([v_lsv.get_id() in majiq_file1_names for v_lsv in majiq_data.get_lsvs()])
                junc_dict = dict(
                    [(rr[0].split('#')[0], int(rr[0].split('#')[1])) for rr in ranks['majiq_' + str(count_pairs)][-1]])
                ranks['majiq_' + str(count_pairs)].append(
                    rank_majiq(np.array(majiq_data.get_lsvs())[exp1_index].tolist(), args.V, args.absolute, args.filter,
                               args.E, args.ranknochange, args.complex_lsvs, shrink=args.shrink,
                               junc_selection=junc_dict))
                n1['majiq_' + str(count_pairs)] = [np.count_nonzero(exp1_index), np.count_nonzero(exp1_index)]
        names_majiq_exp1 = [m[1] for m in majiq_file1_names]
    
    mats_n_orig = None

    if args.mats_files:
        if args.mats_n:
            majiq_N = [None]
        for file_str in args.mats_files:
            mmats_n, mats_rank = rank_mats_original(file_str, args.filter, args.ranknochange, majiq_n=majiq_N)
            if not mats_n_orig:
                mats_n_orig = mmats_n
            if len(ranks['mats']) > 0:
                names_mats_exp1 = [mats_info[0] for mats_info in ranks['mats'][0]]
                names_mats = [mats_info[0] for mats_info in mats_rank]
                exp1_index = np.array([name in names_mats_exp1 for name in names_mats])
                ranks['mats'].append(np.array(mats_rank)[exp1_index])
                n1['mats'].append(np.count_nonzero(exp1_index))
            else:
                ranks['mats'].append(array(mats_rank))

    if args.miso_files:
        for file_nr, file_str in enumerate(args.miso_files):
            miso_rank = rank_miso(file_str, args.filter, args.ranknochange, args.complex_lsvs)
            if args.type_rank == 'only_exp1':
                if file_nr % 2:
                    names_miso = [miso_info[0] for miso_info in miso_rank]
                    names_miso_exp1 = [miso_info[0] for miso_info in ranks['miso'][0]]
                    exp1_index = np.array([name in names_miso_exp1 for name in names_miso])
                    ranks['miso'].append(array(miso_rank)[exp1_index])
                    n1['miso'].append(np.count_nonzero(exp1_index))
                else:
                    ranks['miso'].append(array(miso_rank))

    if args.naive_files:
        names_naive_exp1 = None
        for file_str in args.naive_files:
            naive_data = pickle.load(open(file_str, 'r'))
            naive_rank = rank_naive(naive_data[1], naive_data[0], args.V, args.absolute, args.E, args.ranknochange)
            if args.type_rank == 'only_exp1':
                # Use only MAJIQ selected events for experiment 1
                names_naive = [naive_info[0] for naive_info in naive_rank]
                if not names_naive_exp1:
                    names_naive_exp1 = names_naive
                exp1_index = np.array([name in names_naive_exp1 for name in names_naive])
                ranks['naive'].append(array(naive_rank)[exp1_index])
                n1['naive'].append(np.count_nonzero(exp1_index))
            else:
                ranks['naive'].append(array(naive_rank))

    only_exp1_ranks = []
    for method_name, ranks_pair in ranks.items():
        print "Ranking %s...." % method_name
        rank1, rank2 = ranks_pair
        print "Num events", len(rank1), len(rank2)
        print "Calculating the ratios..."
        # calculate the ratios
        ratios = []
        events = []

        max_events = min(args.max, len(rank1))
        if majiq_N[0] and not args.only_events:
            max_events = min(max_events, majiq_N[0])

        fdr = []
        if args.proximity or args.fullrank:
            # Using proximity or full rank window
            if args.proximity:
                print "Using proximity window of %s..." % args.proximity
            else:
                print "Using full rank2 for all events %s..." % max_events
            found = 0
            fdr = [0]  #zero to avoid using "first" flag for first element
            v_values = []
            for i in xrange(max_events):
                if args.proximity:
                    min_chunk = max(0, i - args.proximity / 2)
                    max_chunk = min_chunk + args.proximity

                elif args.fullrank:  #check in the whole set instead of subsets
                    min_chunk = 0
                    max_chunk = max_events

                if i % 20 == 0:
                    print "Event rank1 n=%s. Window rank2: %s-%s" % (i, min_chunk, max_chunk)

                #check if event1 is inside the window of rank2
                is_hit, rank2_exp = _is_in_chunk(rank1[i], list(rank2[min_chunk:max_chunk]), report_rank2_expec=True)
                found += is_hit
                if args.fdr:
                    v_values.append(rank1[i][1])
                    fdr.append(fdr[-1] + v_values[-1])

                ratios.append(float(found))
                events.append([rank1[i], is_hit, rank2_exp, majiq_N[0] if 'mats' not in method_name else mats_n_orig])

            fdr.pop(0)  #remove now useless first item
            #normalize ratios
            ratios = array(ratios)
            ratios /= ratios.shape[0]
            if args.fdr:  #normalize fdr if we are calculating it
                fdr = array(fdr)
                fdr /= fdr.shape[0]

        else:  # "equalrank" chunks of same n size in both ranks
            import sys

            for i in xrange(max_events):
                chunk1 = list(rank1[0:i + 1])
                chunk2 = list(rank2[0:i + 1])
                # check if event1 is into chunk2
                found = 0
                for event1 in chunk1:
                    found += _is_in_chunk(event1, chunk2)
                    if i == max_events - 1:
                        events.append([event1, _is_in_chunk(event1, chunk2)])
                ratios.append(float(found) / max_events)
                if i % 20 == 0:
                    print "%s..." % i,
                    sys.stdout.flush()

            ratios = array(ratios)

        print "RESULT:", ratios[0:10], "...", ratios[-10:], "length", ratios.shape
        print "Saving in %s" % args.output

        if not os.path.exists(args.output):
            os.makedirs(args.output)

        pickle.dump(ratios,
                    open(args.output + "/ratios.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name),
                         'w'))

        print "Saving events... in %s " % args.output
        pickle.dump(events,
                    open(args.output + "/events.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name),
                         'w'))

        print "Saving N1 size... in %s " % args.output
        pickle.dump(n1[method_name],
                    open(args.output + "/n1.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name), 'w'))

        if args.fdr:
            # print "FDR:", fdr[0:10], "...", fdr[-10:], "length", fdr.shape
            pickle.dump(fdr,
                        open("%s/fdr.%s.%s.pickle" % (args.output, method_name, str(args.type_rank).replace('-', '_')),
                             'w'))
            pickle.dump(v_values, open(
                "%s/fdr.%s.%s_v.pickle" % (args.output, method_name, str(args.type_rank).replace('-', '_')), 'w'))
            plot_fdr(args.output, method_name, fdr)

        if "majiq" in method_name:
            only_exp1_ranks.append(ratios)

    if args.create_restrict_plot:
        create_restrict_plot(ranks)

    print "Done!"


if __name__ == '__main__':
    main()
