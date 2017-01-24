"""
Rank MAJIQ, MISO or MATS events to test delta PSI reproducibility
"""

import matplotlib as mplot
mplot.use('Agg')
import scripts.utils
from collections import defaultdict
import argparse
from pylab import *
try:
    A = np.array([])
    del A
except:
    import numpy as np
    from numpy import *
try:
    import cPickle as pickle
except:
    import pickle
import h5py

import voila
print voila.__path__

from voila.io_voila import VoilaInput as VI
from voila.vlsv import VoilaLsv
from voila.splice_graphics import JunctionGraphic, LsvGraphic, ExonGraphic
import traceback


class VoilaInput(VI):
    def __init__(self, bins_hdf5, graphs):
        self.lsvs = []
        self.metainfo = None
        psi_grps = []
        last_tb = ""
        keys = sorted(bins_hdf5.keys())
        for key in keys:
            if key == 'metainfo':
                continue
            elif type(bins_hdf5[key]) is h5py.Group:
                lsv_grp = bins_hdf5[key]
            elif len(bins_hdf5[key].shape) == 3:
                dpsi_grp = (key, bins_hdf5[key])
            else:
                psi_grps.append((key, bins_hdf5[key]))
        for lsv_id, lsv_attrs in lsv_grp.items():
            try:
                bins_list = self.__get_bins_list(lsv_attrs, *dpsi_grp)
                psi1 = self.__get_bins_list(lsv_attrs, *psi_grps[0])
                psi2 = self.__get_bins_list(lsv_attrs, *psi_grps[1])
                cur_gfx = graphs[lsv_id]
            except:
                last_tb = traceback.format_exc(0).split('\n')[1]
                print 'ERROR: Unable to load %s (%s)' % (lsv_id, last_tb)
            else:
                self.lsvs.append(VoilaLsv(bins_list, cur_gfx, psi1 = psi1, psi2 = psi2))

    @staticmethod
    def __get_bins_list(lsv_attrs, kk, gg):
        return [junc for junc in gg[lsv_attrs.attrs[kk]]]


def get_ratio_nu(rank1, rank2, N, max_N=2000):
    LSVs_1 = map(lambda xx: xx[0], rank1)
    LSVs_2 = map(lambda xx: xx[0], rank2)
    frac_max_N = 1.0/max_N
    return N, np.array([intersect1d(LSVs_1[:n+1], LSVs_2[:n+1]).size *
                        frac_max_N for n in range(max_N)])


def legacy_lsv_junctions_new_voila(graphs_hdf5):
    output = {}
    for lsv_id, lsv_attrs in graphs_hdf5['lsvs'].items():
        visual = lsv_attrs
        exons = []
        for e_id, exon in visual['exons'].items():
            E_id = int(e_id)
            while E_id >= len(exons):
                exons.append(None)
            kwargs = dict(exon.attrs)
            # if 'start' in kwargs:
            #     kwargs['coords'] = [kwargs.pop('start'), kwargs.pop('end')]
            if 'type_exon' not in kwargs:
                kwargs['exon_type_list'] = [0]
            exons[E_id] = ExonGraphic(**kwargs)
        junctions = []
        for j_id, junction in visual['junctions'].items():
            J_id = int(j_id)
            while J_id >= len(junctions):
                junctions.append(None)
            kwargs = dict(junction.attrs)
            kwargs['nreads'] = kwargs.pop('num_reads', 0)
            kwargs['clean_nreads'] = kwargs.pop('num_clean_reads', 0)
            if 'start' in kwargs:
                kwargs['coords'] = [kwargs.pop('start'), kwargs.pop('end')]
            if 'intron_retention' in kwargs:
                kwargs['ir'] = kwargs.pop('intron_retention')
            if 'type_junction' not in kwargs:
                kwargs['type_junction'] = 0
            junctions[J_id] = JunctionGraphic(**kwargs)
        lsv_type = lsv_attrs.attrs.get('lsv_type')
        coords = [lsv_attrs.attrs.get('start'), lsv_attrs.attrs.get('end')]
        output[lsv_id] = LsvGraphic(lsv_type, coords, id = lsv_id, exons = exons, junctions = junctions)
    return output


def legacy_lsv_junctions(graphs_hdf5):
    output = {}
    for lsv_id, lsv_attrs in graphs_hdf5['LSVs'].items():
        visual = lsv_attrs['visual']
        exons = []
        for e_id, exon in visual['exons'].items():
            E_id = int(e_id)
            while E_id >= len(exons):
                exons.append(None)
            kwargs = dict(exon.attrs)
            if 'start' in kwargs:
                kwargs['coords'] = [kwargs.pop('start'), kwargs.pop('end')]
            if 'type_exon' not in kwargs:
                kwargs['type_exon'] = 0
            exons[E_id] = ExonGraphic(**kwargs)
        junctions = []
        for j_id, junction in visual['junctions'].items():
            J_id = int(j_id)
            while J_id >= len(junctions):
                junctions.append(None)
            kwargs = dict(junction.attrs)
            kwargs['nreads'] = kwargs.pop('num_reads', 0)
            kwargs['clean_nreads'] = kwargs.pop('num_clean_reads', 0)
            if 'start' in kwargs:
                kwargs['coords'] = [kwargs.pop('start'), kwargs.pop('end')]
            if 'intron_retention' in kwargs:
                kwargs['ir'] = kwargs.pop('intron_retention')
            if 'type_junction' not in kwargs:
                kwargs['type_junction'] = 0
            junctions[J_id] = JunctionGraphic(**kwargs)
        lsv_type = lsv_attrs.attrs.get('type')
        coords = lsv_attrs.attrs.get('coords')
        output[lsv_id] = LsvGraphic(lsv_type, coords, id = lsv_id, exons = exons, junctions = junctions)
    return output

RANK_TYPES = ['all', 'intersected', 'only_exp1', 'exp1_and_exp2']

def collapse_matrix(matrix):
    xbins, ybins = matrix.shape
    assert xbins == ybins
    DIAG = [matrix.diagonal(offset = xx).sum() for xx in range(1-xbins, xbins)]
    return np.array(DIAG)

def _find_delta_border(V, numbins):
    """Finds the border index to which a V corresponds in its delta_space given the number of bins the matrix will have"""
    delta_space = linspace(-1, 1, num = numbins + 1)[1:]
    return (V > delta_space).sum()

def matrix_area(matrix, V=0.2, absolute=True, collapsed_mat=False):
    """Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute"""
    collapse = matrix if collapsed_mat else collapse_matrix(matrix)
    # get the delta psi histogram borders based on the size of 'collapse'
    # grab the values inside the area of interest
    nbins = collapse.shape[0]
    delta_space = linspace(-1, 1, num = nbins + 1)[1:]
    if absolute:
        delta_space = abs(delta_space)
    if V < 0:
        border = abs(V) > delta_space
    else:
        border = V < delta_space
    area = collapse[border].sum()
    return area


def v_sum(matrix, absolute = True):
    """
    Calculate sum_v v*P(Delta PSI > V)
    """
    vs = arange(0, 1, 0.1)
    excluded = ['matrix', 'absolute', 'collapsed_mat']
    MA = vectorize(matrix_area, [float64], excluded=excluded)
    ret = MA(matrix, V=vs, absolute=absolute).dot(vs)
    return ret


def get_expected_psi(bins):
    bins = np.array(bins)
    step = 1.0 / len(bins)
    return arange(step/2, 1, step).dot(bins)


def expected_dpsi(matrix, collapsed_mat=False, absolute = True):
    """
    Calculate sum_dpsi=Prob(dpsi)*dpsi == sum_v = v*P(Delta PSI)
    """
    collapsed = matrix if collapsed_mat else collapse_matrix(matrix)
    xbins = linspace(-1, 1, num=collapsed.size+1)[:-1] + 1./collapsed.size
    if absolute: xbins = abs(xbins)
    return collapsed.dot(xbins)


def rank_majiq(vlsv_list, V=0.2, absolute=True, dofilter=True, ranknochange=False, 
               prior=None, shrink=True, junc_selection=None, majiq_n=None):

    print "Num of LSVs in majiq: %d" % len(vlsv_list)
    covered_exons = []
    rank = []
    for i, vlsv in enumerate(vlsv_list):
        lsv_bins = vlsv.get_bins()
        if 'i' in vlsv.get_type():
            if max(get_expected_psi(vlsv.psi1[-1]), get_expected_psi(vlsv.psi2[-1])) > 0.1 or len(lsv_bins)<2:
                continue
            lsv_bins = lsv_bins[:-1]

        # Filtering out lsvs that have exons shared with an already added lsv
        lsv_exon_coords = [int(coord) for coord in vlsv.get_id().split(':')[1].split('-')]
        if np.any([ee.get_coords() in covered_exons for ee in vlsv.lsv_graphic.get_exons() if list(ee.get_coords()) <> lsv_exon_coords ]):
            continue
        covered_exons.extend([ee.get_coords() for ee in vlsv.lsv_graphic.get_exons() if list(ee.get_coords()) <> lsv_exon_coords ])

        if len(lsv_bins) > 2:
            if junc_selection:
                junc_n = junc_selection[vlsv.get_id()]
            else:
                e_dpsi = [abs(expected_dpsi(np.array(junc_bins), collapsed_mat=True)) for junc_bins in lsv_bins]
                junc_n = np.argmax(e_dpsi)
        else:
            junc_n = 0
        bins_selected = lsv_bins[junc_n]

        if ranknochange:
            dmatrix = np.exp(np.log(bins_selected))
            dmatrix /= sum(dmatrix)
        else:
            dmatrix = np.array(bins_selected)

        v_expected = expected_dpsi(dmatrix, collapsed_mat=True)
        area = 1.0 - matrix_area(dmatrix, V, absolute, collapsed_mat=True)  # P(Delta PSI < V) = 1 - P(Delta PSI > V)
        rank.append(["%s#%d" % (vlsv.get_id(), junc_n), v_expected, area, int(abs(v_expected)>=V and area<=0.05)])

    expected_mask = np.array([abs(r[1]) >= V for r in rank])
    fdr_mask = np.array([r[2] <= 0.05 for r in rank])
    expected_fdr_mask = np.logical_and(expected_mask, fdr_mask)
    if not majiq_n[0]:
        majiq_n[0] = np.count_nonzero(expected_fdr_mask)

    print "#FDR < 0.05: %d" % np.count_nonzero(fdr_mask)
    print "#E(Delta(PSI))>%.2f: %d" % (V, np.count_nonzero(expected_mask))
    print "#E(Delta(PSI))>%.2f and FDR<0.05: %d" % (V, np.count_nonzero(expected_fdr_mask))

    rank.sort(key=lambda x: (x[1], x[2]), reverse=True)

    return np.array(rank)


def rank_naive(bins_list, names, V=0.2, absolute=True, E=False, ranknochange=False, complex_lsvs=False):
    """Similar to MAJIQ files with the difference that the matrix is already collapsed"""
    rank = []

    covered_exons = []
    print "Num of LSVs in naive_bootstrapping: %d" % len(bins_list)
    for i, lsv_bins in enumerate(bins_list):

        # Filtering out lsvs that have exons shared with an already added lsv
        lsv_exon_coords = [int(coord) for coord in names[i][1].split(':')[1].split('-')]
        if np.any( [ee.get_coords() in covered_exons for ee in names[i][4].get_exons() if list(ee.get_coords()) <> lsv_exon_coords ] ):
            continue
        covered_exons.extend([ee.get_coords() for ee in names[i][4].get_exons() if list(ee.get_coords()) <> lsv_exon_coords ])
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
    return rank, np.count_nonzero(expected_fdr_mask)


def rank_miso(path, dofilter=True, ranknochange=False, complex_lsvs=False):
    rank = scripts.utils.miso_delta_reader(path, dofilter=dofilter, complex_lsvs=complex_lsvs)
    if ranknochange:
        rank.sort(key=lambda x: (abs(x[1]), x[2]))
        # sort first by smallest delta PSI, then by bayes factor
    else:
        rank.sort(key=lambda x: (abs(x[1]), x[2]), reverse=True)
        # sort first by biggest delta PSI, then by inverse bayes factor
    return rank


def rank_mats_original(mats_file, dofilter=True, ranknochange=False, majiq_n=None, V=0.2):
    """Rank Splicing Events detected as differentially expressed by MATS. Compute the FDR for downstream analysis

    :param mats_file:
    :param dofilter:
    :param ranknochange:
    :return: list of events. Events are represented as a list where the 1st pos is the ID and the 2nd the deltapsi.
             NOTE: Since MATS generates IDs that are not unique among runs, we create it as the concatenation of all
             exonic coordinates (3 exons, except special case of MXE with 4).
    """

    # ID        GeneID        geneSymbol        chr        strand        longExonStart_0base        longExonEnd        shortES        shortEE        flankingES        flankingEE        ID        IC_SAMPLE_1        SC_SAMPLE_1        IC_SAMPLE_2        SC_SAMPLE_2        IncFormLen        SkipFormLen        PValue        FDR        IncLevel1        IncLevel2        IncLevelDifference
    # 1292        "ENSMUSG00000032366"        "Tpm1"        9        -        67022592        67024565        67022592        67023441        67031028        67031098        1292        1484,2687,1928        7,12,10        17,11,18        24,23,39        1171        61        0.0        0.0        0.917,0.921,0.909        0.036,0.024,0.023        0.888
    # 194        "ENSMUSG00000025199"        "Chuk"        19        -        44078888        44079266        44078888        44078985        44081945        44081995        194        59,87,60        3,8,2        34,51,61        57,70,38        310        43        0.0        0.0        0.732,0.601,0.806        0.076,0.092,0.182        0.596

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
            pass_thres =int(abs(delta_psi)>=V and fdr<=0.05)
            mats_nn += pass_thres
#            rank.append([geneID, delta_psi, pvalue, fdr, pass_thres])
#            if majiq_n[0]:
#                rank.append([geneID, delta_psi, pvalue, fdr, pass_thres])
#            else:
            if pvalue <= 0.05:
                rank.append([geneID, delta_psi, pvalue, fdr, pass_thres])

    expected_mask = np.array([abs(r[1]) >= V for r in rank])
    fdr_cutoff = 0.05
    if majiq_n[0]:
        while np.count_nonzero(np.logical_and(expected_mask, np.array([r[3] <= fdr_cutoff for r in rank]))) < majiq_n[0] and fdr_cutoff <= 1.0:
            print np.count_nonzero(np.logical_and(expected_mask, np.array([r[3] <= fdr_cutoff for r in rank])))
            fdr_cutoff += 0.05

    else:
        majiq_n[0] = np.count_nonzero(np.array([r[-1] for r in rank]))
    fdr_mask = np.array([r[3] <= fdr_cutoff for r in rank])

    print "MATS:"
    print "#FDR < 0.05: %d" % np.count_nonzero(np.array([r[3] <= 0.05 for r in rank]))
    print "#FDR < %.2f: %d" % (fdr_cutoff, np.count_nonzero(fdr_mask))
    print "#E(Delta(PSI))>%.2f: %d" % (V, np.count_nonzero(expected_mask))
    print "#E(Delta(PSI))>%.2f and FDR<0.05: %d" % (V, np.count_nonzero(np.logical_and(expected_mask, fdr_mask)))

    rank.sort(key=lambda x: (abs(float(x[1])), -float(x[3])))  # biggest delta PSI first, small p-value

    return mats_nn, rank


def rank_dexseq(dexseq_file, dofilter=True, ranknochange=False, majiq_n=None, V=0.2):
    rank = []
    dexseq_nn = 0

    for line in open(dexseq_file).readlines()[1:]:

        sline = line.strip().split()
        if 'NA' in sline:
                continue
        log2fold = float(sline[9])
        padj = float(sline[6])
        ev_name = sline[0]
        pass_thres = int(padj<0.05 and log2fold>4)
        dexseq_nn += pass_thres
        if padj<0.1:
            rank.append([ev_name, log2fold, padj, pass_thres])

    if not majiq_n[0]:
        majiq_n[0] = np.count_nonzero(np.array([r[3] for r in rank]))

    print "DexSeq:"
    print "N:", majiq_n[0]

    rank.sort(key=lambda x: (abs(float(x[1])), -float(x[2])))  # biggest delta PSI first, small p-value

    return dexseq_nn, rank


def rank_suppa2(suppa_file, dofilter=True, ranknochange=False, majiq_n=None, V=0.2):
    rank = []
    suppa_nn = 0
    print "LL"
    for line in open(suppa_file).readlines()[1:]:
        sline = line.strip().split()
        if sline[1] == 'nan':
            continue
        dpsi = float(sline[1])
        pval = float(sline[2])
        ev_name = sline[0]
        pass_thres = int(abs(dpsi)>=V and pval<=0.05)
        suppa_nn += pass_thres
#        if pval<=0.05:
#            rank.append([ev_name, dpsi, pval, pass_thres])
#        if majiq_n[0]:
#            if abs(dpsi)< V:
#                continue
#            rank.append([ev_name, dpsi, pval, pass_thres])
#        else:
        if pval <= 0.05:
                rank.append([ev_name, dpsi, pval, pass_thres])

    expected_mask = np.array([abs(r[1]) >= V for r in rank])
    pval_cutoff = 0.05
    if not majiq_n[0]:
        majiq_n[0] = np.count_nonzero(np.array([r[3] for r in rank]))
    fdr_mask = np.array([r[2] <= pval_cutoff for r in rank])

    print "Suppa2:"
    print "#FDR < 0.05: %d" % np.count_nonzero(fdr_mask)

    print "#FDR < %.2f: %d" % (pval_cutoff, np.count_nonzero(fdr_mask))
    print "#E(Delta(PSI))>%.2f: %d" % (V, np.count_nonzero(expected_mask))
    print "#E(Delta(PSI))>%.2f and FDR<0.05: %d" % (V, np.count_nonzero(np.logical_and(expected_mask, fdr_mask)))

    rank.sort(key=lambda x: (abs(float(x[1])), -float(x[2])))  # biggest delta PSI first, small p-value

    return suppa_nn, rank


def _is_in_chunk(event1, chunk, report_rank2_expec=False):
    for event2 in chunk:
        if event1[0] == event2[0]:  # event[0] is the name of the event
            if report_rank2_expec:
                return [1, event2[1]]
            return 1
    if report_rank2_expec:
        return [0, '']
    return 0


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
    parser.add_argument('--dexseq-files', dest='dexseq_files', nargs=2,
                        help='Dexseq files with paired events to analyze')
    parser.add_argument('--suppa-files', dest='suppa_files', nargs=2, help='SUPPA2 files with paired events to analyze')
    parser.add_argument('--naive-files', dest='naive_files', nargs=2,
                        help='Naive Bootstrapping files with paired events to analyze')
    parser.add_argument('-o', '--output', required=True, help='Output file_str path')
    parser.add_argument('--max', default=1000, type=int, help="Max number of events to analyze")
    parser.add_argument('--V', default=0.2, type=float, help="Steps of best events to take")
    parser.add_argument('--E', default=False, action="store_true", help="Deprecated")
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
    parser.add_argument('--splice-graphics', type = h5py.File, help = "HDF5 file containing relevant splice graphics")
    parser.add_argument('--majiq_N', type = int, nargs = 1, default = [None])
    args = parser.parse_args()

    print args

    if args.splice_graphics:
        print "Parsing splice graphics into legacy format"
        if 'LSVs' in args.splice_graphics.keys():
            graphs = legacy_lsv_junctions(args.splice_graphics)
        else:
            graphs = legacy_lsv_junctions_new_voila(args.splice_graphics)

    print "Calculating ranks..."
    ranks = defaultdict(list)
    n1 = defaultdict(list)
    method_N = {}

    if args.majiq_files:
        method_N['majiq'] = [args.majiq_N]
        count_pairs = 0
        for file_nr, file_str in enumerate(args.majiq_files):
            try:
                majiq_data = pickle.load(open(file_str, 'r'))
            except:
                with h5py.File(file_str, 'r') as IN:
                    majiq_data = VoilaInput(IN, graphs)
            if file_nr % 2 == 0:
                count_pairs += 1
                majiq_rank = rank_majiq(majiq_data.get_lsvs(), args.V, args.absolute, args.filter, args.ranknochange,
                                        args.complex_lsvs, shrink=args.shrink, majiq_n=method_N['majiq'])
                ranks['majiq_' + str(count_pairs)].append(majiq_rank)
                majiq_file1_names = [a[0].split('#')[0] for a in ranks['majiq_1'][0]]
                print "MAJIQ 1st delta(psi) ALL pair %d" % len(majiq_data.get_lsvs())
                print "MAJIQ 1st delta(psi) AFTER RANK pair %d" % len(majiq_file1_names)
                continue

            if args.type_rank == 'only_exp1':
                # Select events from experiment 1
                exp1_index = np.array([v_lsv.get_id() in majiq_file1_names for v_lsv in majiq_data.get_lsvs()])
                junc_dict = dict([(rr[0].split('#')[0], int(rr[0].split('#')[1])) for rr in ranks['majiq_' + str(count_pairs)][-1]])
                majiq_rank = rank_majiq(np.array(majiq_data.get_lsvs())[exp1_index].tolist(), args.V, args.absolute,
                                        args.filter, args.ranknochange, args.complex_lsvs, shrink=args.shrink,
                                        junc_selection=junc_dict)
                ranks['majiq_' + str(count_pairs)].append(majiq_rank)
                n1['majiq_' + str(count_pairs)] = [np.count_nonzero(exp1_index), np.count_nonzero(exp1_index)]
                print "MAJIQ 2nd delta(psi) pair %d" % np.count_nonzero(exp1_index)
        names_majiq_exp1 = [m[1] for m in majiq_file1_names]

    mats_n_orig = None

    if args.mats_files:
        method_N['mats'] = [None]
        for file_str in args.mats_files:
            mmats_n, mats_rank = rank_mats_original(file_str, args.filter, args.ranknochange, majiq_n=method_N['mats'],
                                                    V=args.V)#majiq_N)
            if not mats_n_orig:
                mats_n_orig = mmats_n
            if len(ranks['mats']) > 0:
                names_mats_exp1 = [mats_info[0] for mats_info in ranks['mats'][0]]
                names_mats = [mats_info[0] for mats_info in mats_rank]
                exp1_index = np.array([name in names_mats_exp1 for name in names_mats])
                ranks['mats'].append(np.array(mats_rank)[exp1_index])
                n1['mats'].append(np.count_nonzero(exp1_index))
                print "MATS 2nd delta(psi) pair: %d" % np.count_nonzero(exp1_index)
            else:
                print "MATS 1st delta(psi) pair: %d" % len(mats_rank)
                ranks['mats'].append(array(mats_rank))

    if args.suppa_files:

        method_N['suppa'] = [None]
        for file_str in args.suppa_files:
            msuppa_n, suppa_rank = rank_suppa2(file_str, args.filter, args.ranknochange,
                                               majiq_n=method_N['suppa'], V=args.V)
            if not mats_n_orig:
                mats_n_orig = msuppa_n
            if len(ranks['suppa']) > 0:
                names_suppa_exp1 = [suppa_info[0] for suppa_info in ranks['suppa'][0]]
                names_suppa = [suppa_info[0] for suppa_info in suppa_rank]
                exp1_index = np.array([name in names_suppa_exp1 for name in names_suppa])
                ranks['suppa'].append(np.array(suppa_rank)[exp1_index])
                n1['suppa'].append(np.count_nonzero(exp1_index))
                print "SUPPA 2nd delta(psi) pair: %d" % np.count_nonzero(exp1_index)
            else:
                print "SUPPA 1st delta(psi) pair: %d" % len(suppa_rank)
                ranks['suppa'].append(array(suppa_rank))

    if args.dexseq_files:

        method_N['dexseq'] = [None]
        for file_str in args.dexseq_files:
            mdexseq_n, dexseq_rank = rank_dexseq(file_str, args.filter,
                                                 args.ranknochange,
                                                 majiq_n=method_N['dexseq'], V=args.V)
            if not mats_n_orig:
                mats_n_orig = mdexseq_n
            if len(ranks['dexseq']) > 0:
                names_dexseq_exp1 = [dexseq_info[0] for dexseq_info in
                                     ranks['dexseq'][0]]
                names_dexseq = [dexseq_info[0] for dexseq_info in dexseq_rank]
                exp1_index = np.array([name in names_dexseq_exp1 for name in
                                       names_dexseq])
                ranks['dexseq'].append(np.array(dexseq_rank)[exp1_index])
                n1['dexseq'].append(np.count_nonzero(exp1_index))
                print "DexSeq 2nd delta(psi) pair: %d" % np.count_nonzero(exp1_index)
            else:
                print "DexSeq 1st delta(psi) pair: %d" % len(dexseq_rank)
                ranks['dexseq'].append(array(dexseq_rank))

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
        naive_N = None
        for file_str in args.naive_files:
            naive_data = pickle.load(open(file_str, 'r'))
            naive_rank, naive_n_tmp = rank_naive(naive_data[1], naive_data[0], args.V, args.absolute, args.E, args.ranknochange)
            if not naive_N:
                naive_N = naive_n_tmp
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
        max_events = min(max_events, method_N[method_name][0])

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
                events.append([rank1[i], is_hit, rank2_exp, method_N[method_name][0] if 'mats' not in method_name else mats_n_orig])

            fdr.pop(0)  #remove now useless first item
            #normalize ratios
            ratios = array(ratios)
            ratios /= ratios.shape[0]
            if args.fdr:  #normalize fdr if we are calculating it
                fdr = array(fdr)
                fdr /= fdr.shape[0]

        else:  # "equalrank" chunks of same n size in both ranks
            N, ratios = get_ratio_nu(rank1, rank2, max_events)

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


        print "Saving rank1... in %s " % args.output
        pickle.dump(rank1,
                    open(args.output + "/rank1.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name),
                         'w'))

        print "Saving rank2... in %s " % args.output
        pickle.dump(rank2,
                    open(args.output + "/rank2.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name),
                         'w'))


        print "Saving ratios_nu... in %s " % args.output
        pickle.dump((N, ratios),
                    open(args.output + "/ratios_nu.%s.%s.pickle" % (str(args.type_rank).replace('-', '_'), method_name),
                         'w'))

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
