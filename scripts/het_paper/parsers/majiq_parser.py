from voila.io_voila import VoilaInput
import numpy as np


def _find_delta_border(V, numbins):
    """Finds the border index to which a V corresponds in its delta_space given the number of bins the matrix will have"""
    delta_space = list(np.linspace(-1, 1, num=numbins + 1))
    delta_space.pop(0)  # first border to the left is -1, and we are not interested in it
    # get the index position that corresponds to the V threshold
    for i, value in enumerate(delta_space):
        if value > V:
            return i
    # if nothing hit, V = 1
    return numbins


def prob_mat_changing_above_threshold(matrix, V=0.2, absolute=True):
    """Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute"""

    # get the delta psi histogram borders based on the size of 'collapse'
    border = _find_delta_border(V, matrix.shape[0])
    # grab the values inside the area of interest
    area = []
    if V < 0:
        area.append(matrix[0:border + 1])
        if absolute:  # if absolute V, pick the other side of the array
            area.append(matrix[-border - 1:])
    else:
        area.append(matrix[border:])
        if absolute and border != 0:  # if absolute V, pick the other side of the array
            area.append(matrix[0:len(matrix) - border])

    return sum(area)


def rank_majiq(vlsv_list, V=0.2, prob=True, complex_lsvs=False):

    MINTHRESHOLD = 0.95
    if E:
        MINTHRESHOLD = 0.20
    rank = []

    print "Num of LSVs in old_majiq: %d" % len(vlsv_list)
    covered_exons = []
    for i, vlsv in enumerate(vlsv_list):
        lsv_bins = vlsv.bins
        junc_n = 0
        if 'i' in vlsv.lsv_graphic.lsv_type:
            if max(vlsv.means_psi1[-1], vlsv.means_psi2[-1]) > 0.1 or len(lsv_bins) < 2:
                continue
            lsv_bins = lsv_bins[:-1]

        # Filtering out lsvs that have exons shared with an already added lsv
        lsv_exon_coords = [int(coord) for coord in vlsv.lsv_graphic.gene_id.split(':')[1].split('-')]
        if np.any([(ee.start, ee.end) in covered_exons for ee in vlsv.lsv_graphic.exons if list((ee.start, ee.end)) <> lsv_exon_coords]):
            continue
        covered_exons.extend([(ee.start, ee.end) for ee in vlsv.lsv_graphic.exons if list((ee.start, ee.end)) <> lsv_exon_coords ])
        if len(lsv_bins) > 2:
            most_change = 0
            for jj, junc_bins in enumerate(lsv_bins):
                if abs(vlsv.means[jj]) > most_change:
                    bins_selected = junc_bins
                    junc_n = jj
        else:
            bins_selected = lsv_bins[0]

        dmatrix = np.array(bins_selected)

        v_expected = vlsv.means[0]
        area = 1.0 - prob_mat_changing_above_threshold(dmatrix, V)
        # P(Delta PSI < V) = 1 - P(Delta PSI > V)

    rank.append(["%s#%d" % (vlsv.lsv_graphic.gene_id, junc_n), v_expected, area, int(abs(v_expected)>=V and area<=0.05)])

    expected_mask = np.array([abs(r[1]) >= V for r in rank])
    fdr_mask = np.array([r[2] <= 0.05 for r in rank])
    expected_fdr_mask = np.logical_and(expected_mask, fdr_mask)
    # rank = np.array(rank)[expected_fdr_mask].tolist()
    if majiq_n:
        majiq_n[0] = np.count_nonzero(expected_fdr_mask)

    print "#FDR < 0.05: %d" % np.count_nonzero(fdr_mask)
    print "#E(Delta(PSI))>%.2f: %d" % (V, np.count_nonzero(expected_mask))
    print "#E(Delta(PSI))>%.2f and FDR<0.05: %d" % (V, np.count_nonzero(expected_fdr_mask))

    rank.sort(key=lambda x: (x[1], x[2]), reverse=True)

    return np.array(rank)


def parse(filen1, filen2, threshold):

    majiq_data1 = VoilaInput.from_hdf5_file(hdf5_filename=filen1)
    majiq_data2 = VoilaInput.from_hdf5_file(hdf5_filename=filen2)

    ranks1 = rank_majiq(majiq_data1.get_lsvs(), threshold, shrink=True)
    ranks2 = rank_majiq(majiq_data1.get_lsvs(), threshold, shrink=True)


if __name__ == '__main__':

    import sys
    parse(sys.argv[1])
