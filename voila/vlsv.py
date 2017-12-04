import os
from collections import defaultdict

import numpy
import numpy as np
import scipy.interpolate as sinter

from voila.hdf5 import BinsDataSet, Psi1DataSet, Psi2DataSet
from voila.hdf5 import HDF5
from voila.splice_graphics import LsvGraphic


class OrphanJunctionException(Exception):
    def __init__(self, m):
        self.message = m


class Het(HDF5):
    def __init__(self):
        """
        Het stat data. 
        :param experiment_names: List of experiment names. 
        """
        super(Het, self).__init__()
        self.groups = []
        self.junction_stats = []

    def add_group(self, expected_psi, median):
        """
        Add per LSV het group stats. 
        :param expected_psi: Expected psi 2d array, where x axis is experiments and y axis is junctions
        :param median: Median psi bins
        :param group_name: The name of the group being added
        :return: 
        """
        self.groups.append(HetGroup(expected_psi, median))

    def add_junction_stats(self, stats):
        """
        Add per junction het stat data.
        :param stats: 2d array of group comparison stats 
        :param stat_name: Name of stat
        :param junction_id: Junction id
        :return: 
        """
        self.junction_stats.append(stats)

    def cls_dict(self):
        return {'groups': HetGroup}

    def __iter__(self):
        for k, v in self.__dict__.items():
            yield k, v


class HetGroup(HDF5):
    def __init__(self, expected_psi, median):
        super(HetGroup, self).__init__()
        self.expected_psi = expected_psi
        self.median = median

    def get_psi(self, experiment_index, junction_index):
        return self.expected_psi[experiment_index][junction_index]

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None).from_hdf5(h)

    def __iter__(self):
        for k, v in self.__dict__.items():
            yield k, v


class VoilaLsv(LsvGraphic):
    def __init__(self, lsv_graphic, bins_list=None, means=None, means_psi1=None, psi1=None, means_psi2=None, psi2=None,
                 het=None):
        """
        The lsv data voila needs to visualize LSVs.
        :param bins_list: bins matrix
        :param lsv_graphic: lsv graphic object which is used to create this object
        :param means: means data
        :param means_psi1: means data for psi1
        :param psi1: psi1 matrix
        :param means_psi2: means data for psi2
        :param bin: psi2 matrix
        """

        if lsv_graphic:
            super(VoilaLsv, self).__init__(
                lsv_type=lsv_graphic.lsv_type,
                start=lsv_graphic.start,
                end=lsv_graphic.end,
                lsv_id=lsv_graphic.lsv_id,
                name=lsv_graphic.name,
                strand=lsv_graphic.strand,
                exons=lsv_graphic.exons,
                junctions=lsv_graphic.junctions,
                chromosome=lsv_graphic.chromosome
            )
        else:
            super(VoilaLsv, self).__init__(None, None, None, None)

        # For binary LSVs, we only store one of the junctions... use their property for the full value.
        self.trunc_bins = bins_list
        self.trunc_psi1 = psi1
        self.trunc_psi2 = psi2
        self.trunc_means = means
        self.trunc_means_psi1 = means_psi1
        self.trunc_means_psi2 = means_psi2

        self.het = het

        self.categories = None
        self.excl_incl = []

        if lsv_graphic:
            # Store collapsed matrix to save some space
            if self.is_delta_psi():
                self.trunc_bins = [collapse_matrix(bins) for bins in self.trunc_bins]

            # excl_incl is used to calculate the expected psi
            if self.is_delta_psi() and self.means:
                for mean in self.means:
                    if mean < 0:
                        self.excl_incl.append([-mean, 0])
                    else:
                        self.excl_incl.append([0, mean])

            self.init_categories()

    @staticmethod
    def _extend_means(means):
        if means is not None:
            means = numpy.array(means)
            if numpy.size(means, 0) == 1:
                # means.append(1 - means[0])
                means = numpy.append(means, numpy.array(1 - means[0]))
        return means

    @staticmethod
    def _extend_bins(bins):
        if bins is not None:
            bins = numpy.array(bins)
            if numpy.size(bins, 0) == 1:
                # bins.append(bins[-1][::-1])
                bins = numpy.append(bins, [numpy.flip(bins[-1], 0)], axis=0)
        return bins

    @property
    def bins(self):
        return self._extend_bins(self.trunc_bins)

    @property
    def psi1(self):
        return self._extend_bins(self.trunc_psi1)

    @property
    def psi2(self):
        return self._extend_bins(self.trunc_psi2)

    @property
    def means(self):
        ms = self._extend_means(self.trunc_means)
        if self.bins is not None:
            if ms is None and self.bins.size > 0:
                if self.is_delta_psi():
                    ms = [get_expected_dpsi(b) for b in self.bins]
                else:
                    ms = [get_expected_psi(b) for b in self.bins]

        return ms

    @property
    def means_psi1(self):
        return self._extend_means(self.trunc_means_psi1)

    @property
    def means_psi2(self):
        return self._extend_means(self.trunc_means_psi2)

    @property
    def variances(self):
        variances = []
        if self.bins is not None and self.bins.size > 0:
            for bin in self.bins:
                step_bins = 1.0 / len(bin)
                projection_prod = bin * np.arange(step_bins / 2, 1, step_bins) ** 2
                variances.append(np.sum(projection_prod) - self.means[-1] ** 2)
        return variances

    def to_dict(self, ignore_keys=()):
        ignore_keys = ('trunc_bins', 'trunc_psi1', 'trunc_psi2', 'trunc_means', 'trunc_means_psi1', 'trunc_means_psi2')
        d = super(VoilaLsv, self).to_dict(ignore_keys)
        d.update({
            'bins': self.bins,
            'psi1': self.psi1,
            'psi2': self.psi2,
            'means': self.means,
            'means_psi1': self.means_psi1,
            'means_psi2': self.means_psi2,
            'variances': self.variances
        })
        return d

    def is_delta_psi(self):
        """
        Return true if this LSV is delta psi.
        :return: bool
        """
        # we have to check for None, because they might be numpy objects
        return self.trunc_psi1 is not None and self.trunc_psi2 is not None

    def get_extension(self):
        """
        Returns start and stop coordinates for this LSV.
        :return:
        """
        return [self.exons[0].start, self.exons[-1].end]

    def njuncs(self):
        """
        Return the number of junctions from categories. This is used in rendering some HTML files.
        :return: int
        """
        return self.categories['njuncs']

    def nexons(self):
        """
        Return the number of exons from categories.  This is used in rendering some HTML files.
        :return: int
        """
        return self.categories['nexons']

    def categories2css(self):
        """
        Return css class names used when rendering some HTML files.
        :return: str
        """
        css_cats = []
        for c in self.categories:
            if type(self.categories[c]) == bool and self.categories[c]:
                css_cats.append(c)
        return ' '.join(css_cats)

    @staticmethod
    def is_lsv_changing(means, threshold):
        """
        Return true if lsv is changing based on threshold.
        :param threshold: lsv threshold value
        :return: bool
        """
        means = np.array(means)
        return max(means[means > 0].sum(), abs(means[means < 0].sum())) >= threshold

    def exclude(self):
        return ['categories', 'trunc_bins', 'trunc_psi1', 'trunc_psi2', 'het']

    def to_hdf5(self, h, use_id=True):
        if use_id:
            h = h.create_group(os.path.join('/lsvs', self.lsv_id.split(':')[0], self.lsv_id))

        super(VoilaLsv, self).to_hdf5(h, use_id)

        # categories
        cat_grp = h.create_group('categories')
        for key in self.categories:
            cat_grp.attrs[key] = self.categories[key]

        # bins
        BinsDataSet(h, self.trunc_bins).encode_list()

        if self.is_delta_psi():
            # psi1
            Psi1DataSet(h, self.trunc_psi1).encode_list()

            # psi2
            Psi2DataSet(h, self.trunc_psi2).encode_list()

        # het groups
        # TODO: Revisit this solution
        if hasattr(self, 'het') and self.het:
            het_groups_grp = h.create_group('het')
            self.het.to_hdf5(het_groups_grp)

    def from_hdf5(self, h):
        # categories
        self.categories = {}
        cat_attrs = h['categories'].attrs
        for key, value in cat_attrs.items():
            self.categories[key] = value

        # bins
        self.trunc_bins = BinsDataSet(h).decode_list()

        # psi1
        self.trunc_psi1 = Psi1DataSet(h).decode_list()

        # psi2
        self.trunc_psi2 = Psi2DataSet(h).decode_list()

        try:
            self.het = Het.easy_from_hdf5(h['het'])
        except KeyError:
            self.het = None

        return super(VoilaLsv, self).from_hdf5(h)

    def to_gff3(self):
        """
        Format class data as GFF3.
        :return: str
        """

        def find_exon_a5(lexonG, jidx):
            for eG in lexonG:
                if jidx in eG.a5:
                    return eG
            raise OrphanJunctionException("Orphan junction %s in lsv %s." % (
                repr(self.junctions[jidx].coords), lsvId))

        def find_exon_a3(lexonG, jidx):
            for eG in lexonG:
                if jidx in eG.a3:
                    return eG
            raise OrphanJunctionException("Orphan junction %s in lsv %s." % (
                repr(self.junctions[jidx].coords), lsvId))

        # fields = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        trans = []
        lexons = self.exons
        lsvId = self.lsv_id
        chrom = self.chromosome
        strand = self.strand
        gStart = lexons[0].coords()[0]
        gEnd = lexons[-1].coords()[1]
        first = 0
        last = 1
        if strand == '-':
            gStart = lexons[-1].coords()[1]
            gEnd = lexons[0].coords()[0]
            first, last = last, first

        gene_str = '\t'.join([chrom, 'old_majiq', 'gene', str(gStart), str(gEnd), '.', strand, '0',
                              'Name=%s;ID=%s' % (lsvId, lsvId)])

        trans.append(gene_str)
        for jid, junc in enumerate(self.junctions):
            mrna = '%s\told_majiq\tmRNA\t' % chrom
            mrna_id = '%s.%d' % (lsvId, jid)
            ex1 = '%s\told_majiq\texon\t' % chrom
            ex2 = '%s\told_majiq\texon\t' % chrom

            ex1G = find_exon_a5(lexons, jid)
            ex2G = find_exon_a3(lexons, jid)

            if strand == '-':
                ex1G, ex2G = ex2G, ex1G

            mrna += '%d\t%d\t' % (ex1G.coords()[first], ex2G.coords()[last])

            if self.lsv_type.startswith('t'):
                ex1G, ex2G = ex2G, ex1G
                ex1 += '%d\t%d\t' % (junc.coords()[last], ex1G.coords()[last])
                ex2 += '%d\t%d\t' % (ex2G.coords()[first], junc.coords()[first])
            else:
                ex1 += '%d\t%d\t' % (ex1G.coords()[first], junc.coords()[first])
                ex2 += '%d\t%d\t' % (junc.coords()[last], ex2G.coords()[last])
            mrna += '.\t%s\t0\tName=%s;Parent=%s;ID=%s' % (strand, mrna_id, lsvId, mrna_id)
            ex1 += '.\t%s\t0\tName=%s.lsv;Parent=%s;ID=%s.lsv' % (strand, mrna_id, mrna_id, mrna_id)
            ex2 += '.\t%s\t0\tName=%s.ex;Parent=%s;ID=%s.ex' % (strand, mrna_id, mrna_id, mrna_id)
            trans.append(mrna)

            trans.append(ex1)
            trans.append(ex2)
        else:
            if strand == '-':
                for ii, t in enumerate(trans):
                    t_fields = t.split()
                    t_fields[3], t_fields[4] = t_fields[4], t_fields[3]
                    trans[ii] = '\t'.join(t_fields)

            lsv_gtf = '\n'.join(trans)

        return lsv_gtf

    def init_categories(self):
        """
        Create categories values for LSV.
        :return: None
        """
        categories = defaultdict()
        juns = self.lsv_type.split('|')
        ir = 'i' in juns

        try:
            juns.remove('i')
        except ValueError:
            # No sign of intron retention
            pass

        ssites = set(int(s[0]) for s in juns[1:])
        exons = defaultdict(list)

        for s in juns[1:]:
            exs = s[1:].split('.')

            try:
                ssite = int(exs[1].split('o')[0])
                exons[exs[0]].append(ssite)
            except IndexError:
                pass

        categories['ES'] = len(exons.keys()) > 1
        categories['prime5'] = len(ssites) > 1
        categories['prime3'] = max([len(exons[e]) for e in exons]) > 1
        categories['njuncs'] = np.sum(['e0' not in junc for junc in juns[1:]])
        categories['nexons'] = len(exons.keys()) + 1
        categories['source'] = juns[0] == 's'
        categories['target'] = juns[0] == 't'
        categories['ir'] = ir
        if juns[0] == 't':
            categories['prime5'], categories['prime3'] = categories['prime3'], categories['prime5']

        self.categories = categories

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None).from_hdf5(h)


def get_expected_dpsi(bins):
    return sum(np.array(bins) * np.arange(-1 + 1. / len(bins), 1., 2. / len(bins)))


def get_expected_psi(bins):
    bins = np.array(bins)
    step = 1.0 / bins.size
    projection_prod = bins * np.arange(step / 2, 1, step)
    return np.sum(projection_prod)


def collapse_matrix(matrix):
    """
    Collapse the diagonals probabilities in 1-D and return them
    """
    collapse = []
    matrix = np.array(matrix)
    matrix_corner = matrix.shape[0]
    for i in range(-matrix_corner + 1, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)


def _find_delta_border(V, numbins):
    """
    Finds the border index to which a V corresponds in its delta_space given the number of bins the matrix will have
    :param V:
    :param numbins:
    :return:
    """
    delta_space = list(np.linspace(-1, 1, num=numbins + 1))
    delta_space.pop(0)  # first border to the left is -1, and we are not interested in it
    # get the index position that corresponds to the V threshold
    for i, value in enumerate(delta_space):
        if value > V:
            return i
    # if nothing hit, V = 1
    return numbins


def matrix_area(matrix, V=0.2, absolute=True, collapsed_mat=False):
    """
    Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute.
    :param collapsed_mat:
    :param V:
    :param absolute:
    :param matrix:
    :return:
    """
    collapse = matrix
    if not collapsed_mat:
        collapse = collapse_matrix(matrix)
    collapse = np.concatenate(([0], collapse))
    collapse = np.cumsum(collapse)
    xbins = np.linspace(0, 1, num=collapse.size)
    if absolute:
        Vabs = abs(V)
        left, right = np.interp([-Vabs, V], xbins, collapse, left=0, right=1)
        area = left + (1 - right)
    else:
        area = np.interp(V, xbins, collapse, left=0, right=1)
        if V >= 0:
            area = 1 - area
    return area


def matrix_area_spline(matrix, V=0.2, absolute=True, collapsed_mat=False, kind=2):
    """
    Returns the probability of an event to be above a certain threshold. The absolute flag describes if the value is absolute.
    :param collapsed_mat: if True, matrix must be 1d, else matrix must be 2d square
    :param V: threshold
    :param absolute: if True, calculate two-tailed area
    :param matrix: list or array of probabilities, must sum to 1
    :param kind: parameter kind from scipy.interpolate.interp1d
    :return: probability of delta psi exceeding a certain threshold
    """
    collapse = matrix
    if not collapsed_mat:
        collapse = collapse_matrix(matrix)
    collapse = np.concatenate(([0], collapse))
    collapse = np.cumsum(collapse)
    xbins = np.linspace(-1, 1, num=collapse.size)
    spline = sinter.interp1d(xbins, collapse, kind=kind, fill_value=(0, 1), bounds_error=False)
    if absolute:
        Vabs = abs(V)
        left, right = spline([-Vabs, V])
        area = left + (1 - right)
    else:
        area = spline(V)
        if V >= 0:
            area = 1 - area
    return float(area)
