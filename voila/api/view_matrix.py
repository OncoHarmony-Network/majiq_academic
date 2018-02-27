import math
from itertools import zip_longest

import numpy

from voila import constants
from voila.api.matrix_hdf5 import DeltaPsi, Psi
from voila.exceptions import NoLsvsFound
from voila.utils.voila_log import voila_log
from voila.vlsv import get_expected_dpsi, matrix_area, is_lsv_changing


def unpack_means(value):
    if numpy.size(value, 0) == 1:
        value = numpy.append(value, numpy.array(1 - value[0]))
    return value


def unpack_bins(value):
    if numpy.size(value, 0) == 1:
        value = numpy.append(value, [numpy.flip(value[-1], 0)], axis=0)
    return value


class ViewPsi(Psi):
    class _ViewPsi(Psi._Psi):
        def get_all(self):
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'reference_exon', tuple(self.reference_exon)
            yield 'gene_id', self.gene_id
            yield 'is_target', self.is_target

            fields = list(self.fields)

            fields.remove('bins')
            yield 'group_bins', dict(self.group_bins)

            fields.remove('means')
            yield 'group_means_rounded', dict(self.group_means_rounded)

            yield from self.get_many(fields)

        @property
        def junctions(self):
            return self.get('junctions')

        @property
        def means(self):
            yield from unpack_means(self.get('means'))

        @property
        def group_means(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.means

        @property
        def group_means_rounded(self):
            for group_name, means in self.group_means:
                yield group_name, numpy.around(tuple(means), decimals=3)

        @property
        def bins(self):
            return unpack_bins(self.get('bins'))

        @property
        def group_bins(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.bins

        @property
        def variances(self):
            def get_expected_psi(bins):
                step = 1.0 / bins.size
                projection_prod = bins * numpy.arange(step / 2, 1, step)
                return numpy.sum(projection_prod)

            for b in self.bins:
                epsi = get_expected_psi(b)
                step_bins = 1.0 / b.size
                projection_prod = b * numpy.arange(step_bins / 2, 1, step_bins) ** 2
                yield numpy.sum(projection_prod) - epsi ** 2

        @property
        def junction_count(self):
            return len(tuple(self.means))

    def psi(self, lsv_id):
        return self._ViewPsi(self, lsv_id)


class ViewDeltaPsi(DeltaPsi):
    class _ViewDeltaPsi(DeltaPsi._DeltaPsi):
        def get_all(self):
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'reference_exon', self.reference_exon
            yield 'excl_incl', self.excl_incl
            yield 'is_target', self.is_target

            fields = list(self.fields)

            fields.remove('group_bins')
            yield 'group_bins', dict(self.group_bins)

            fields.remove('group_means')
            yield 'group_means_rounded', dict(self.group_means_rounded)

            fields.remove('bins')
            yield 'bins', self.bins
            yield 'means_rounded', self.means_rounded

            yield from self.get_many(fields)

        @property
        def junctions(self):
            return self.get('junctions')

        @property
        def bins(self):
            return unpack_bins(self.get('bins'))

        @property
        def group_bins(self):
            group_names = self.matrix_hdf5.group_names
            group_bins = self.get('group_bins')
            for group_name, value in zip(group_names, group_bins):
                yield group_name, unpack_bins(value)

        @property
        def means(self):
            for b in self.bins:
                yield get_expected_dpsi(b)

        @property
        def means_rounded(self):
            return numpy.around(tuple(self.means), decimals=3)

        @property
        def group_means(self):
            for value in self.get('group_means'):
                yield unpack_means(value)

        @property
        def group_means_rounded(self):
            group_names = self.matrix_hdf5.group_names
            for group_name, means in zip(group_names, self.group_means):
                yield group_name, numpy.around(means, decimals=3)

        @property
        def excl_incl(self):
            for mean in self.means:
                if mean < 0:
                    yield [-mean, 0]
                else:
                    yield [0, mean]

        @property
        def junction_count(self):
            return len(tuple(self.means))

    def delta_psi(self, lsv_id):
        return self._ViewDeltaPsi(self, lsv_id)


class ViewMatrix:
    def __init__(self, matrix):
        self.matrix = matrix

    def view_lsv_count(self, args):
        value = len(tuple(self.view_lsv_ids(args)))
        voila_log().info('Found {} LSVs'.format(value))
        if not value:
            raise NoLsvsFound()
        return value

    def view_gene_lsvs(self, args, gene_id):
        yield from self.valid_lsvs(args, self.matrix.lsv_ids([gene_id]))

    def view_lsv_ids(self, args):
        if args.lsv_ids:
            lsv_ids = args.lsv_ids
        elif args.gene_ids:
            lsv_ids = self.matrix.lsv_ids(args.gene_ids)
        else:
            lsv_ids = self.matrix.lsv_ids()

        yield from self.valid_lsvs(args, lsv_ids)

    def view_gene_ids(self, args):
        if args.gene_ids:
            gene_ids = args.gene_ids
        else:
            gene_ids = self.matrix.gene_ids

        for gene_id in gene_ids:
            if any(self.valid_lsvs(args, self.view_gene_lsvs(args, gene_id))):
                yield gene_id

    @property
    def metadata(self):
        metadata = super().metadata
        experiment_names = metadata['experiment_names']
        group_names = super().group_names

        if experiment_names.size > 1:
            metadata['experiment_names'] = [numpy.insert(exps, 0, '{0} Combined'.format(group)) for exps, group in
                                            zip(experiment_names, group_names)]

        return metadata

    def paginated_genes(self, args):
        def grouper(iterable, n, fillvalue=None):
            args = [iter(iterable)] * n
            return zip_longest(*args, fillvalue=fillvalue)

        for page in grouper(self.view_gene_ids(args), constants.MAX_GENES):
            yield tuple(p for p in page if p is not None)

    def page_count(self, args):
        gene_count = len(tuple(self.view_gene_ids(args)))
        return math.ceil(gene_count / constants.MAX_GENES)


class ViewPsiMatrix(ViewMatrix):
    def valid_lsvs(self, args, lsv_ids):
        threshold = None
        percent_threshold = None
        m_a = matrix_area

        if hasattr(args, 'show_all') and not args.show_all:
            threshold = args.threshold
            percent_threshold = args.percent_threshold

        for lsv_id in lsv_ids:
            dpsi = self.matrix.psi(lsv_id)
            if not args.lsv_ids or lsv_id in args.lsv_ids:
                if not threshold or is_lsv_changing(dpsi.means, threshold):
                    if not percent_threshold or any(m_a(b, collapsed_mat=True) >= percent_threshold for b in dpsi.bins):
                        yield lsv_id


class ViewDeltaPsiMatrix(ViewPsiMatrix):
    def valid_lsvs(self, args, lsv_ids):
        threshold = None
        percent_threshold = None
        m_a = matrix_area

        if hasattr(args, 'show_all') and not args.show_all:
            threshold = args.threshold
            percent_threshold = args.percent_threshold

        for lsv_id in lsv_ids:
            dpsi = self.matrix.delta_psi(lsv_id)
            if not args.lsv_ids or lsv_id in args.lsv_ids:
                if not threshold or is_lsv_changing(dpsi.means, threshold):
                    if not percent_threshold or any(m_a(b, collapsed_mat=True) >= percent_threshold for b in dpsi.bins):
                        yield lsv_id
