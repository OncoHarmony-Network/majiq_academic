import math
from abc import ABC, abstractmethod
from itertools import zip_longest

import numpy as np
import scipy.special

from voila import constants
from voila.api.matrix_hdf5 import DeltaPsi, Psi, Heterogen
from voila.exceptions import NoLsvsFound
from voila.utils.voila_log import voila_log
from voila.vlsv import get_expected_dpsi, is_lsv_changing, matrix_area


def unpack_means(value):
    if np.size(value, 0) == 1:
        value = np.append(value, np.array(1 - value[0]))
    return value


def unpack_bins(value):
    if np.size(value, 0) == 1:
        value = np.append(value, [np.flip(value[-1], 0)], axis=0)
    return value


def passes_probability_threshold(bins, probability_threshold, threshold):
    return any(matrix_area(b, threshold=threshold) >= probability_threshold for b in bins)


class ViewMatrix(ABC):
    args = None
    group_names = None
    experiment_names = None
    gene_ids = None

    class _ViewMatrix:
        @property
        def junctions(self):
            return self.get('junctions')

    @property
    def view_metadata(self):
        group_names = self.group_names
        experiment_names = self.experiment_names
        metadata = {'group_names': group_names}

        if experiment_names.size > 1:
            metadata['experiment_names'] = [np.insert(exps, 0, '{0} Combined'.format(group)) for exps, group in
                                            zip(experiment_names, group_names)]
        else:
            metadata['experiment_names'] = experiment_names

        return metadata

    @abstractmethod
    def lsv_ids(self, gene_ids=None):
        return ()

    @abstractmethod
    def valid_lsvs(self, lsv_ids):
        return ()

    def view_lsv_ids(self):
        args = self.args
        if args.lsv_ids:
            lsv_ids = args.lsv_ids
        elif args.gene_ids:
            lsv_ids = self.lsv_ids(args.gene_ids)
        else:
            lsv_ids = self.lsv_ids()

        yield from self.valid_lsvs(lsv_ids)

    def view_lsv_count(self):
        value = len(tuple(self.view_lsv_ids()))
        voila_log().info('Found {} LSVs'.format(value))
        if not value:
            raise NoLsvsFound()
        return value

    def page_count(self):
        gene_count = len(tuple(self.view_gene_ids()))
        return math.ceil(gene_count / constants.MAX_GENES)

    def paginated_genes(self):
        def grouper(iterable, n, fillvalue=None):
            args = [iter(iterable)] * n
            return zip_longest(*args, fillvalue=fillvalue)

        for page in grouper(self.view_gene_ids(), constants.MAX_GENES):
            yield tuple(p for p in page if p is not None)

    def view_gene_ids(self):
        args = self.args
        if args.gene_ids:
            gene_ids = args.gene_ids
        else:
            gene_ids = self.gene_ids

        for gene_id in gene_ids:
            if any(self.valid_lsvs(self.view_gene_lsvs(gene_id))):
                yield gene_id

    def view_gene_lsvs(self, gene_id):
        yield from self.valid_lsvs(self.lsv_ids([gene_id]))


class ViewPsi(Psi, ViewMatrix):
    def __init__(self, args):
        super().__init__(args.voila_file)
        self.args = args

    class _ViewPsi(Psi._Psi, ViewMatrix._ViewMatrix):
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
        def means(self):
            yield from unpack_means(self.get('means'))

        @property
        def group_means(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.means

        @property
        def group_means_rounded(self):
            for group_name, means in self.group_means:
                yield group_name, np.around(tuple(means), decimals=3)

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
                projection_prod = bins * np.arange(step / 2, 1, step)
                return np.sum(projection_prod)

            for b in self.bins:
                epsi = get_expected_psi(b)
                step_bins = 1.0 / b.size
                projection_prod = b * np.arange(step_bins / 2, 1, step_bins) ** 2
                yield np.sum(projection_prod) - epsi ** 2

        @property
        def junction_count(self):
            return len(tuple(self.means))

    def psi(self, lsv_id):
        return self._ViewPsi(self, lsv_id)

    def valid_lsvs(self, lsv_ids):
        args = self.args

        for lsv_id in lsv_ids:
            if not args.lsv_ids or lsv_id in args.lsv_ids:
                yield lsv_id


class ViewDeltaPsi(DeltaPsi, ViewMatrix):
    def __init__(self, args):
        super().__init__(args.voila_file)
        self.args = args

    class _ViewDeltaPsi(DeltaPsi._DeltaPsi, ViewMatrix._ViewMatrix):
        def __init__(self, matrix_hdf5, lsv_id):
            super().__init__(matrix_hdf5, lsv_id)

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
            return np.around(tuple(self.means), decimals=3)

        @property
        def group_means(self):
            for value in self.get('group_means'):
                yield unpack_means(value)

        @property
        def group_means_rounded(self):
            group_names = self.matrix_hdf5.group_names
            for group_name, means in zip(group_names, self.group_means):
                yield group_name, np.around(means, decimals=3)

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

        def is_lsv_changing(self, threshold):
            return is_lsv_changing(self.means, threshold)

        def probability_threshold(self):
            args = self.matrix_hdf5.args
            probability_threshold = args.probability_threshold
            threshold = args.probability_threshold
            return any(matrix_area(b, threshold=threshold) >= probability_threshold for b in self.bins)

        def high_probability_non_changing(self):
            prior = np.log(self.matrix_hdf5.prior[1 if self.intron_retention else 0])
            non_changing_threshold = self.matrix_hdf5.args.non_changing_threshold
            for bin in self.bins:
                A = np.log(bin) - prior
                R = np.exp(A - scipy.special.logsumexp(A))
                yield matrix_area(R, non_changing_threshold, non_changing=True)

    def delta_psi(self, lsv_id):
        return self._ViewDeltaPsi(self, lsv_id)

    def valid_lsvs(self, lsv_ids):
        threshold = None
        probability_threshold = None
        args = self.args

        if not args.show_all:
            threshold = args.threshold
            probability_threshold = args.probability_threshold

        for lsv_id in lsv_ids:
            delta_psi = self.delta_psi(lsv_id)
            if not args.lsv_ids or lsv_id in args.lsv_ids:
                if not threshold or delta_psi.is_lsv_changing(threshold):
                    if not probability_threshold or delta_psi.probability_threshold():
                        yield lsv_id


class ViewHeterogen(Heterogen, ViewMatrix):
    def valid_lsvs(self, lsv_ids):
        args = self.args

        for lsv_id in lsv_ids:
            if not args.lsv_ids or lsv_id in args.lsv_ids:
                yield lsv_id

    def __init__(self, args):
        super().__init__(args.voila_file)
        self.args = args

    class _ViewHeterogen(Heterogen._Heterogen, ViewMatrix._ViewMatrix):
        def __init__(self, matrix_hdf5, lsv_id):
            super().__init__(matrix_hdf5, lsv_id)

        @property
        def mean_psi(self):
            return self.get('mean_psi')

        @property
        def mu_psi(self):
            return self.get('mu_psi')

        @property
        def junction_stats(self):
            return self.get('junction_stats')

    def heterogen(self, lsv_id):
        return self._ViewHeterogen(self, lsv_id)

    @property
    def view_metadata(self):
        return {**super().view_metadata, **{
            'stat_names': self.stat_names
        }}
