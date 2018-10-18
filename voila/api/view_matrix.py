import math
import operator
from abc import ABC, abstractmethod
from functools import reduce
from itertools import zip_longest

import numpy as np
import scipy.special

from voila import constants
from voila.api.matrix_hdf5 import DeltaPsi, Psi, Heterogen
from voila.config import Config
from voila.exceptions import NoLsvsFound, LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile
from voila.utils.voila_log import voila_log
from voila.vlsv import get_expected_dpsi, is_lsv_changing, matrix_area, get_expected_psi


def unpack_means(value):
    if np.size(value, 0) == 1:
        value = np.append(value, np.array(1 - value[0]))
    return value.tolist()


def unpack_bins(value):
    if np.size(value, 0) == 1:
        value = np.append(value, [np.flip(value[-1], 0)], axis=0)
    return value.tolist()


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
    def metadata(self):
        group_names = self.group_names.tolist()
        experiment_names = self.experiment_names.tolist()
        metadata = {'group_names': group_names}

        if len(experiment_names) > 1:
            metadata['experiment_names'] = [[group + ' Combined'] + exps for exps, group in
                                            zip(experiment_names, group_names)]
        else:
            metadata['experiment_names'] = experiment_names

        return metadata

    @abstractmethod
    def lsv_ids(self, gene_ids=None):
        return ()

    @abstractmethod
    def valid_lsvs(self, lsv_ids):
        # todo: There should be a way to filter lsv ids for tsvs
        for lsv_id in lsv_ids:
            yield lsv_id

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
        for gene_id in self.gene_ids:
            if any(self.valid_lsvs(self.view_gene_lsvs(gene_id))):
                yield gene_id

    def view_gene_lsvs(self, gene_id):
        yield from self.valid_lsvs(self.lsv_ids([gene_id]))


class ViewPsi(Psi, ViewMatrix):
    def __init__(self):
        config = Config()
        super().__init__(config.voila_file)

    class _ViewPsi(Psi._Psi, ViewMatrix._ViewMatrix):
        def get_all(self):
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'reference_exon', tuple(self.reference_exon)
            yield 'gene_id', self.gene_id
            yield 'target', self.target

            fields = list(self.fields)

            fields.remove('bins')
            yield 'group_bins', dict(self.group_bins)

            fields.remove('means')
            yield 'group_means', dict(self.group_means)

            yield from self.get_many(fields)

        @property
        def means(self):
            yield from unpack_means(self.get('means'))

        # @property
        # def junctions(self):
        #     return self.get('junctions')

        @property
        def group_means(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], list(self.means)

        # @property
        # def group_means(self):
        #     for group_name, means in self.group_means:
        #         yield group_name, tuple(means)

        @property
        def bins(self):
            return unpack_bins(self.get('bins'))

        @property
        def group_bins(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.bins

        @property
        def variances(self):
            for b in self.bins:
                epsi = get_expected_psi(b)
                # b used to be a numpy array and now it's a list...
                step_bins = 1.0 / len(b)
                projection_prod = b * np.arange(step_bins / 2, 1, step_bins) ** 2
                yield np.sum(projection_prod) - epsi ** 2

        @property
        def junction_count(self):
            return len(tuple(self.means))

    def lsv(self, lsv_id):
        return self._ViewPsi(self, lsv_id)

    def valid_lsvs(self, lsv_ids):
        yield from super().valid_lsvs(lsv_ids)


class ViewDeltaPsi(DeltaPsi, ViewMatrix):
    def __init__(self):
        self.config = Config()
        super().__init__(self.config.voila_file)

    class _ViewDeltaPsi(DeltaPsi._DeltaPsi, ViewMatrix._ViewMatrix):
        def __init__(self, matrix_hdf5, lsv_id):
            self.config = matrix_hdf5.config
            super().__init__(matrix_hdf5, lsv_id)

        def get_all(self):
            yield 'gene_id', self.gene_id
            yield 'lsv_id', self.lsv_id
            yield '_id', self.lsv_id
            yield 'reference_exon', self.reference_exon
            yield 'excl_incl', self.excl_incl
            yield 'target', self.target

            fields = list(self.fields)

            fields.remove('group_bins')
            yield 'group_bins', dict(self.group_bins)

            fields.remove('group_means')
            yield 'group_means', dict(self.group_means)

            fields.remove('bins')
            yield 'bins', self.bins
            yield 'means', self.means

            yield from self.get_many(fields)

        @property
        def bins(self):
            return unpack_bins(self.get('bins'))

        @property
        def group_bins(self):
            group_names = self.matrix_hdf5.group_names
            for group_name, value in zip(group_names, self.get('group_bins')):
                yield group_name, unpack_bins(value)

        @property
        def means(self):
            for b in self.bins:
                yield get_expected_dpsi(b)

        @property
        def group_means(self):
            group_names = self.matrix_hdf5.group_names
            for group_name, means in zip(group_names, self.get('group_means')):
                yield group_name, means.tolist()

        @property
        def excl_incl(self):
            for mean in self.means:
                if mean < 0:
                    yield (-mean, 0)
                else:
                    yield (0, mean)

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
            prior = self.matrix_hdf5.prior[1 if self.intron_retention else 0]
            non_changing_threshold = self.config.non_changing_threshold
            for bin in self.bins:
                A = np.log(bin) - prior
                R = np.exp(A - scipy.special.logsumexp(A))
                yield matrix_area(R, non_changing_threshold, non_changing=True)

    def lsv(self, lsv_id):
        return self._ViewDeltaPsi(self, lsv_id)

    def valid_lsvs(self, lsv_ids):
        threshold = None
        probability_threshold = None
        config = self.config

        if not config.show_all:
            threshold = config.threshold
            probability_threshold = config.probability_threshold

        for lsv_id in super().valid_lsvs(lsv_ids):
            delta_psi = self.lsv(lsv_id)
            if not threshold or delta_psi.is_lsv_changing(threshold):
                if not probability_threshold or delta_psi.probability_threshold():
                    yield lsv_id


class ViewHeterogens:
    def __init__(self, args):
        config = Config()
        self.args = args
        self.view_heterogens = tuple(ViewHeterogen(args, f) for f in config.voila_files)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    class _ViewHeterogens:
        def __init__(self, matrix_hdf5, lsv_id):
            self.matrix_hdf5 = matrix_hdf5
            self.lsv_id = lsv_id
            self.heterogens = tuple(vh.heterogen(lsv_id) for vh in self.matrix_hdf5.view_heterogens)

        def get_all(self):
            yield '_id', self.lsv_id
            yield 'gene_id', self.gene_id
            yield 'mean_psi', tuple(self.mean_psi)
            yield 'reference_exon', self.reference_exon
            yield 'lsv_type', self.lsv_type
            yield 'dpsi', self.dpsi
            print(self.matrix_hdf5.view_metadata['stat_names'])
            for stat_name in self.matrix_hdf5.view_metadata['stat_names']:
                yield stat_name, self.junction_stat(stat_name)
            yield 'mu_psi', self.mu_psi
            yield 'junctions', self.junctions
            yield 'A5SS', self.a5ss
            yield 'A3SS', self.a3ss
            yield 'exon_skipping', self.exon_skipping
            yield 'exon_count', self.exon_count
            yield 'target', self.target
            yield 'binary', self.binary

        def get_attr(self, attr):
            s = set()
            for h in self.heterogens:
                try:
                    s.add(getattr(h, attr))
                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    pass
            assert len(s) == 1
            return s.pop()

        @property
        def gene_id(self):
            return self.get_attr('gene_id')

        @property
        def a5ss(self):
            return self.get_attr('a5ss')

        @property
        def a3ss(self):
            return self.get_attr('a3ss')

        @property
        def exon_skipping(self):
            return self.get_attr('exon_skipping')

        @property
        def exon_count(self):
            return self.get_attr('exon_count')

        @property
        def target(self):
            return self.get_attr('target')

        @property
        def binary(self):
            return self.get_attr('binary')

        @property
        def dpsi(self):
            d = {}

            for het in self.heterogens:
                try:
                    group_names = het.matrix_hdf5.view_metadata['group_names']
                    dpsi_values = list(het.dpsi)
                    try:
                        d[group_names[0]][group_names[1]] = dpsi_values
                    except KeyError:
                        d[group_names[0]] = {group_names[1]: dpsi_values}
                except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                    pass

            return d

        def junction_stat(self, stat_name):
            d = {}
            for het in self.heterogens:
                try:
                    group_names = het.matrix_hdf5.view_metadata['group_names']
                    stat_names = het.matrix_hdf5.view_metadata['stat_names']
                    stat_idx = list(stat_names).index(stat_name)
                    trans_junc_stats = het.junction_stats.T
                    stats_values = trans_junc_stats[stat_idx].tolist()
                    try:
                        d[group_names[0]][group_names[1]] = stats_values
                    except KeyError:
                        d[group_names[0]] = {group_names[1]: stats_values}
                except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                    pass
            return d

            #
            # group_names = self.matrix_hdf5.view_metadata['group_names']
            # ir = int(self.lsv_type.endswith('i'))
            # arr = np.empty((len(self.junctions) + ir, len(group_names), len(group_names)), dtype=np.float)
            # arr.fill(-1)
            #
            # for het in self.heterogens:
            #     try:
            #         stat_names = het.matrix_hdf5.view_metadata['stat_names']
            #         stat_idx = list(stat_names).index(stat_name)
            #         trans_junc_stats = het.junction_stats.T
            #         for junc_idx, stat in enumerate(trans_junc_stats[stat_idx]):
            #             z, y = sorted(group_names.index(g) for g in het.matrix_hdf5.view_metadata['group_names'])
            #             assert arr[junc_idx][y][z] == -1
            #             arr[junc_idx][y][z] = stat
            #
            #     except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
            #         pass
            #
            # return self.compact_array(arr)

        @staticmethod
        def compact_array(arr):
            mask = np.ones(arr.shape[1:], dtype=bool)
            triangle_mask = np.tril(mask, k=-1)
            return arr[:, triangle_mask]

        @property
        def reference_exon(self):
            ref_ex = {tuple(h.reference_exon) for h in self.heterogens}
            assert len(ref_ex) == 1
            return ref_ex.pop()

        def psi_data(self, data_type):
            d = {}
            for het in self.heterogens:
                try:
                    for psi_d, sample in zip(getattr(het, data_type), het.matrix_hdf5.view_metadata['group_names']):
                        if sample in d:
                            test_shape = min(d[sample].shape, psi_d.shape)
                            arr1 = d[sample][0:test_shape[0], 0:test_shape[1]]
                            arr2 = psi_d[0:test_shape[0], 0:test_shape[1]]
                            arr_d = arr2 - arr1
                            threshold = 0.001
                            if (np.abs(arr_d) > threshold).any():
                                voila_log().warning('PSI distribution from {} containing the group "{}" '
                                                    'don\'t agree within a {} threshold.  Double check the '
                                                    'files containing this group. {}.'.format(
                                    het.matrix_hdf5.h.filename, sample, threshold, data_type))
                        else:
                            d[sample] = psi_d
                except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                    pass
                except IndexError:
                    print(data_type)
                    print(psi_d)
                    exit()

            for sample in self.matrix_hdf5.view_metadata['group_names']:
                try:
                    yield d[sample]
                except KeyError:
                    yield None

        @property
        def mu_psi(self):
            meta = self.matrix_hdf5.view_metadata
            has_combined = meta['experiment_names'][0][0].endswith('Combined')
            exp_count = max(len(x) for x in meta['experiment_names']) - int(has_combined)
            grp_count = len(self.matrix_hdf5.view_metadata['group_names'])
            junc_count = len(self.junctions) + int(self.lsv_type.endswith('i'))
            arr = np.empty((grp_count, exp_count, junc_count))

            arr.fill(-1)

            for idx, x in enumerate(self.psi_data('mu_psi')):
                if x is not None:
                    arr[idx][0:x.shape[0], 0:x.shape[1]] = x

            for xs in arr.T:
                yield [[] if np.all(x == -1) else x for x in xs.T]

        @property
        def mean_psi(self):
            yield from self.psi_data('mean_psi')

        @property
        def lsv_type(self):
            return self.get_attr('lsv_type')

        @property
        def intron_retention(self):
            return self.get_attr('intron_retention')

        @property
        def junctions(self):
            juncs = None
            source_file = None

            for h in self.heterogens:
                try:
                    if juncs is None:
                        juncs = h.junctions
                        source_file = h.matrix_hdf5.h.filename

                    if not np.array_equal(juncs, h.junctions):
                        filename = h.matrix_hdf5.h.filename
                        voila_log().warning('For junctions in LSV {}, {} and {} do not agree.'.format(self.lsv_id,
                                                                                                      filename,
                                                                                                      source_file))
                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    pass

            return juncs

        @property
        def junction_stats(self):
            hets = self.heterogens
            one_het = len(hets) == 1
            for het in hets:
                try:
                    meta = het.matrix_hdf5.view_metadata
                    groups = '_'.join(meta['group_names'])
                    for name, stat in zip(meta['stat_names'], het.junction_stats.T):
                        if one_het:
                            yield name, stat
                        else:
                            yield '{} {}'.format(groups, name), stat
                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    pass

    @property
    def one_heterogen(self):
        return len(self.view_heterogens) == 1

    @property
    def junction_stats_column_names(self):
        one_het = self.one_heterogen
        for vh in self.view_heterogens:
            meta = vh.view_metadata
            groups = '_'.join(meta['group_names'])
            for name in meta['stat_names']:
                if one_het:
                    yield name
                else:
                    yield '{} {}'.format(groups, name)

    @property
    def analysis_type(self):
        analysis_types = {vh.analysis_type for vh in self.view_heterogens}
        if len(analysis_types) == 1:
            return analysis_types.pop()

    @property
    def view_metadata(self):
        group_names = tuple(sorted({a for vh in self.view_heterogens for a in vh.view_metadata['group_names']}))
        stat_names = tuple(sorted({a for vh in self.view_heterogens for a in vh.view_metadata['stat_names']}))
        experiment_names = {}

        for vh in self.view_heterogens:
            for group, exp in zip(vh.view_metadata['group_names'], vh.view_metadata['experiment_names']):
                exp = sorted(e for e in exp if e)
                if group in experiment_names:
                    assert np.array_equal(experiment_names[group], exp), '{} {}'.format(experiment_names[group], exp)
                else:
                    experiment_names[group] = exp

        experiment_names = tuple(experiment_names[group] for group in group_names)

        return {
            'experiment_names': experiment_names,
            'group_names': group_names,
            'stat_names': stat_names,
            '_id': 'metadata'
        }

    def view_gene_ids(self):
        yield from sorted({g for vh in self.view_heterogens for g in vh.view_gene_ids()})

    def view_gene_lsvs(self, gene_id):
        yield from {l for vh in self.view_heterogens for l in vh.view_gene_lsvs(gene_id)}

    def heterogen(self, lsv_id):
        return self._ViewHeterogens(self, lsv_id)

    def paginated_genes(self):
        def grouper(iterable, n, fillvalue=None):
            args = [iter(iterable)] * n
            return zip_longest(*args, fillvalue=fillvalue)

        for page in grouper(self.view_gene_ids(), constants.MAX_GENES):
            yield tuple(p for p in page if p is not None)

    def page_count(self):
        gene_count = len(tuple(self.view_gene_ids()))
        return math.ceil(gene_count / constants.MAX_GENES)

    @property
    def gene_ids(self):
        yield from {g for vh in self.view_heterogens for g in vh.gene_ids}

    def view_lsv_ids(self):
        yield from set(lsv_id for vh in self.view_heterogens for lsv_id in vh.view_lsv_ids())


class ViewHeterogen(Heterogen, ViewMatrix):
    def __init__(self, args, voila_file):
        super().__init__(voila_file)
        self.args = args

    def valid_lsvs(self, lsv_ids):
        yield from super().valid_lsvs(lsv_ids)

    class _ViewHeterogen(Heterogen._Heterogen, ViewMatrix._ViewMatrix):
        def __init__(self, matrix_hdf5, lsv_id):
            super().__init__(matrix_hdf5, lsv_id)

        def get_all(self):
            yield '_id', self.lsv_id
            yield 'mean_psi', self.mean_psi
            yield 'dpsi', list(self.dpsi)
            for idx, stat_name in enumerate(self.matrix_hdf5.stat_names):
                yield stat_name, self.junction_stat(idx)

        def junction_stat(self, idx):
            return self.junction_stats.T[idx]

        @property
        def dpsi(self):
            for bins in zip(*self.mean_psi):
                yield abs(reduce(operator.__sub__, (get_expected_psi(b) for b in bins)))

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
        return {**super().view_metadata, **{'stat_names': self.stat_names}}
