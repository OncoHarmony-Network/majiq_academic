import operator
from abc import ABC
from functools import reduce
from itertools import chain

import numpy as np
import scipy.special

from voila.api.matrix_hdf5 import DeltaPsi, Psi, Heterogen
from voila.config import Config
from voila.exceptions import LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile
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
    group_names = None
    experiment_names = None
    gene_ids = None

    class _ViewMatrix:
        @property
        def junctions(self):
            return self.get('junctions')

    # @property
    # def metadata(self):
    #     group_names = self.group_names.tolist()
    #     experiment_names = self.experiment_names.tolist()
    #     metadata = {'group_names': group_names}
    #
    #     if len(experiment_names) > 1:
    #         metadata['experiment_names'] = [[group + ' Combined'] + exps for exps, group in
    #                                         zip(experiment_names, group_names)]
    #     else:
    #         metadata['experiment_names'] = experiment_names
    #
    #     return metadata


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

        @property
        def group_means(self):
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], list(self.means)

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


class ViewHeterogens:
    def __init__(self):
        self._group_names = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    class _ViewHeterogens:
        def __init__(self, matrix_hdf5, lsv_id):
            self.matrix_hdf5 = matrix_hdf5
            self.lsv_id = lsv_id

        def get_all(self):
            yield '_id', self.lsv_id
            yield 'gene_id', self.gene_id
            yield 'mean_psi', tuple(self.mean_psi)
            yield 'reference_exon', self.reference_exon
            yield 'lsv_type', self.lsv_type
            yield 'dpsi', self.dpsi
            for stat_name in self.matrix_hdf5.stat_names:
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
            voila_files = Config().voila_files
            s = set()
            for f in voila_files:
                with ViewHeterogen(f) as m:
                    try:
                        het = m.heterogen(self.lsv_id)
                        s.add(getattr(het, attr))
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
        def reference_exon(self):
            ref_ex = {tuple(h.reference_exon) for h in self.heterogens}
            assert len(ref_ex) == 1
            return ref_ex.pop()

        @property
        def dpsi(self):
            d = {}

            for het in self.heterogens:
                try:
                    group_names = het.matrix_hdf5.group_names
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
                    group_names = het.matrix_hdf5.group_names
                    stat_names = het.matrix_hdf5.stat_names
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

        @property
        def mu_psi(self):
            voila_files = Config().voila_files
            group_names = self.matrix_hdf5.group_names
            experiment_names = self.matrix_hdf5.experiment_names
            exps_len = max(len(e) for e in experiment_names)
            juncs_len = len(self.junctions)
            grps_len = len(group_names)
            mu_psi = np.empty((grps_len, juncs_len, exps_len))

            mu_psi.fill(-1)

            for f in voila_files:
                with ViewHeterogen(f) as m:
                    try:

                        het = m.lsv(self.lsv_id)

                        mus = het.mu_psi
                        mus = mus.transpose((1, 0, 2))

                        for mu, grp in zip(mus, m.group_names):
                            idx = group_names.index(grp)
                            mu_shp = mu.shape
                            mu_psi[idx][0:mu_shp[0], 0:mu_shp[1]] = mu

                    except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                        pass

            mu_psi = mu_psi.transpose((1, 0, 2))
            return mu_psi.tolist()

        @property
        def mean_psi(self):
            voila_files = Config().voila_files
            group_names = self.matrix_hdf5.group_names
            juncs_len = len(self.junctions)
            grps_len = len(group_names)
            mean_psi = np.empty((grps_len, juncs_len, 40))
            mean_psi.fill(-1)

            for f in voila_files:
                with ViewHeterogen(f) as m:
                    try:
                        het = m.lsv(self.lsv_id)

                        means = het.mean_psi
                        means = means.transpose((1, 0, 2))

                        for mn, grp in zip(means, m.group_names):
                            idx = group_names.index(grp)
                            mn_shp = mn.shape
                            mean_psi[idx][0:mn_shp[0], 0:mn_shp[1]] = mn
                    except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                        pass

            mean_psi = mean_psi.transpose((1, 0, 2))
            return mean_psi.tolist()

        @property
        def lsv_type(self):
            return self.get_attr('lsv_type')

        @property
        def intron_retention(self):
            return self.get_attr('intron_retention')

        @property
        def junctions(self):
            config = Config()
            juncs = None
            for f in config.voila_files:
                with ViewHeterogen(f) as m:
                    if juncs is None:
                        try:
                            juncs = m.lsv(self.lsv_id).junctions
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                            pass
                    if juncs is not None:
                        return juncs

        @property
        def junction_stats(self):
            config = Config()
            voila_files = config.voila_files
            for f in voila_files:
                with ViewHeterogen(f) as m:
                    het = m.lsv(self.lsv_id)
                    groups = '_'.join(m.group_names)
                    stat_names = m.stat_names
                    try:
                        for name, stat in zip(stat_names, het.junction_stats.T):
                            if len(voila_files) == 1:
                                yield name, stat
                            else:
                                yield groups + ' ' + name, stat
                    except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                        pass

    @property
    def junction_stats_column_names(self):
        voila_files = Config().voila_files

        for f in voila_files:
            with ViewHeterogen(f) as m:
                groups = '_'.join(m.group_names)
                for name in m.stat_names:
                    if len(voila_files) == 1:
                        yield name
                    else:
                        yield groups + ' ' + name

    @property
    def analysis_type(self):
        analysis_types = {vh.analysis_type for vh in self.view_heterogens}
        if len(analysis_types) == 1:
            return analysis_types.pop()

    @property
    def group_names(self):
        if self._group_names is None:
            config = Config()
            group_names = set()
            for f in config.voila_files:
                with ViewHeterogen(f) as m:
                    for grp in m.group_names:
                        group_names.add(grp)

            self._group_names = list(sorted(group_names))

        return self._group_names

    @property
    def experiment_names(self):
        config = Config()
        exp_names = {}
        for f in config.voila_files:
            with ViewHeterogen(f) as m:
                for exp, grp in zip(m.experiment_names, m.group_names):
                    exp_names[grp] = exp

        return [exp_names[grp] for grp in self.group_names]

    def lsv(self, lsv_id):
        return self._ViewHeterogens(self, lsv_id)

    @property
    def gene_ids(self):
        voila_files = Config().voila_files
        vhs = [ViewHeterogen(f) for f in voila_files]
        yield from set(chain(*(v.gene_ids for v in vhs)))
        for v in vhs:
            v.close()

    def lsv_ids(self, gene_ids=None):
        voila_files = Config().voila_files
        vhs = [ViewHeterogen(f) for f in voila_files]
        yield from set(chain(*(v.lsv_ids(gene_ids) for v in vhs)))
        for v in vhs:
            v.close()


class ViewHeterogen(Heterogen, ViewMatrix):
    def __init__(self, voila_file):
        try:
            super().__init__(voila_file)
        except OSError:
            print(voila_file)
            raise

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

    def lsv(self, lsv_id):
        return self._ViewHeterogen(self, lsv_id)
