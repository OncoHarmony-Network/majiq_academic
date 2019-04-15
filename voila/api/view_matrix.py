import operator
from abc import ABC
from functools import reduce
from itertools import chain

import numpy as np

from voila.api.matrix_hdf5 import DeltaPsi, Psi, Heterogen, MatrixType
from voila.api.matrix_utils import unpack_means, unpack_bins, generate_excl_incl, generate_means, \
    generate_high_probability_non_changing, generate_variances
from voila.config import ViewConfig
from voila.exceptions import LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile, LsvIdNotFoundInAnyVoilaFile
from voila.vlsv import is_lsv_changing, matrix_area, get_expected_psi


class ViewMatrix(ABC):
    group_names = None
    experiment_names = None
    gene_ids = None

    def lsv_ids(self, gene_ids=None):
        raise NotImplementedError()

    def lsv(self, lsv_id):
        raise NotImplementedError()

    def lsvs(self, gene_id=None):
        """
        Get a generator for all the lsvs.  If gene id is specified, then only return lsv for that gene.
        :param gene_id: gene id
        :return: generator of lsv objects
        """
        if gene_id:
            for lsv_id in self.lsv_ids(gene_ids=[gene_id]):
                yield self.lsv(lsv_id)
        else:
            for lsv_id in self.lsv_ids():
                yield self.lsv(lsv_id)


class ViewMatrixType(MatrixType):
    def __init__(self, matrix_hdf5, lsv_id, fields):
        super().__init__(matrix_hdf5, lsv_id, fields)

    @property
    def means(self):
        raise NotImplementedError()

    @property
    def bins(self):
        """
        Get bins data from voila file.
        :return: list
        """
        return unpack_bins(self.get('bins'))

    @property
    def junction_count(self):
        """
        Get count of all junctions in this lsv.
        :return: integer
        """
        return len(tuple(self.means))


class ViewMulti:
    """
    View for set of Voila  files.  This is used in creation of tsv and html files.
    Base class
    """
    def __init__(self, view_class):
        """
        :param view_class: class for view of single item (ex, ViewPsi, ViewHeterogen)
        """
        self._group_names = None
        self.view_class = view_class

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    @property
    def experiment_names(self):
        """
        Experiment names for this set of het voila files.
        :return: List
        """
        config = ViewConfig()
        exp_names = {}
        for f in config.voila_files:
            with self.view_class(f) as m:
                for exp, grp in zip(m.experiment_names, m.group_names):
                    exp_names[grp] = exp

        return [exp_names[grp] for grp in self.group_names]

    @property
    def group_names(self):
        """
        Group names for this set of het voila files.
        :return: list
        """
        config = ViewConfig()
        grp_names = []
        for f in config.voila_files:
            with self.view_class(f) as m:
                for grp in m.group_names:
                    if not grp in grp_names:
                        grp_names.append(grp)

        return grp_names

    @property
    def splice_graph_experiment_names(self):
        """
        experiment names for the splice graph drop down.
        :return: List
        """

        config = ViewConfig()
        exp_names = {}
        for f in config.voila_files:
            with self.view_class(f) as m:
                for exp, grp in zip(m.splice_graph_experiment_names, m.group_names):
                    exp_names[grp] = exp

        return [exp_names[grp] for grp in self.group_names]

    @property
    def gene_ids(self):
        """
        Get a set of gene ids from all het voila files.
        :return: generator
        """

        voila_files = ViewConfig().voila_files
        vhs = [self.view_class(f) for f in voila_files]
        yield from set(chain(*(v.gene_ids for v in vhs)))
        for v in vhs:
            v.close()

    def lsv_ids(self, gene_ids=None):
        """
        Get a set of lsv ids from all voila files for specified gene ids. If gene ids is None, then get all lsv ids.
        :param gene_ids: list of gene ids
        :return:
        """

        voila_files = ViewConfig().voila_files
        vhs = [self.view_class(f) for f in voila_files]
        yield from set(chain(*(v.lsv_ids(gene_ids) for v in vhs)))
        for v in vhs:
            v.close()

    def lsv(self, lsv_id):
        raise NotImplementedError()

    def lsvs(self, gene_id=None):
        """
        Get all lsvs for set of het voila files.
        :param gene_id: gene id
        :return:
        """

        if gene_id:
            gene_ids = [gene_id]
        else:
            gene_ids = None

        for lsv_id in self.lsv_ids(gene_ids):
            yield self.lsv(lsv_id)



class _ViewMulti:
    """
    Base class for .lsv() call return of other multi objects
    """
    def __init__(self, matrix_hdf5, lsv_id, view_class):
        """
        :param view_class: class for view of single item (ex, ViewPsi, ViewHeterogen)
        """
        self.matrix_hdf5 = matrix_hdf5
        self.lsv_id = lsv_id
        self.view_class = view_class

    def _get_prop(self, prop, cast=None):
        """
        Look for the first input voila file with the property which exists. This does NOT validate the property
        is the same across all files where lsv id exists.

        cast should be specified when trying to retrieve generator style properties, because the generator
        will be processed inside this function, and attempting to iterate the generator outside of
        will cause an unrelated error because you have exited the "with" block.
        :return: property value
        """
        config = ViewConfig()
        propval = None
        for f in config.voila_files:
            with ViewPsi(f) as m:
                if propval is None:
                    try:

                        propval = getattr(m.lsv(self.lsv_id), prop)
                        if cast:
                            propval = cast(propval)

                    except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                        pass
                if propval is not None:
                    return propval

    def _get_prop_multi(self, prop):
        """
        Look for some attribute in all input files, and return a group:property dict
        it is assummed that the property itself will be read as a group:property dict
        from each individual file.
        :return: property value
        """
        config = ViewConfig()
        propval = None
        groups_to_props = {}
        for f in config.voila_files:
            with ViewPsi(f) as m:
                try:
                    propval = dict(getattr(m.lsv(self.lsv_id), prop))

                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    pass
                if propval is not None:
                    for key in propval:
                        groups_to_props[key] = propval[key]

        return groups_to_props


    def get_attr(self, attr):
        """
        For attributes that exist is each het file, this will get all values and confirm they're all equal.
        :param attr: attribute found in het voila file.
        :return: attribute value
        """

        voila_files = ViewConfig().voila_files
        s = set()
        for f in voila_files:
            with self.view_class(f) as m:
                try:
                    inner = m.lsv(self.lsv_id)
                    s.add(getattr(inner, attr))
                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    pass
        if len(s) == 0:
            raise LsvIdNotFoundInAnyVoilaFile
        assert len(s) == 1, s
        return s.pop()

    @property
    def reference_exon(self):
        return self.get_attr('reference_exon')

    @property
    def target(self):
        return self.get_attr('target')

    @property
    def source(self):
        return self.get_attr('source')

    @property
    def binary(self):
        return self.get_attr('binary')

    @property
    def complex(self):
        return self.get_attr('complex')

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
    def lsv_type(self):
        return self.get_attr('lsv_type')

    @property
    def intron_retention(self):
        return self.get_attr('intron_retention')


    @property
    def junctions(self):
        """
        Finds first file with junctions list for this specific lsv id. This does NOT validate the junctions list
        is the same across all files where lsv id exists.
        :return: numpy array
        """
        return self._get_prop('junctions')

class ViewPsis(ViewMulti):

    def __init__(self):
        super().__init__(ViewPsi)

    class _ViewPsis(_ViewMulti):
        def __init__(self, matrix_hdf5, lsv_id):
            super().__init__(matrix_hdf5, lsv_id, ViewPsi)

        @property
        def all_group_means(self):
            """
            Get group means from all files as a dict of group_name:means
            """
            return self._get_prop_multi('group_means')

        @property
        def group_means(self):
            """
            Finds first file with this specific lsv id and gets group means. This does NOT consider any other
            input files upon finding the first one with the lsv id
            """
            return self._get_prop('group_means', dict)

        @property
        def group_bins(self):
            """
            Finds first file with this specific lsv id and gets group means. This does NOT consider any other
            input files upon finding the first one with the lsv id
            """
            return self._get_prop('group_bins', dict)

    def lsv(self, lsv_id):
        """
        Get view heterogens object for this lsv id.
        :param lsv_id: lsv id
        :return: view heterogens object
        """
        return self._ViewPsis(self, lsv_id)





class ViewPsi(Psi, ViewMatrix):
    def __init__(self, voila_file=None):
        """
        This represents a single Psi voila file.  ViewPsis uses this class to retrieve data from the individual
        files.
        :param voila_file: voila file name
        """
        if not voila_file:
            voila_file = ViewConfig().voila_file
        super().__init__(voila_file)

    class _ViewPsi(Psi._Psi, ViewMatrixType):
        @property
        def means(self):
            """
            Get means data from voila file.
            :return: list
            """
            return list(unpack_means(self.get('means')))

        @property
        def group_bins(self):
            """
            Get bins in a dictionary where the key in the name of the group it belongs to.
            :return: generator of key, value
            """
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], self.bins

        @property
        def group_means(self):
            """
            Get means data from voila file.
            :return: generator
            """
            group_names = self.matrix_hdf5.group_names
            yield group_names[0], list(self.means)

        @property
        def variances(self):
            """
            Create variance data of bins data.
            :return: list
            """
            return generate_variances(self.bins)

    def lsv(self, lsv_id):
        """
        Get lsv object by lsv id.
        :param lsv_id: lsv id
        :return: lsv object
        """
        return self._ViewPsi(self, lsv_id)


class ViewDeltaPsi(DeltaPsi, ViewMatrix):
    def __init__(self):
        """
        View for delta psi matrix.  This is used in creation of tsv and html files.
        """
        self.config = ViewConfig()
        super().__init__(self.config.voila_file)

    class _ViewDeltaPsi(DeltaPsi._DeltaPsi, ViewMatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            self.config = matrix_hdf5.config
            super().__init__(matrix_hdf5, lsv_id)

        @property
        def group_bins(self):
            """
            Get dictionary of bins by group name.
            :return: generator of key, value
            """
            group_names = self.matrix_hdf5.group_names
            for group_name, value in zip(group_names, self.get('group_bins')):
                yield group_name, unpack_bins(value)

        @property
        def means(self):
            """
            Create mean data from bins data.
            :return: list
            """
            return generate_means(self.bins)

        @property
        def group_means(self):
            """
            Get dictionary of mean by group name.
            :return: generator of key, value
            """
            group_names = self.matrix_hdf5.group_names
            for group_name, means in zip(group_names, self.get('group_means')):
                yield group_name, means.tolist()

        @property
        def excl_incl(self):
            """
            Using means data, create exclude/include list.
            :return: list
            """
            return generate_excl_incl(self.means)

        def is_lsv_changing(self, threshold):
            """
            Is lsv changing based on threshold supplied by user.
            :param threshold: threshold value, usually 0.2
            :return: boolean
            """
            return is_lsv_changing(self.means, threshold)

        def probability_threshold(self):
            """
            Get probability that the junction in an LSV are changing.

            :return: list of p values
            """
            args = self.matrix_hdf5.args
            probability_threshold = args.probability_threshold
            threshold = args.probability_threshold
            return any(matrix_area(b, threshold=threshold) >= probability_threshold for b in self.bins)

        def high_probability_non_changing(self):
            """
            Get probability that junctions in an lsv aren't changing.
            :return: list
            """
            return generate_high_probability_non_changing(self.intron_retention, self.matrix_hdf5.prior,
                                                          self.config.non_changing_threshold, self.bins)

    def lsv(self, lsv_id):
        """
        Get delta psi object by lsv id.
        :param lsv_id: lsv id
        :return: delta psi object
        """
        return self._ViewDeltaPsi(self, lsv_id)






class ViewHeterogens(ViewMulti):

    def __init__(self):
        super().__init__(ViewHeterogen)

    class _ViewHeterogens(_ViewMulti):
        def __init__(self, matrix_hdf5, lsv_id):
            super().__init__(matrix_hdf5, lsv_id, ViewHeterogen)

        @property
        def group_bins(self):
            """
            Associate means values with it's experiment group name.
            :return: Generator key/value
            """

            mean_psi = self.mean_psi
            mean_psi = np.array(mean_psi)
            mean_psi = mean_psi.transpose((1, 0, 2))
            for group_name, mean in zip(self.matrix_hdf5.group_names, mean_psi):
                yield group_name, mean.tolist()

        def junction_heat_map(self, stat_name, junc_idx):
            voila_files = ViewConfig().voila_files
            hets_grps = self.matrix_hdf5.group_names
            hets_grps_len = len(hets_grps)
            s = np.ndarray((hets_grps_len, hets_grps_len))

            s.fill(-1)

            for f in voila_files:
                with ViewHeterogen(f) as m:

                    het = m.lsv(self.lsv_id)
                    try:

                        stat_idx = het.matrix_hdf5.stat_names
                        stat_idx = list(stat_idx)
                        stat_idx = stat_idx.index(stat_name)

                        stat_value = het.junction_stats
                        stat_value = stat_value.T
                        stat_value = stat_value[stat_idx][junc_idx]

                        dpsi_value = het.dpsi_signed
                        dpsi_value = dpsi_value[junc_idx]

                        grp_names = m.group_names
                        grp_names.sort()

                        x = hets_grps.index(grp_names[0])
                        y = hets_grps.index(grp_names[1])

                        s[x][y] = dpsi_value
                        s[y][x] = stat_value

                    except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                        pass

            return s.tolist()

        @property
        def mu_psi(self):
            """
            Find mu_psi in all het voila files and create a matrix that matches the new unified set of group/experiment
            names.  In the case where the experiments are all the same number, this will fill the empty space in the
            matrix with -1.
            :return: List
            """
            voila_files = ViewConfig().voila_files
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
            """
             Find mean_psi in all het voila files and create a matrix that matches the new unified set of group/experiment
            names.  In the case where the experiments are all the same number, this will fill the empty space in the
            matrix with -1.
            :return: List
            """
            voila_files = ViewConfig().voila_files
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
        def junction_stats(self):
            """
            This gets associates stat test names with their values.
            :return: generator key/value
            """
            config = ViewConfig()
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
    def stat_names(self):
        """
        Gets list of stat test names.
        :return: List
        """
        names = set()
        voila_files = ViewConfig().voila_files
        for f in voila_files:
            with ViewHeterogen(f) as m:
                for s in m.stat_names:
                    names.add(s)

        return list(sorted(names))

    @property
    def junction_stats_column_names(self):
        """
        Stat column names for tsv output.
        :return: generator
        """
        voila_files = ViewConfig().voila_files

        for f in voila_files:
            with ViewHeterogen(f) as m:
                groups = '_'.join(m.group_names)
                for name in m.stat_names:
                    if len(voila_files) == 1:
                        yield name
                    else:
                        yield groups + ' ' + name

    def lsv(self, lsv_id):
        """
        Get view heterogens object for this lsv id.
        :param lsv_id: lsv id
        :return: view heterogens object
        """

        return self._ViewHeterogens(self, lsv_id)





class ViewHeterogen(Heterogen, ViewMatrix):
    def __init__(self, voila_file):
        """
        This represents a single het voila file.  ViewHeterogens uses this class to retrieve data from the individual
        files.
        :param voila_file: voila file name
        """
        super().__init__(voila_file)

    class _ViewHeterogen(Heterogen._Heterogen, ViewMatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            super().__init__(matrix_hdf5, lsv_id)

        @property
        def dpsi(self):
            """
            Calculated the absolute difference in psi for heat map.
            :return: list
            """
            return [abs(reduce(operator.__sub__, (get_expected_psi(b) for b in bs))) for bs in self.mean_psi]

        @property
        def dpsi_signed(self):
            """
            Calculated the difference in psi for heat map. (with negative values possible)
            :return: list
            """
            return [reduce(operator.__sub__, (get_expected_psi(b) for b in bs)) for bs in self.mean_psi]

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
