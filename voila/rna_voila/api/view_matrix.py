import operator
from abc import ABC
from functools import reduce
from itertools import chain

import numpy as np
import scipy.stats

from rna_voila.api.matrix_hdf5 import DeltaPsi, Psi, Heterogen, MatrixType
from rna_voila.api.matrix_utils import unpack_means, unpack_bins, generate_excl_incl, generate_means, \
    generate_high_probability_non_changing, generate_variances, generate_standard_deviations
from rna_voila.config import ViewConfig
from rna_voila.exceptions import LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile, LsvIdNotFoundInAnyVoilaFile
from rna_voila.vlsv import is_lsv_changing, matrix_area, get_expected_psi
from multiprocessing import Pool
from itertools import combinations

class ViewMatrix(ABC):
    group_names = None
    experiment_names = None
    gene_ids = None

    def check_group_consistency(self):
        config = ViewConfig()
        groups_defined = {}
        for f in config.voila_files:
            with ViewPsi(f) as m:

                for i, group_name in enumerate(m.group_names):
                    experiment_names = [exp for exp in m.experiment_names[i] if exp != '']
                    experiment_names.sort()
                    if group_name not in groups_defined:
                        groups_defined[group_name] = [experiment_names]
                    else:
                        groups_defined[group_name].append(experiment_names)

        warnings = []
        for group_name, experiment_groups in groups_defined.items():
            if not all(e == experiment_groups[0] for e in experiment_groups):
                warnings.append((group_name, experiment_groups))

        return warnings


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
        Get bins data from rna_voila file.
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

    def __init__(self, voila_file=None):
        if voila_file != None:
            print("Warning, view multipsi calling with specific voila file not supported, using all voila inputs")
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
            Get means data from rna_voila file.
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
            Get means data from rna_voila file.
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

        @property
        def standard_deviations(self):
            """
            Create variance data of bins data.
            :return: list
            """
            return generate_standard_deviations(self.bins)

    def lsv(self, lsv_id):
        """
        Get lsv object by lsv id.
        :param lsv_id: lsv id
        :return: lsv object
        """
        return self._ViewPsi(self, lsv_id)


class ViewDeltaPsi(DeltaPsi, ViewMatrix):
    def __init__(self, voila_file=None):
        """
        View for delta psi matrix.  This is used in creation of tsv and html files.
        """
        self.config = ViewConfig()
        if voila_file is None:
            super().__init__(self.config.voila_file)
        else:
            super().__init__(voila_file)

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

        def changing(self, changing_threshold, probability_changing_threshold, junc_i=None):
            """
            Yet another changing function wrapper, because I can't figure out how the legacy ones were even being used
            This is used for creating the changing column in voila classifier. It functions similar in API to the
            non changing function below

            """
            for _junc_i in range(len(self.bins)) if junc_i is None else [junc_i]:
                junc_bins = self.bins[_junc_i]
                if matrix_area(junc_bins, threshold=changing_threshold) >= probability_changing_threshold:
                    return True

            return False

        def high_probability_non_changing(self, non_changing_threshold=None, junc_i=None):
            """
            Get probability that junctions in an lsv aren't changing.
            :return: list
            """
            if non_changing_threshold is None:
                non_changing_threshold = self.config.non_changing_between_group_dpsi

            bins = self.bins if junc_i is None else [self.bins[junc_i]]

            result = generate_high_probability_non_changing(self.intron_retention, self.matrix_hdf5.prior,
                                                          non_changing_threshold, bins)

            return result if junc_i is None else result[0]

    def lsv(self, lsv_id):
        """
        Get delta psi object by lsv id.
        :param lsv_id: lsv id
        :return: delta psi object
        """
        return self._ViewDeltaPsi(self, lsv_id)






class ViewHeterogens(ViewMulti):

    def __init__(self, voila_file=None, group_order_override=None):
        if voila_file != None:
            print("Warning, view heterogen calling with specific voila file not supported, using all voila inputs")
        self.group_order_override = group_order_override
        super().__init__(ViewHeterogen)

    class _ViewHeterogens(_ViewMulti):
        def __init__(self, matrix_hdf5, lsv_id):
            super().__init__(matrix_hdf5, lsv_id, ViewHeterogen)
            self.matrix_hdf5 = matrix_hdf5
            self.lsv_id = lsv_id

        @property
        def groups_quantified(self):

            # mu_psi has shape (njunc, ngroups, max(nexperiments))
            mu_psi = np.array(self.mu_psi)[0]
            num_quant = (mu_psi >= 0).sum(axis=-1)
            return num_quant

        def get_attr(self, attr):
            """
            For attributes that exist is each het file, this will get all values and confirm they're all equal.
            :param attr: attribute found in het voila file.
            :return: attribute value
            """
            config = ViewConfig()
            if config.strict_indexing:
                s = set()
                for f in config.voila_files:
                    with ViewHeterogen(f) as m:
                        try:
                            het = m.lsv(self.lsv_id)
                            s.add(getattr(het, attr))
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                            pass
                assert len(s) == 1, s
                return s.pop()
            else:
                for f in config.voila_files:
                    with ViewHeterogen(f) as m:
                        try:
                            het = m.lsv(self.lsv_id)
                            return getattr(het, attr)
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                            pass

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

                        dpsi_value = het.dpsi_median_signed
                        dpsi_value = dpsi_value[junc_idx]

                        grp_names = m.group_names
                        grp_names.sort()

                        x = hets_grps.index(grp_names[0])
                        y = hets_grps.index(grp_names[1])

                        if x - y > 0:
                            s[x][y] = stat_value
                            s[y][x] = dpsi_value
                        else:
                            s[x][y] = dpsi_value
                            s[y][x] = stat_value

                    except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                        pass

            return s.tolist()

        @property
        def median_psi(self):
            return self._median_psi()

        def quantile_psi(self, quantile):
            return self._median_psi(quantile=quantile)

        def _median_psi(self, quantile=None):
            """ Find median_psi per group in all het voila files

            Returns
            -------
            np.array(shape=(num_groups, num_junctions))
                The median PSI per junction for each group. -1 if not quantified
            """
            voila_files = ViewConfig().voila_files
            group_names = self.matrix_hdf5.group_names
            juncs_len = len(self.junctions)
            grps_len = len(group_names)
            # fill result
            median_psi = np.empty((grps_len, juncs_len))
            median_psi.fill(-1)

            for f in voila_files:
                with ViewHeterogen(f) as m:
                    try:
                        het = m.lsv(self.lsv_id)

                        if quantile is None:
                            medians = het.median_psi()
                        else:
                            medians = het.quantile_psi(quantile)
                        # get index into medians (second axis) per group
                        for ndx, grp in enumerate(m.group_names):
                            idx = group_names.index(grp)  # index to median_psi
                            median_psi[idx, :] = medians[:, ndx]  # fill median_psi

                    except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                        pass
            return median_psi

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
        def junction_psisamples_stats(self):
            """
            Get stats from psisamples quantiles with values
            :return: generator key/value
            """
            config = ViewConfig()
            voila_files = config.voila_files
            for f in voila_files:
                with ViewHeterogen(f) as m:
                    het = m.lsv(self.lsv_id)
                    groups = '-'.join(m.group_names)
                    stat_names = m.stat_names
                    try:
                        for name, stat in zip(stat_names, het.junction_psisamples_stats.T):
                            if len(voila_files) == 1:
                                yield f"{name}_quantile", stat
                            else:
                                yield f"{groups} {name}_quantile", stat
                    except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                        pass

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
                    groups = '-'.join(m.group_names)
                    stat_names = m.stat_names
                    try:
                        for name, stat in zip(stat_names, het.junction_stats.T):
                            if len(voila_files) == 1:
                                yield name, stat
                            else:
                                yield groups + ' ' + name, stat
                    except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                        pass

        def changing(
            self,
            pvalue_threshold: float = 0.05,
            between_group_dpsi: float = 0.05
        ):
            """ Boolean of heuristic for changing heterogeneous events

            Parameters
            ----------
            pvalue_threshold: float
                Minimum p-value for which an LSV/junction can return true. Uses
                minimum p-value from all tests provided
            between_group_dpsi: float
                Maximum absolute difference in median values of PSI for which
                an LSV/junction can return true

            Returns
            -------
            Generator yielding group names and boolean array per junction
            """
            config = ViewConfig()
            voila_files = config.voila_files
            for f in voila_files:
                with ViewHeterogen(f) as m:
                    het = m.lsv(self.lsv_id)
                    groups = '-'.join(m.group_names)
                    try:
                        changing = het.changing(
                            pvalue_threshold=pvalue_threshold,
                            between_group_dpsi=between_group_dpsi
                        )
                        if len(voila_files) == 1:
                            yield "changing", changing
                        else:
                            yield f"{groups} changing", changing
                    except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                        pass

        def nonchanging(
            self,
            pvalue_threshold: float = 0.05,
            within_group_iqr: float = 0.10,
            between_group_dpsi: float = 0.05
        ):
            """ Boolean of heuristic for nonchanging heterogeneous events

            We define highly-confident non-changing events from MAJIQ Heterogen
            as being (1) above a nominal p-value threshold, (2) within-group
            variance is sufficiently low as measured by IQR, (3) between-group
            dPSI is sufficiently low as measured by difference in medians. We
            accept that between-group dPSI threshold may be redundant in
            combination with the other two thresholds.

            Parameters
            ----------
            pvalue_threshold: float
                Minimum p-value for which an LSV/junction can return true. Uses
                minimum p-value from all tests provided
            within_group_iqr: float
                Maximum IQR within a group for which an LSV/junction can return
                true
            between_group_dpsi: float
                Maximum absolute difference in median values of PSI for which
                an LSV/junction can return true

            Returns
            -------
            Generator yielding group names and boolean array per junction
            """
            config = ViewConfig()
            voila_files = config.voila_files
            for f in voila_files:
                with ViewHeterogen(f) as m:
                    het = m.lsv(self.lsv_id)
                    groups = '-'.join(m.group_names)
                    try:
                        nonchanging = het.nonchanging(
                            pvalue_threshold=pvalue_threshold,
                            within_group_iqr=within_group_iqr,
                            between_group_dpsi=between_group_dpsi
                        )
                        if len(voila_files) == 1:
                            yield "nonchanging", nonchanging
                        else:
                            yield f"{groups} nonchanging", nonchanging
                    except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                        pass

        @property
        def junction_scores(self):
            """
            This gets associates stat test score names with their values. (if any exist)
            :return: generator key/value
            """
            config = ViewConfig()
            voila_files = config.voila_files
            for f in voila_files:
                with ViewHeterogen(f) as m:
                    groups = '-'.join(m.group_names)
                    score_names = ('tnom_score',) if 'TNOM' in m.stat_names else ()
                    try:
                        try:
                            score_vals = m.get(self.lsv_id, 'tnom_score').T
                        except KeyError:
                            score_names = ()

                        for score_name in score_names:
                            if len(voila_files) == 1:
                                yield score_name, score_vals
                            else:
                                yield groups + ' ' + score_name, score_vals

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
    def psi_samples_summary(self):
        """
        Summarize number of psi_samples used for posterior samples statistics

        :return: str indicating:
            + unknown: when all samples are missing this value in their metadata
            + inconsistent: different values found (or some missing/not)
            + value shared/consistent for all comparisons
        """
        values = set()  # type: Set[str]
        voila_files = ViewConfig().voila_files
        for f in voila_files:
            with ViewHeterogen(f) as m:
                try:
                    values.add(str(m.psi_samples))
                except KeyError:
                    values.add(
                        "missing metadata (older version of MAJIQ)"
                    )
        if len(values) > 1:
            return f"inconsistent ({sorted(values)})"
        elif values:  # only 1
            return values.pop()
        else:
            return f"no metadata (input files: {voila_files})"

    @property
    def test_percentile_summary(self):
        """
        Summarize percentile taken over psi-sample statistics

        :return: str indicating:
            + unknown: when all samples are missing this value in their metadata
            + inconsistent: different values found (or some missing/not)
            + value shared/consistent for all comparisons
        """
        values = set()  # type: Set[str]
        voila_files = ViewConfig().voila_files
        for f in voila_files:
            with ViewHeterogen(f) as m:
                try:
                    values.add(str(m.test_percentile))
                except KeyError:
                    values.add(
                        "missing metadata (older version of MAJIQ)"
                    )
        if len(values) > 1:
            return f"inconsistent ({sorted(values)})"
        elif values:  # only 1
            return values.pop()
        else:
            return f"no metadata (input files: {voila_files})"

    @property
    def junction_psisamples_stats_column_names(self):
        return (f"{x}_quantile" for x in self.junction_stats_column_names)

    @property
    def junction_stats_column_names(self):
        """
        Stat column names for tsv output.
        :return: generator
        """
        voila_files = ViewConfig().voila_files

        for f in voila_files:
            with ViewHeterogen(f) as m:
                groups = '-'.join(m.group_names)
                for name in m.stat_names:
                    if len(voila_files) == 1:
                        yield name
                    else:
                        yield groups + ' ' + name

    def _per_group_column_names(self, prefix):
        voila_files = ViewConfig().voila_files

        if len(voila_files) == 1:
            yield prefix
        else:
            for f in voila_files:
                with ViewHeterogen(f) as m:
                    groups = '-'.join(m.group_names)
                    yield f"{groups} {prefix}"

    @property
    def changing_column_names(self):
        """ Column names associated with changing values
        """
        return self._per_group_column_names('changing')

    @property
    def nonchanging_column_names(self):
        """ Column names associated with nonchanging values
        """
        return self._per_group_column_names('nonchanging')

    @property
    def junction_scores_column_names(self):
        """
        Stat column names for tsv output.
        :return: generator
        """
        voila_files = ViewConfig().voila_files

        for f in voila_files:
            with ViewHeterogen(f) as m:
                groups = '-'.join(m.group_names)
                score_names = ('tnom_score',) if 'TNOM' in m.stat_names else ()
                for name in score_names:
                    if len(voila_files) == 1:
                        yield name
                    else:
                        yield groups + ' ' + name

    @property
    def experiment_names(self):
        """
        Experiment names for this set of het voila files.
        :return: List
        """
        config = ViewConfig()
        exp_names = {}
        for f in config.voila_files:
            with ViewHeterogen(f) as m:
                for exp, grp in zip(m.experiment_names, m.group_names):
                    exp_names[grp] = exp

        return [exp_names[grp] for grp in self.group_names]

    @property
    def group_names(self):
        """
        Group names for this set of het voila files.
        :return: list
        """
        if self.group_order_override:
            return self.group_order_override
        config = ViewConfig()
        grp_names = []
        for f in config.voila_files:
            with ViewHeterogen(f) as m:
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
            with ViewHeterogen(f) as m:
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
        vhs = [ViewHeterogen(f) for f in voila_files]
        yield from set(chain(*(v.gene_ids for v in vhs)))
        for v in vhs:
            v.close()

    @staticmethod
    def pair_merge(pairs):
        if len(pairs) == 1:
            if type(pairs[0]) is set:
                return pairs[0]
            else:
                with ViewHeterogen(pairs[0]) as v:
                    return set(v.lsv_ids())
        else:
            if type(pairs[0]) is set:
                s1 = pairs[0]
            else:
                with ViewHeterogen(pairs[0]) as v:
                    s1 = set(v.lsv_ids())
            if type(pairs[1]) is set:
                s2 = pairs[1]
            else:
                with ViewHeterogen(pairs[1]) as v:
                    s2 = set(v.lsv_ids())
            return s1 | s2


    def lsv_ids(self, gene_ids=None):
        """
        Get a set of lsv ids from all voila files for specified gene ids. If gene ids is None, then get all lsv ids.
        :param gene_ids: list of gene ids
        :return:
        """
        if gene_ids or len(ViewConfig().voila_files) == 1:
            voila_files = ViewConfig().voila_files
            vhs = [ViewHeterogen(f) for f in voila_files]
            yield from set(chain(*(v.lsv_ids(gene_ids) for v in vhs)))
            for v in vhs:
                v.close()
        else:
            vhs = ViewConfig().voila_files
            p = Pool(ViewConfig().nproc)
            while len(vhs) > 1:
                vhs = [vhs[i:i + 2] for i in range(0, len(vhs), 2)]
                vhs = p.map(self.pair_merge, vhs)

            yield from vhs[0]

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
        def dpsi_median_signed(self):
            """
            Calculated the difference in median psi for heat map. (with negative values possible)
            :return: list
            """
            return [reduce(operator.__sub__, bs) for bs in self.median_psi()]

        @property
        def mean_psi(self):
            return self.get('mean_psi')

        @property
        def mu_psi(self):
            return self.get('mu_psi')

        @property
        def mu_psi_nanmasked(self):
            """ Mask missing/unquantified values in mu_psi with nan
            """
            mu_psi = self.mu_psi  # load mu psi ndarray into memory
            return np.where(mu_psi >= 0, mu_psi, np.nan)  # mask unquantified

        @property
        def junction_stats(self):
            """ Junction statistics computed on posterior means
            """
            return self.get('junction_stats')

        @property
        def junction_psisamples_stats(self):
            """ Quantile from junction statistics on posterior samples
            """
            return self.get('junction_psisamples_stats')

        def median_psi(self, mu_psi=None):
            """ Get group medians of psi_mean (per junction if present)

            Parameters
            ----------
            mu_psi: Optional[np.array(shape=(..., 2, num_experiments))]
                Values to compute median per group, with unquantified values
                masked as nan. If not specified, use self.mu_psi_nanmasked.
                Preceding axes typically correspond to junctions.

            Returns
            -------
            np.array(shape=(..., 2))
                Median value of quantified values of PSI per group/junction
            """
            if mu_psi is None:
                mu_psi = self.mu_psi_nanmasked
            return np.nanmedian(mu_psi, axis=-1)

        def quantile_psi(self, quantile):
            return np.nanquantile(self.mu_psi_nanmasked, quantile, axis=-1)

        def iqr_psi(self, mu_psi=None):
            """ Get group IQRs of psi_mean (per junction if present)

            Parameters
            ----------
            mu_psi: Optional[np.array(shape=(..., 2, num_experiments))]
                Values to compute median per group, with unquantified values
                masked as nan. If not specified, use self.mu_psi_nanmasked.
                Preceding axes typically correspond to junctions.

            Returns
            -------
            np.array(shape=(..., 2))
                IQR of PSI per group/junction
            """
            if mu_psi is None:
                mu_psi = self.mu_psi_nanmasked
            return scipy.stats.iqr(mu_psi, axis=-1, nan_policy="omit")


        def dpsi_median(self,
                         junc_i: int = None):

            if junc_i is None:
                junc_i = slice(None)  # vectorize over all junctions

            # get masked mean values of psi per group/experiment
            mu_psi = self.mu_psi_nanmasked[junc_i]
            # use to compute median psi per group
            median_psi = self.median_psi(mu_psi)
            # did difference in medians pass threshold?
            dpsi_median = np.abs(median_psi[..., 1] - median_psi[..., 0])
            return dpsi_median

        def changing(
            self,
            pvalue_threshold: float = 0.05,
            between_group_dpsi: float = 0.2,
            junc_i: int = None
        ):


            if junc_i is None:
                junc_i = slice(None)  # vectorize over all junctions

            # all statistics must be less than p-value threshold
            pvalue_passed = np.nanmax(self.junction_stats[junc_i], axis=-1) <= pvalue_threshold

            dpsi_passed = self.dpsi_median(junc_i) >= between_group_dpsi

            # pvalue and dpsi thresholds must all pass

            return pvalue_passed & dpsi_passed

        def nonchanging(
            self,
            pvalue_threshold: float = 0.05,
            within_group_iqr: float = 0.10,
            between_group_dpsi: float = 0.05,
            junc_i: int = None
        ):
            """ Boolean of heuristic for nonchanging heterogeneous events

            We define highly-confident non-changing events from MAJIQ Heterogen
            as being (1) above a nominal p-value threshold, (2) within-group
            variance is sufficiently low as measured by IQR, (3) between-group
            dPSI is sufficiently low as measured by difference in medians. We
            accept that between-group dPSI threshold may be redundant in
            combination with the other two thresholds.

            Parameters
            ----------
            pvalue_threshold: float
                Minimum p-value for which an LSV/junction can return true. Uses
                minimum p-value from all tests provided
            within_group_iqr: float
                Maximum IQR within a group for which an LSV/junction can return
                true
            between_group_dpsi: float
                Maximum absolute difference in median values of PSI for which
                an LSV/junction can return true
            junc_i: int
                Only calculate result for one junction from the LSV at the specified index
                None (default) returns a list of
            """
            if junc_i is None:
                junc_i = slice(None)  # vectorize over all junctions

            # all statistics must be greater than p-value threshold
            pvalue_passed = np.nanmin(self.junction_stats[junc_i], axis=-1) >= pvalue_threshold
            # mean values of PSI per group/experiment for requested junctions
            # missing quantifications are masked as nan
            mu_psi = self.mu_psi_nanmasked[junc_i]
            # use to compute IQR and median
            iqr_psi = self.iqr_psi(mu_psi)
            median_psi = self.median_psi(mu_psi)
            # did IQR pass threshold?
            iqr_passed = (iqr_psi <= within_group_iqr).all(axis=-1)  # both groups
            # did difference in medians pass threshold?
            dpsi_passed = np.abs(median_psi[..., 1] - median_psi[..., 0]) <= between_group_dpsi

            # pvalue, iqr, and dpsi thresholds must all pass
            return pvalue_passed & iqr_passed & dpsi_passed

    def lsv(self, lsv_id):
        return self._ViewHeterogen(self, lsv_id)
