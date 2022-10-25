import statistics
from rna_voila.config import ClassifyConfig
import numpy as np
from rna_voila.vlsv import get_expected_psi, matrix_area
from itertools import combinations
from operator import itemgetter
from rna_voila.api import Matrix
from rna_voila.api.view_matrix import ViewHeterogen
from rna_voila import constants
from rna_voila.exceptions import GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile
from rna_voila.api.matrix_utils import generate_variances
from rna_voila.api import view_matrix
from collections import OrderedDict
from rna_voila.api import SpliceGraph
from statistics import median, StatisticsError
from math import ceil

def fRound(x):
    try:
        x = iter(x)
        ret = ';'.join(f'{_x:0.3e}' for _x in x)
    except TypeError as te:
        ret = f'{x:0.3e}'
    return ret

class QuantificationWriter:

    def __init__(self):
        """

        :param avg_multival: if true, when finding quantifications with multiple matches, instead of semicolon, avg them
        """

        self.config = ClassifyConfig()
        self.dpsi_quant_idxs = []
        self.quantifications_int = self.quantification_intersection()
        self.avg_multival = False


    @staticmethod
    def semicolon(value_list):
        return ';'.join(str(x) for x in value_list)

    def quantifications(self, module, parity=None, edge=None, node=None):
        """
        Edge / Parity is used to find LSVs
        Node is used to filter lsvs to specific node (the node that has THAT lsv)
        :return:
        """


        lsvs = self.parity2lsv(module, parity, edge=edge)


        out = []
        for field in self.quantifications_int:

            quantification_vals = []
            for lsv_id in lsvs:
                 #print(self.quantifications_int[field](lsv_id, edge))

                quants = self.quantifications_int[field][0](*self.quantifications_int[field][1:])(lsv_id, edge)
                if quants is None:
                    quantification_vals.append('')
                else:
                    for val in quants:
                        quantification_vals.append(val)

            out.append(self.semicolon(quantification_vals))


        return out

    def parity2lsv_node(self, module, parity, node=None):

        if parity == 's':
            lsvs = module.source_lsv_ids
            if node:
                lsvs = set(filter(lambda lsv: lsv.endswith(node.untrimmed_range_str()), lsvs))
        elif parity == 't':
            lsvs = module.target_lsv_ids
            if node:
                lsvs = set(filter(lambda lsv: lsv.endswith(node.untrimmed_range_str()), lsvs))
        else:
            lsvs = module.target_lsv_ids.union(module.source_lsv_ids)
        return lsvs

    def parity2lsv(self, module, parity, edge=None):

        if parity == 's':
            lsv_ids_mod = module.source_lsv_ids
        elif parity == 't':
            lsv_ids_mod = module.target_lsv_ids
        else:
            lsv_ids_mod = module.target_lsv_ids.union(module.source_lsv_ids)

        if edge is not None:
            lsvs = set()
            edges = [edge] if type(edge) is not list else edge
            for _edge in edges:
                lsvs = lsvs.union(set([lsv_id for lsv_id in _edge.lsvs if lsv_id in lsv_ids_mod]))
            return lsvs

        return set(lsv_ids_mod)

    def edge_quant(self, module, edge, field):
        """
        Get one quantification number for a specific edge
        """
        lsvs = self.parity2lsv(module, None, edge=edge)
        to_avg = []
        for lsv_id in lsvs:

            vals = self.quantifications_int[field][0](*self.quantifications_int[field][1:])(lsv_id, edge)

            for val in vals:
                to_avg.append(val)
        # print(to_avg)
        if to_avg:
            return np.mean(to_avg)
        else:
            return -1

    def _filter_edges(self, edge, lsv):
        if type(edge) != list:
            edge = [edge]
        for _edge in edge:
            # loop through junctions to find one matching range of edge
            try:
                for j, junc in enumerate(lsv.get('junctions')):
                    if junc[0] == _edge.absolute_start and junc[1] == _edge.absolute_end:
                        return j
                else:
                    # junction not quantified by majiq
                    pass
            except:
                pass
        return slice(None)

    def quantification_intersection(self):
        """
        Look at all psi and dpsi quant headers and find the appropriate intersection
        we need to then define which function should be called for each resulting column

        This is likely a very confusing section so what is going on requires some context.

        Basically, for each combination group name + stat, we only want to show one column for it in the TSV
        However, when looking at all files, we may come across it twice. In order to determine if we have
        the information for it in a specific voila file, we need to follow a specific algorithm to get the data from
        the voila file in the first place.

        So, we loop over all possible voila files associated with that stat, and once we find a valid value, we
        return it.

        The bottom half of this function loops over all voila files to build up a list of stats keys to the functions
        that will be run for each event to get the required data for that event. (all of the "_" functions inside
        this function, return functions)

        This is an efficiency compromise, because we can build the list of functions once for a gene, and only need
        to open and read the voila files again when the quantification function is called.

        :return:
        """

        def _inner_edge_aggregate(lsv, all_quants, edge, _round=True):
            if edge:
                edges = [edge] if not type(edge) is list else edge
                vals = []
                for _edge in edges:
                    edge_idx = self._filter_edges(_edge, lsv)
                    if edge_idx is None:
                        continue
                    else:
                        vals.append(all_quants[edge_idx])

                return (fRound(x) for x in vals) if _round else vals
            else:
                if self.avg_multival and all_quants:
                    return np.mean(all_quants)
                return (fRound(x) for x in all_quants) if _round else all_quants

        def _psi_psi(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with Matrix(voila_file) as m:
                        try:
                            lsv = m.psi(lsv_id)
                            return _inner_edge_aggregate(lsv, lsv.get('means'), edge)
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                return None
            return f

        def _psi_var(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with Matrix(voila_file) as m:
                        try:
                            lsv = m.psi(lsv_id)
                            return _inner_edge_aggregate(lsv, generate_variances([lsv.get('bins')][0]), edge)
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                return None
            return f

        def _het_psi(voila_files, group_idxs):
            def f(lsv_id, edge=None):
                found_psis = []
                for voila_file, group_idx in zip(voila_files, group_idxs):
                    with ViewHeterogen(voila_file) as m:
                        try:
                            lsv = m.lsv(lsv_id)
                            edge_idx = self._filter_edges(edge, lsv)
                            if edge_idx is None:
                                continue
                            medians = lsv.median_psi()
                            psi = medians[edge_idx][group_idx]
                            found_psis.append(psi)
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                if not found_psis:
                    return None
                else:
                    return [fRound(sum(found_psis) / len(found_psis))]
            return f

        def _het_stats(voila_files, stat_idx):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with view_matrix.ViewHeterogen(voila_file) as m:
                        try:
                            lsv = m.heterogen(lsv_id)
                            return _inner_edge_aggregate(lsv, [x[stat_idx] for x in m.lsv(lsv_id).junction_stats], edge)
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                return None
            return f

        def _het_dpsi(voila_files, group_idxs1, group_idxs2):
            def f(lsv_id, edge):
                for voila_file, group_idx1, group_idx2 in zip(voila_files, group_idxs1, group_idxs2):
                    with ViewHeterogen(voila_file) as m:
                        try:
                            # for this one the _inner_edge_aggregate is not general enough - I had to do it manually
                            lsv = m.lsv(lsv_id)

                            edges = [edge] if not type(edge) is list else edge
                            vals = []

                            for _edge in edges:
                                edge_idx = self._filter_edges(_edge, lsv)
                                if edge_idx is None:
                                    continue
                                else:
                                    medians = lsv.median_psi()

                                    psi_g1 = medians[edge_idx][group_idx1]
                                    psi_g2 = medians[edge_idx][group_idx2]

                                    vals.append(psi_g1-psi_g2)

                            return (fRound(x) for x in vals)

                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                return None
            return f


        def _dpsi_psi(voila_files, group_idx):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with Matrix(voila_file) as m:
                        try:
                            lsv = m.delta_psi(lsv_id)
                            return _inner_edge_aggregate(lsv, lsv.get('group_means')[group_idx], edge)
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                return None
            return f

        def _dpsi_dpsi(voila_files):
            def f(lsv_id, edge):
                for voila_file in voila_files:
                    with view_matrix.ViewDeltaPsi(voila_file) as m:
                        try:
                            # for this one the _inner_edge_aggregate is not general enough - I had to do it manually
                            lsv = m.lsv(lsv_id)
                            bins = lsv.get('group_bins')

                            edges = [edge] if not type(edge) is list else edge
                            vals = []
                            for _edge in edges:
                                edge_idx = self._filter_edges(_edge, lsv)
                                if edge_idx is None:
                                    continue
                                else:
                                    vals.append(lsv.excl_incl[edge_idx][1] - lsv.excl_incl[edge_idx][0])

                            return (fRound(x) for x in vals)

                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                return None
            return f

        def _dpsi_p_change(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with view_matrix.ViewDeltaPsi(voila_file) as m:
                        try:
                            # for this one the _inner_edge_aggregate is not general enough - I had to do it manually
                            lsv = m.lsv(lsv_id)
                            bins = lsv.bins
                            if edge:
                                edges = [edge] if not type(edge) is list else edge
                                vals = []
                                for _edge in edges:
                                    edge_idx = self._filter_edges(_edge, lsv)
                                    if edge_idx is None:
                                        continue
                                    else:
                                        vals.append(matrix_area(bins[edge_idx], self.config.changing_between_group_dpsi))
                                return (fRound(x) for x in vals)
                            else:
                                if self.avg_multival:
                                    return np.mean((matrix_area(b, self.config.changing_between_group_dpsi) for b in bins))
                                return (
                                            fRound(matrix_area(b, self.config.changing_between_group_dpsi)) for b in bins
                                        )
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                return None
            return f

        def _dpsi_p_nonchange(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with view_matrix.ViewDeltaPsi(voila_file) as m:
                        try:
                            lsv = m.lsv(lsv_id)
                            return _inner_edge_aggregate(lsv, lsv.high_probability_non_changing(), edge)
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                            continue
                return None
            return f

        def _junction_changing(voila_files):

            # should only have one edge specified --
            # iterate through voila files, if we fine any case where junction is changing,
            # return true

            def f(lsv_id, edge=None):

                found_quant = False
                edges = [edge] if not type(edge) is list else edge

                for voila_file in voila_files:
                    try:
                        with Matrix(voila_file) as m1:
                            analysis_type = m1.analysis_type

                        if edge:

                            for _edge in edges:

                                if analysis_type == constants.ANALYSIS_HETEROGEN:
                                    with view_matrix.ViewHeterogen(voila_file) as m:
                                        lsv = m.lsv(lsv_id)

                                        edge_idx = self._filter_edges(_edge, lsv)
                                        if edge_idx is None:
                                            continue
                                        else:

                                            is_changing = lsv.changing(
                                                         self.config.changing_pvalue_threshold,
                                                         self.config.changing_between_group_dpsi,
                                                         edge_idx)
                                            found_quant = True

                                            if bool(is_changing):
                                                return [True]

                                elif analysis_type == constants.ANALYSIS_DELTAPSI:
                                    with view_matrix.ViewDeltaPsi(voila_file) as m:
                                        lsv = m.lsv(lsv_id)

                                        edge_idx = self._filter_edges(_edge, lsv)
                                        if edge_idx is None:
                                            continue
                                        else:
                                            is_changing = lsv.changing(
                                                self.config.changing_between_group_dpsi,
                                                self.config.probability_changing_threshold,
                                                edge_idx)
                                            found_quant = True

                                            if bool(is_changing):
                                                return [True]


                    except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                        continue

                return [False] if found_quant else ""

            return f

        def _reads(splice_graph_file, gene_id, _experiment_names):
            def f(lsv_id, edge=None):
                with SpliceGraph(splice_graph_file) as sg:

                    try:
                        junc = {'start': edge.start, 'end': edge.end, 'gene_id': gene_id}
                        reads = ceil(median((x['reads'] for x in sg.junction_reads_exp(junc, _experiment_names))))
                    except:
                        reads = ''

                return [reads]

            return f

        tmp = OrderedDict()
        self.types2headers = {'psi':[], 'dpsi':[]}

        tmp['junction_changing'] = (_junction_changing, self.config.voila_files)

        # junc {'gene_id': 'ENSMUSG00000001419', 'start': 88168458, 'end': 88168632, 'has_reads': 1, 'annotated': 1, 'is_simplified': 0, 'is_constitutive': 0}



        for voila_file in self.config.voila_files:

            with Matrix(voila_file) as m:
                analysis_type = m.analysis_type
                group_names = m.group_names
                experiment_names = m.experiment_names
                if analysis_type == constants.ANALYSIS_HETEROGEN:
                    stat_names = m.stat_names
                else:
                    stat_names = None



            if analysis_type == constants.ANALYSIS_PSI:
                for group, experiments in zip(group_names, experiment_names[:-1]):

                    if self.config.show_read_counts:

                        header = f'{group}_median_reads'
                        tmp[header] = (_reads, self.config.splice_graph_file, self.graph.gene_id if self.graph else None, experiments)


                for key in ("median_psi", "var_psi",):
                        header = "%s_%s" % (group, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                        else:
                            if key == "median_psi":
                                tmp[header] = (_psi_psi, [voila_file])
                            elif key == "var_psi":
                                tmp[header] = (_psi_var, [voila_file])


            elif analysis_type == constants.ANALYSIS_HETEROGEN:

                group_idxs = {}
                for i, group in enumerate(group_names):
                    group_idxs[group] = i
                    for key in ("median_psi",):
                        header = "%s_het_%s" % (group, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                            tmp[header][2].append(i)
                        else:
                            if key == "median_psi":
                                tmp[header] = (_het_psi, [voila_file], [i])



                for group1, group2 in combinations(group_names, 2):
                    for key in ("median_dpsi",):
                        header = "%s-%s_het_%s" % (group1, group2, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                            tmp[header][2].append(group_idxs[group1])
                            tmp[header][3].append(group_idxs[group2])
                        else:
                            if key == "median_dpsi":
                                self.dpsi_quant_idxs.append(len(tmp))
                                tmp[header] = (_het_dpsi, [voila_file], [group_idxs[group1]], [group_idxs[group2]])

                    for j, key in enumerate(stat_names):
                        header = "%s-%s_het_%s" % (group1, group2, key.lower())
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                        else:
                            tmp[header] = (_het_stats, [voila_file], j)


            else:
                for i, group in enumerate(group_names):
                    for key in ("median_psi",):
                        header = "%s_%s" % (group, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                        else:
                            if key == "median_psi":
                                tmp[header] = (_dpsi_psi, [voila_file], i)
                        self.types2headers['psi'].append(header)

                changing_thresh_key = "probability_changing"
                non_changing_thresh_key = "probability_non_changing"
                for key in ("median_dpsi", changing_thresh_key, non_changing_thresh_key):
                    header = "%s_%s" % ('-'.join(reversed(group_names)), key)
                    if header in tmp:
                        tmp[header][1].append(voila_file)
                    else:
                        if key == "median_dpsi":
                            self.dpsi_quant_idxs.append(len(tmp))
                            tmp[header] = (_dpsi_dpsi, [voila_file])
                            self.types2headers['dpsi'].append(header)

                        elif key == changing_thresh_key:
                            tmp[header] = (_dpsi_p_change, [voila_file])
                        elif key == non_changing_thresh_key:
                            tmp[header] = (_dpsi_p_nonchange, [voila_file])

        return tmp


class MultiQuantWriter(QuantificationWriter):

    def __init__(self):

        super().__init__()
        self.config = ClassifyConfig()

    def gen_lsvs_list(self, module, quant_identifiers):
        lsvs = []
        missing_any = False
        for parity, edge in quant_identifiers:
            lsv_ids = self.parity2lsv(module, parity, edge=edge)
            if lsv_ids:
                lsvs.append((lsv_ids.pop(), edge,))
            else:
                missing_any = True
        return lsvs, missing_any

    def _reads(self, splice_graph_file, gene_id, _experiment_names, edge):

        with SpliceGraph(splice_graph_file) as sg:

            junc = {'start': edge.start, 'end': edge.end, 'gene_id': gene_id}
            try:
                reads = ceil(median((x['reads'] for x in sg.junction_reads_exp(junc, _experiment_names))))
            except statistics.StatisticsError:
                return None
        return reads

    def event_changing(self, module, quant_identifiers):

        # should only have one edge specified --
        # iterate through voila files, if we fine any case where junction is changing,
        # return true
        lsvs, missing_any = self.gen_lsvs_list(module, quant_identifiers)
        if missing_any:
            return False
        found_quant = False

        for voila_file in self.config.voila_files:

            with Matrix(voila_file) as m1:
                analysis_type = m1.analysis_type

            junc_results = []

            for lsv_id, _edge in lsvs:

                try:

                    if analysis_type == constants.ANALYSIS_HETEROGEN:
                        with view_matrix.ViewHeterogen(voila_file) as m:
                            lsv = m.lsv(lsv_id)

                            edge_idx = self._filter_edges(_edge, lsv)

                            if edge_idx is None:
                                break
                            else:

                                is_changing = lsv.changing(
                                    self.config.changing_pvalue_threshold,
                                    self.config.changing_between_group_dpsi,
                                    edge_idx)

                                junc_results.append(is_changing)

                                is_changing_secondary = lsv.changing(
                                    1.0,
                                    self.config.changing_between_group_dpsi_secondary,
                                    edge_idx)

                                found_quant = True

                                if not is_changing_secondary:
                                    break

                    elif analysis_type == constants.ANALYSIS_DELTAPSI:
                        with view_matrix.ViewDeltaPsi(voila_file) as m:
                            lsv = m.lsv(lsv_id)

                            edge_idx = self._filter_edges(_edge, lsv)
                            if edge_idx is None:
                                continue
                            else:
                                is_changing = lsv.changing(
                                    self.config.changing_between_group_dpsi,
                                    self.config.probability_changing_threshold,
                                    edge_idx)

                                junc_results.append(is_changing)

                                is_changing_secondary = lsv.changing(
                                    self.config.changing_between_group_dpsi_secondary,
                                    self.config.probability_changing_threshold,
                                    edge_idx)

                                found_quant = True

                                if not is_changing_secondary:
                                    break

                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                    break

            else:
                # this executes is we never broke, in other words, all secondary filters passed
                # so we check if any primary ones passed, if this is not true, we move to the next
                # voila file...
                if any(bool(x) is True for x in junc_results):
                    return True

        return False if found_quant else ''


    def event_non_changing(self, module, quant_identifiers):

        lsvs, missing_any = self.gen_lsvs_list(module, quant_identifiers)
        if missing_any:
            return False

        junc_results = []

        for voila_file in self.config.voila_files:

            with Matrix(voila_file) as m1:
                analysis_type = m1.analysis_type
                experiments = m1.experiment_names[:-1][0]

            for lsv_id, _edge in lsvs:

                if self.config.non_changing_median_reads_threshold:
                    if not _edge:
                        continue
                    mean_reads = self._reads(self.config.splice_graph_file, self.graph.gene_id if self.graph else None, experiments, _edge)
                    if not mean_reads or mean_reads < self.config.non_changing_median_reads_threshold:
                        junc_results.append(False)
                        continue

                try:
                    if analysis_type == constants.ANALYSIS_HETEROGEN:
                        with view_matrix.ViewHeterogen(voila_file) as m:
                            lsv = m.lsv(lsv_id)

                            edge_idx = self._filter_edges(_edge, lsv)
                            if edge_idx is None:
                                continue
                            else:
                                is_non_changing = lsv.nonchanging(self.config.non_changing_pvalue_threshold,
                                                        self.config.non_changing_within_group_iqr,
                                                        self.config.non_changing_between_group_dpsi,
                                                        edge_idx)
                                junc_results.append(is_non_changing)

                    elif analysis_type == constants.ANALYSIS_DELTAPSI:
                        with view_matrix.ViewDeltaPsi(voila_file) as m:
                            lsv = m.lsv(lsv_id)

                            edge_idx = self._filter_edges(_edge, lsv)
                            if edge_idx is None:
                                continue
                            else:
                                non_changing_quant = lsv.high_probability_non_changing(
                                    self.config.non_changing_between_group_dpsi, edge_idx)

                                is_non_changing = non_changing_quant >= self.config.probability_non_changing_threshold

                                junc_results.append(is_non_changing)

                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                    continue

        if not junc_results:
            return ''

        # bool() needed here because they are of type "numpy._bool" by default
        return all(bool(x) is True for x in junc_results)

