from rna_voila.config import ClassifyConfig
import numpy as np
from rna_voila.vlsv import get_expected_psi, matrix_area
from itertools import combinations
from operator import itemgetter
from rna_voila.api import Matrix
from rna_voila import constants
from rna_voila.exceptions import GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile
from rna_voila.api.matrix_utils import generate_variances
from rna_voila.api import view_matrix
from collections import OrderedDict

SIG_FIGS = 3

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

    def quantifications(self, module, parity=None, edge=None, node=None, prev_quants=None):
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
                try:

                    #print(self.quantifications_int[field](lsv_id, edge))

                    quants = self.quantifications_int[field][0](*self.quantifications_int[field][1:])(lsv_id, edge)
                    if quants is None:
                        quantification_vals.append('')
                    else:
                        for val in quants:
                            quantification_vals.append(val)


                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                    quantification_vals.append('')


            out.append(self.semicolon(quantification_vals))

        if prev_quants is not None and len(prev_quants) > 1:

            if prev_quants[0] is True or out[0] is True:
                out[0] = True
            elif prev_quants[0] == '' and out[0] == '':
                out[0] = ''
            else:
                out[0] = False

            if prev_quants[1] is False or out[1] is False:
                out[1] = False
            elif prev_quants[1] == '' and out[1] == '':
                out[1] = ''
            else:
                out[1] = True

        return out

    def parity2lsv_node(self, module, parity, node=None):

        if parity == 's':
            lsvs = module.source_lsv_ids
            if node:
                lsvs = set(filter(lambda lsv: lsv.endswith(node.untrimmed_range_str(replace_w_na=True)), lsvs))
        elif parity == 't':
            lsvs = module.target_lsv_ids
            if node:
                lsvs = set(filter(lambda lsv: lsv.endswith(node.untrimmed_range_str(replace_w_na=True)), lsvs))
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


        def _filter_edges(edge, lsv):
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

        def _inner_edge_aggregate(lsv, all_quants, edge):
            if edge:
                edges = [edge] if not type(edge) is list else edge
                vals = []
                for _edge in edges:
                    edge_idx = _filter_edges(_edge, lsv)
                    if edge_idx is None:
                        continue
                    else:
                        vals.append(all_quants[edge_idx])

                return (round(x, SIG_FIGS) for x in vals)
            else:
                if self.avg_multival and all_quants:
                    return np.mean(all_quants)
                return (round(x, SIG_FIGS) for x in all_quants)

        def _psi_psi(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with Matrix(voila_file) as m:
                        lsv = m.psi(lsv_id)
                        return _inner_edge_aggregate(lsv, lsv.get('means'), edge)
                return None
            return f

        def _psi_var(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with Matrix(voila_file) as m:
                        lsv = m.psi(lsv_id)
                        return _inner_edge_aggregate(lsv, generate_variances([lsv.get('bins')][0]), edge)
                return None
            return f

        def _het_psi(voila_files, group_idx):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with Matrix(voila_file) as m:
                        lsv = m.heterogen(lsv_id)
                        return _inner_edge_aggregate(lsv, [get_expected_psi(x) for x in np.array(list(lsv.get('mean_psi'))).transpose((1, 0, 2))[group_idx]], edge)
                return None
            return f


        def _het_stats(voila_files, stat_idx):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with view_matrix.ViewHeterogen(voila_file) as m:
                        lsv = m.heterogen(lsv_id)
                        return _inner_edge_aggregate(lsv, [x[stat_idx] for x in m.lsv(lsv_id).junction_stats], edge)
                return None
            return f

        def _het_dpsi(voila_files, group_idx1, group_idx2):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with Matrix(voila_file) as m:
                        # for this one the _inner_edge_aggregate is not general enough - I had to do it manually
                        lsv = m.heterogen(lsv_id)
                        if edge:
                            edges = [edge] if not type(edge) is list else edge
                            vals = []

                            for _edge in edges:
                                edge_idx = _filter_edges(_edge, lsv)
                                if edge_idx is None:
                                    continue
                                else:
                                    arr = np.array(list(lsv.get('mean_psi'))).transpose((1, 0, 2))
                                    psi_g1 = get_expected_psi(arr[group_idx1][edge_idx])
                                    psi_g2 = get_expected_psi(arr[group_idx2][edge_idx])
                                    vals.append(psi_g1-psi_g2)

                            return (round(x, SIG_FIGS) for x in vals)
                        else:
                            group_means = []
                            arr = np.array(list(lsv.get('mean_psi'))).transpose((1, 0, 2))
                            psis_g1 = arr[group_idx1]
                            psis_g2 = arr[group_idx2]
                            for psi_g1, psi_g2 in zip(psis_g1, psis_g2):
                                group_means.append(get_expected_psi(psi_g1) - get_expected_psi(psi_g2))
                            if self.avg_multival:
                                return np.mean(group_means)
                            return (round(x, SIG_FIGS) for x in group_means)
                return None
            return f


        def _dpsi_psi(voila_files, group_idx):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with Matrix(voila_file) as m:
                        lsv = m.delta_psi(lsv_id)
                        return _inner_edge_aggregate(lsv, lsv.get('group_means')[group_idx], edge)
                return None
            return f

        def _dpsi_dpsi(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with view_matrix.ViewDeltaPsi(voila_file) as m:
                        # for this one the _inner_edge_aggregate is not general enough - I had to do it manually
                        lsv = m.lsv(lsv_id)
                        bins = lsv.get('group_bins')

                        if edge:
                            edges = [edge] if not type(edge) is list else edge
                            vals = []
                            for _edge in edges:
                                edge_idx = _filter_edges(_edge, lsv)
                                if edge_idx is None:
                                    continue
                                else:
                                    vals.append(lsv.excl_incl[edge_idx][1] - lsv.excl_incl[edge_idx][0])

                            return (round(x, SIG_FIGS) for x in vals)
                        else:
                            if self.avg_multival:
                                return np.mean((lsv.excl_incl[i][1] - lsv.excl_incl[i][0] for i in range(np.size(bins, 0))))
                            return (
                                        round(x, SIG_FIGS) for x in ((lsv.excl_incl[i][1] - lsv.excl_incl[i][0] for i in
                                        range(np.size(bins, 0))))
                                    )
                return None
            return f

        def _dpsi_p_change(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with view_matrix.ViewDeltaPsi(voila_file) as m:
                        # for this one the _inner_edge_aggregate is not general enough - I had to do it manually
                        lsv = m.lsv(lsv_id)
                        bins = lsv.bins
                        if edge:
                            edges = [edge] if not type(edge) is list else edge
                            vals = []
                            for _edge in edges:
                                edge_idx = _filter_edges(_edge, lsv)
                                if edge_idx is None:
                                    continue
                                else:
                                    vals.append(matrix_area(bins[edge_idx], self.config.changing_between_group_dpsi))
                            return (round(x, SIG_FIGS) for x in vals)
                        else:
                            if self.avg_multival:
                                return np.mean((matrix_area(b, self.config.changing_between_group_dpsi) for b in bins))
                            return (
                                        round(matrix_area(b, self.config.changing_between_group_dpsi), SIG_FIGS) for b in bins
                                    )
                return None
            return f

        def _dpsi_p_nonchange(voila_files):
            def f(lsv_id, edge=None):
                for voila_file in voila_files:
                    with view_matrix.ViewDeltaPsi(voila_file) as m:
                        lsv = m.lsv(lsv_id)
                        return _inner_edge_aggregate(lsv, lsv.high_probability_non_changing(), edge)
                return None
            return f

        def _global_changing(voila_files):
            def f(lsv_id, edge=None):

                junc_results = []
                edges = [edge] if not type(edge) is list else edge

                for voila_file in voila_files:
                    with Matrix(voila_file) as m1:
                        analysis_type = m1.analysis_type

                    if edge:

                        for _edge in edges:

                            if analysis_type == constants.ANALYSIS_HETEROGEN:
                                with view_matrix.ViewHeterogen(voila_file) as m:
                                    lsv = m.lsv(lsv_id)

                                    edge_idx = _filter_edges(_edge, lsv)
                                    if edge_idx is None:
                                        continue
                                    else:

                                        is_changing = lsv.changing(
                                                     self.config.changing_pvalue_threshold,
                                                     self.config.changing_between_group_dpsi,
                                                     edge_idx)

                                        junc_results.append(is_changing)

                            elif analysis_type == constants.ANALYSIS_DELTAPSI:
                                with view_matrix.ViewDeltaPsi(voila_file) as m:
                                    lsv = m.lsv(lsv_id)

                                    edge_idx = _filter_edges(_edge, lsv)
                                    if edge_idx is None:
                                        continue
                                    else:
                                        is_changing = lsv.changing(
                                            self.config.changing_between_group_dpsi,
                                            self.config.probability_changing_threshold,
                                            edge_idx)

                                        junc_results.append(is_changing)

                if not junc_results:
                    return ''

                # bool() needed here because they are of type "numpy._bool" by default
                return [any(bool(x) is True for x in junc_results)]

            return f

        def _global_non_changing(voila_files):
            def f(lsv_id, edge=None):

                junc_results = []
                edges = [edge] if not type(edge) is list else edge

                for voila_file in voila_files:
                    with Matrix(voila_file) as m1:
                        analysis_type = m1.analysis_type

                    if edge:

                        for _edge in edges:

                            if analysis_type == constants.ANALYSIS_HETEROGEN:
                                with view_matrix.ViewHeterogen(voila_file) as m:
                                    lsv = m.lsv(lsv_id)

                                    edge_idx = _filter_edges(_edge, lsv)
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

                                    edge_idx = _filter_edges(_edge, lsv)
                                    if edge_idx is None:
                                        continue
                                    else:
                                        non_changing_quant = lsv.high_probability_non_changing(
                                            self.config.non_changing_threshold, edge_idx)

                                        is_non_changing = non_changing_quant >= self.config.probability_non_changing_threshold

                                        junc_results.append(is_non_changing)

                if not junc_results:
                    return ''

                # bool() needed here because they are of type "numpy._bool" by default
                return [all(bool(x) is True for x in junc_results)]

            return f


        tmp = OrderedDict()
        self.types2headers = {'psi':[], 'dpsi':[]}

        tmp['event_changing'] = (_global_changing, self.config.voila_files)
        tmp['event_non_changing'] = (_global_non_changing, self.config.voila_files)

        for voila_file in self.config.voila_files:

            with Matrix(voila_file) as m:
                analysis_type = m.analysis_type
                group_names = m.group_names
                if analysis_type == constants.ANALYSIS_HETEROGEN:
                    stat_names = m.stat_names
                else:
                    stat_names = None


            if analysis_type == constants.ANALYSIS_PSI:
                for group in group_names:
                    for key in ("E(PSI)", "Var(E(PSI))",):
                        header = "%s_%s" % (group, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                        else:
                            if key == "E(PSI)":
                                tmp[header] = (_psi_psi, [voila_file])
                            elif key == "Var(E(PSI))":
                                tmp[header] = (_psi_var, [voila_file])


            elif analysis_type == constants.ANALYSIS_HETEROGEN:

                group_idxs = {}
                for i, group in enumerate(group_names):
                    group_idxs[group] = i
                    for key in ("E(PSI)",):
                        header = "%s_HET_%s" % (group, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                        else:
                            if key == "E(PSI)":
                                tmp[header] = (_het_psi, [voila_file], i)



                for group1, group2 in combinations(group_names, 2):
                    for key in ("E(dPSI)",):
                        header = "%s-%s_HET_%s" % (group1, group2, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                        else:
                            if key == "E(dPSI)":
                                self.dpsi_quant_idxs.append(len(tmp))
                                tmp[header] = (_het_dpsi, [voila_file], group_idxs[group1], group_idxs[group2])

                    for j, key in enumerate(stat_names):
                        header = "%s-%s_HET_%s" % (group1, group2, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                        else:
                            tmp[header] = (_het_stats, [voila_file], j)


            else:
                for i, group in enumerate(group_names):
                    for key in ("E(PSI)",):
                        header = "%s_%s" % (group, key)
                        if header in tmp:
                            tmp[header][1].append(voila_file)
                        else:
                            if key == "E(PSI)":
                                tmp[header] = (_dpsi_psi, [voila_file], i)
                        self.types2headers['psi'].append(header)

                changing_thresh_key = "P(|dPSI|>=%.2f)" % self.config.changing_between_group_dpsi
                non_changing_thresh_key = "P(|dPSI|<=%.2f)" % self.config.non_changing_threshold
                for key in ("E(dPSI)", changing_thresh_key, non_changing_thresh_key):
                    header = "%s_%s" % ('-'.join(reversed(group_names)), key)
                    if header in tmp:
                        tmp[header][1].append(voila_file)
                    else:
                        if key == "E(dPSI)":
                            self.dpsi_quant_idxs.append(len(tmp))
                            tmp[header] = (_dpsi_dpsi, [voila_file])
                            self.types2headers['dpsi'].append(header)

                        elif key == changing_thresh_key:
                            tmp[header] = (_dpsi_p_change, [voila_file])
                        elif key == non_changing_thresh_key:
                            tmp[header] = (_dpsi_p_nonchange, [voila_file])

        # else:
        #     # in the case of training data, the output needs to be explicitly for each input file
        #     # we need a kind of format like {'input file name': {'psi1': psi1 function
        #     tmp = OrderedDict()
        #     for voila_file in self.config.voila_files:
        #         with Matrix(voila_file) as m:
        #             analysis_type = m.analysis_type
        #             group_names = m.group_names
        #
        #         for i, group in enumerate(group_names):
        #             for key in ("E(PSI)",):
        #                 header = "%s_%s" % (group, key)
        #                 if header in tmp:
        #                     tmp[header][1].append(voila_file)
        #                 else:
        #                     if key == "E(PSI)":
        #                         tmp[header] = (_dpsi_psi, [voila_file], i)
        #
        #         changing_thresh_key = "P(|dPSI|>=%.2f)" % self.config.changing_threshold
        #         non_changing_thresh_key = "P(|dPSI|<=%.2f)" % self.config.non_changing_threshold
        #         for key in ("E(dPSI)", changing_thresh_key, non_changing_thresh_key):
        #             header = "%s_%s" % ('-'.join(group_names), key)
        #             if header in tmp:
        #                 tmp[header][1].append(voila_file)
        #             else:
        #                 if key == "E(dPSI)":
        #                     tmp[header] = (_dpsi_dpsi, [voila_file])
        #                 elif key == changing_thresh_key:
        #                     tmp[header] = (_dpsi_p_change, [voila_file])
        #                 elif key == non_changing_thresh_key:
        #                     tmp[header] = (_dpsi_p_nonchange, [voila_file])

        return tmp