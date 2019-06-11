import csv
import os
from voila.api import Matrix
from voila import constants
from voila.exceptions import GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile
from voila.api.matrix_utils import generate_variances
from voila.api import view_matrix
from collections import OrderedDict
from voila.config import ClassifyConfig
import multiprocessing
import numpy as np
from voila.vlsv import get_expected_psi, matrix_area


def semicolon(value_list):
    return ';'.join(str(x) for x in value_list)


summaryVars2Headers = {
    'cassette_exon': 'Cassette',
    'tandem_cassette': 'Tandem Cassette',
    'alt3ss': 'Alt 3',
    'alt5ss': 'Alt 5',
    'p_alt3ss': 'P_Alt 3',
    'p_alt5ss': 'P_Alt 5',
    'alt3and5ss': 'Alt 3 and Alt 5',
    'mutually_exclusive': 'MXE',
    'alternative_intron': 'Intron Retention',
    'ale': 'ALE',
    'afe': 'AFE',
    'p_ale': 'P_ALE',
    'p_afe': 'P_AFE',
    'orphan_junction': 'Orphan Junction',
    'constitutive': 'Constitutive Junction',
    'constitutive_intron': 'Constitutive Intron',
    'multi_exon_spanning': 'Multi Exon Spanning',
    'exitron': 'Exitron',
}

class TsvWriter:
    """
    Output AS data from one gene
    """

    def __init__(self, graph, gene_id):
        """

        :param output_path: The folder where all output TSV files will be written under
        :param graph: the Graph object of the gene
        """

        self.common_headers = ['Module ID', 'Gene ID', 'Gene Name', 'Chr', 'Strand', 'LSV ID(s)']
        self.graph = graph
        self.gene_id = gene_id
        self.config = ClassifyConfig()
        self.quantifications_int = self.quantification_intersection()
        self.pid = multiprocessing.current_process().pid

        # we could do some crazy thing to yield to all of the different output types at once (across each method)
        # (in order to save memory) But for now we just save modules in a list. Will ammend later if memory use
        # becomes an issue.
        if self.graph:
            self.modules = self.graph.modules()

            # self.as_types = {x.idx: x.as_types() for x in self.modules}
            self.as_types = {x.idx: x.as_types() for x in self.modules}

    def quantification_intersection(self):
        """
        Look at all psi and dpsi quant headers and find the appropriate intersection
        we need to then define which function should be called for each resulting column
        :return:
        """
        SIG_FIGS = 3

        def _filter_edges(edge, lsv):
            if type(edge) != list:
                edge = [edge]
            for _edge in edge:
                # loop through junctions to find one matching range of edge
                for j, junc in enumerate(lsv.get('junctions')):
                    if junc[0] == _edge.start and junc[1] == _edge.end:
                        return j
                else:
                    # junction not quantified by majiq
                    pass

        def _psi_psi(voila_file):
            def f(lsv_id, edge=None):
                with Matrix(voila_file) as m:
                    lsv = m.psi(lsv_id)
                    if edge:
                        edge_idx = _filter_edges(edge, lsv)
                        if edge_idx is None:
                            return ''
                        else:
                            return round(lsv.get('means')[edge_idx], SIG_FIGS)
                    else:
                        return (round(x, SIG_FIGS) for x in lsv.get('means'))
            return f

        def _psi_var(voila_file):
            def f(lsv_id, edge=None):
                with Matrix(voila_file) as m:
                    lsv = m.psi(lsv_id)
                    if edge:
                        edge_idx = _filter_edges(edge, lsv)
                        if edge_idx is None:
                            return ''
                        else:
                            return round(generate_variances([lsv.get('bins')][0])[edge_idx], SIG_FIGS)
                    else:
                        return (round(x, SIG_FIGS) for x in generate_variances([lsv.get('bins')[0]]))
            return f

        def _dpsi_psi(voila_file, group_idx):
            def f(lsv_id, edge=None):
                with Matrix(voila_file) as m:
                    lsv = m.delta_psi(lsv_id)
                    if edge:
                        edge_idx = _filter_edges(edge, lsv)
                        if edge_idx is None:
                            return ''
                        else:
                            return round(lsv.get('group_means')[group_idx][edge_idx], SIG_FIGS)
                    else:
                        return (round(x, SIG_FIGS) for x in lsv.get('group_means')[group_idx])
            return f

        def _dpsi_dpsi(voila_file):
            def f(lsv_id, edge=None):
                with view_matrix.ViewDeltaPsi(voila_file) as m:
                    lsv = m.lsv(lsv_id)
                    bins = lsv.get('group_bins')
                    if edge:
                        edge_idx = _filter_edges(edge, lsv)
                        if edge_idx is None:
                            return ''
                        else:
                            return round(lsv.excl_incl[edge_idx][1] - lsv.excl_incl[edge_idx][0], SIG_FIGS)
                    else:
                        return (
                                    round(x, SIG_FIGS) for x in ((lsv.excl_incl[i][1] - lsv.excl_incl[i][0] for i in
                                    range(np.size(bins, 0))))
                                )
            return f

        def _dpsi_p_thresh(voila_file):
            def f(lsv_id, edge=None):
                with view_matrix.ViewDeltaPsi(voila_file) as m:
                    lsv = m.lsv(lsv_id)

                    bins = lsv.bins
                    if edge:
                        edge_idx = _filter_edges(edge, lsv)
                        if edge_idx is None:
                            return ''
                        else:
                            return round(matrix_area(bins[edge_idx], self.config.threshold), SIG_FIGS)
                    else:
                        return (
                                    round(matrix_area(b, self.config.threshold), SIG_FIGS) for b in bins
                                )
            return f

        def _dpsi_p_nonchange(voila_file):
            def f(lsv_id, edge=None):
                #lsv = m.delta_psi(lsv_id)
                with view_matrix.ViewDeltaPsi(voila_file) as m:
                    lsv = m.lsv(lsv_id)
                    if edge:
                        edge_idx = _filter_edges(edge, lsv)
                        if edge_idx is None:
                            return ''
                        else:
                            return round(lsv.high_probability_non_changing()[edge_idx], SIG_FIGS)
                    else:
                        return (round(x, SIG_FIGS) for x in lsv.high_probability_non_changing())

            return f

        tmp = OrderedDict()
        for voila_file in self.config.voila_files:

            with Matrix(voila_file) as m:
                analysis_type = m.analysis_type
                group_names = m.group_names


            if analysis_type == constants.ANALYSIS_PSI:
                for group in group_names:
                    for key in ("E(PSI)", "Var(E(PSI))",):
                        header = "%s_%s" % (group, key)
                        if not header in tmp:
                            if key == "E(PSI)":
                                tmp[header] = _psi_psi(voila_file)
                            elif key == "Var(E(PSI))":
                                tmp[header] = _psi_var(voila_file)
                        else:
                            pass
                            #print("found duplicate key %s" % header)
            else:
                for i, group in enumerate(group_names):
                    for key in ("E(PSI)",):
                        header = "%s_%s" % (group, key)
                        if not header in tmp:
                            if key == "E(PSI)":
                                tmp[header] = _dpsi_psi(voila_file, i)
                        else:
                            pass
                            #print("found duplicate key %s" % header)

                thresh_key = "P(|dPSI|>=%.2f)" % self.config.threshold
                high_prob_thresh_key = "P(|dPSI|<=%.2f)" % self.config.non_changing_threshold
                for key in ("E(dPSI)", thresh_key, high_prob_thresh_key):
                    header = "%s_%s" % ('-'.join(group_names), key)
                    if not header in tmp:
                        if key == "E(dPSI)":
                            tmp[header] = _dpsi_dpsi(voila_file)
                        elif key == thresh_key:
                            tmp[header] = _dpsi_p_thresh(voila_file)
                        elif key == high_prob_thresh_key:
                            tmp[header] = _dpsi_p_nonchange(voila_file)
                    else:
                        pass
                        #print("found duplicate key %s" % header)

                # tmp.append('%s_E(dPSI) per LSV junction' % trunc_name)
                # tmp.append('%s_E(dPSI) per LSV junction' % trunc_name)
                # tmp.append('%s_E(dPSI) per LSV junction' % trunc_name)
                # tmp.append('%s_E(dPSI) per LSV junction' % trunc_name)
                # tmp.append('%s_E(dPSI) per LSV junction' % trunc_name)

        return tmp

    @property
    def quantification_headers(self):
        return list(self.quantifications_int.keys())


    @staticmethod
    def tsv_names():
        if ClassifyConfig().putative_multi_gene_regions:
            return ['p_multi_gene_region.tsv']
        names = ['summary.tsv', 'cassette.tsv', 'alt3prime.tsv', 'alt5prime.tsv', 'alt3and5prime.tsv',
                 'mutually_exclusive.tsv', 'alternate_last_exon.tsv', 'alternate_first_exon.tsv',
                 'alternative_intron.tsv', 'p_alt5prime.tsv', 'p_alt3prime.tsv', 'multi_exon_spanning.tsv',
                 'tandem_cassette.tsv', 'exitron.tsv', 'p_alternate_last_exon.tsv', 'p_alternate_first_exon.tsv']
        if ClassifyConfig().keep_constitutive:
            names.append('constitutive.tsv')
        return names

    @staticmethod
    def delete_tsvs():
        for tsv_file in TsvWriter.tsv_names():
            config = ClassifyConfig()
            path = os.path.join(config.directory, tsv_file)
            if os.path.exists(path):
                os.remove(path)

    def parity2lsv(self, module, parity, edge=None, node=None):

        if parity == 's':
            lsvs = module.source_lsv_ids
            if edge or node:
                if not node:
                    node = self.graph.start_node(edge)
                lsvs = set(filter(lambda lsv: lsv.endswith(node.untrimmed_range_str()), lsvs))

        elif parity == 't':
            lsvs = module.target_lsv_ids
            if edge or node:
                if not node:
                    node = self.graph.end_node(edge)
                lsvs = set(filter(lambda lsv: lsv.endswith(node.untrimmed_range_str()), lsvs))
        else:
            lsvs = module.target_lsv_ids.union(module.source_lsv_ids)
        return lsvs

    def common_data(self, module, parity=None, edge=None, node=None):
        """
        Extract the certain cols from the CSV which are generally similar across all outputs,

        """
        lsvs = self.parity2lsv(module, parity, edge, node)

        return ["%s_%d" % (self.gene_id, module.idx), self.gene_id, self.graph.gene_name,
                self.graph.chromosome, self.graph.strand, semicolon(lsvs)]

    def quantifications(self, module, parity=None, edge=None, node=None):
        """
        Edge / Parity is used to find LSVs
        Node is used to filter lsvs to specific node (the node that has THAT lsv)
        :return:
        """

        lsvs = self.parity2lsv(module, parity, node=node)

        out = []

        for field in self.quantifications_int:
            quantification_vals = []
            for lsv_id in lsvs:
                try:

                    #print(self.quantifications_int[field](lsv_id, edge))
                    quantification_vals.append(self.quantifications_int[field](lsv_id, edge))
                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile) as e:
                    quantification_vals.append('')
                    #print(e)
            out.append(semicolon(quantification_vals))


        return out


    def start_headers(self, headers, filename):
        """
        Start a tsv file with the required headers, only if it does not yet exist

        """
        if not os.path.exists(os.path.join(self.config.directory, filename)):
            with open(os.path.join(self.config.directory, filename), 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
                writer.writerow(headers)

    def start_all_headers(self):

        if self.config.putative_multi_gene_regions:
            headers = ['Gene ID_Region', 'Gene ID', 'Gene Name', 'Chr', 'Strand', 'First Exon Start coord',
                       'First Exon End coord', 'Last Exon Start coord', "Last Exon End coord"]
            self.start_headers(headers, 'p_multi_gene_region.tsv')
            return

        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Junction Name',
                                         'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'cassette.tsv')
        self.start_headers(headers, 'alt3prime.tsv')
        self.start_headers(headers, 'alt5prime.tsv')
        self.start_headers(headers, 'p_alt5prime.tsv')
        self.start_headers(headers, 'p_alt3prime.tsv')
        self.start_headers(headers, 'alt3and5prime.tsv')
        self.start_headers(headers, 'mutually_exclusive.tsv')
        self.start_headers(headers, 'alternate_last_exon.tsv')
        self.start_headers(headers, 'alternate_first_exon.tsv')
        self.start_headers(headers, 'p_alternate_last_exon.tsv')
        self.start_headers(headers, 'p_alternate_first_exon.tsv')
        self.start_headers(headers, 'alternative_intron.tsv')
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Tandem Exon Coordinates', 'Num_Tandem_Exons',
                                         'Junction Name', 'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'multi_exon_spanning.tsv')
        self.start_headers(headers, 'tandem_cassette.tsv')
        headers = self.common_headers + ['Exon coordinate', 'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'exitron.tsv')
        headers = self.common_headers + ["Cassette", "Tandem Cassette",
                                         "Alt 3", "Alt 5", "P_Alt 3", "P_Alt 5", "Alt 3 and Alt 5", "MXE",
                                         "Alternative Intron", "ALE", "AFE",
                                         "P_ALE", "P_AFE", "Orphan Junction"]
        if self.config.keep_constitutive:
            headers.append("Constitutive Junction")
            headers.append("Constitutive Intron")
        headers += ["Multi Exon Spanning", "Exitron", "Complex", "Number of Events", "Collapsed Event Name"]
        self.start_headers(headers, 'summary.tsv')
        if self.config.keep_constitutive:
            headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                             'Exon Spliced With Coordinate', 'Junction Name',
                                             'Junction Coordinate', 'Is Intron',
                                             'Collapsed Event Name'] + self.quantification_headers
            self.start_headers(headers, 'constitutive.tsv')

    def cassette(self):
        with open(os.path.join(self.config.directory, 'cassette.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'cassette_exon':
                            src_common = self.common_data(module, 's', event['C1'])
                            trg_common = self.common_data(module, 't', event['C2'])

                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2',
                                   event['Skip'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Skip'], event['C1']))

                            row = [event['C1'].range_str(), 'A', event['A'].range_str(), 'C1_A',
                                   event['Include1'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1'], event['C1']))

                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1',
                                   event['Skip'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Skip'], event['C2']))

                            row = [event['C2'].range_str(), 'A', event['A'].range_str(), 'C2_A',
                                   event['Include2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2'], event['C2']))

    def alt3prime(self):
        with open(os.path.join(self.config.directory, 'alt3prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'alt3ss':
                            src_common = self.common_data(module, 's', node=event['E1'])
                            trg_common = self.common_data(module, 't', node=event['E2'])

                            if src_common[5]:
                                row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
                                       event['Proximal'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Proximal']))
                                row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
                                       event['Distal'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Distal']))
                            if trg_common[5]:
                                row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Proximal',
                                       event['Proximal'].range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', event['Proximal']))
                                row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Distal',
                                       event['Distal'].range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', event['Distal']))

    def alt5prime(self):
        with open(os.path.join(self.config.directory, 'alt5prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'alt5ss':
                            src_common = self.common_data(module, 's', node=event['E1'])
                            trg_common = self.common_data(module, 't', node=event['E2'])
                            if src_common[5]:
                                row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
                                       event['Proximal'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Proximal']))
                                row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
                                       event['Distal'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Distal']))
                            if trg_common[5]:
                                row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Proximal',
                                       event['Proximal'].range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', event['Proximal']))
                                row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Distal',
                                       event['Distal'].range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', event['Distal']))

    def p_alt5prime(self):
        with open(os.path.join(self.config.directory, 'p_alt5prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'p_alt5ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')

                            row = [event['C1'].range_str(), 'E2', event['A'].range_str(), 'E3_E2_Proximal',
                                   event['Include2'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include2']))

                            row = [event['C1'].range_str(), 'E1', event['C2'].range_str(), 'E3_E1_Distal',
                                   event['Skip'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Skip']))

                            event['Include1'].junc['start'] += 1
                            event['Include1'].junc['end'] -= 1

                            row = [event['C2'].range_str(), 'E2', event['A'].range_str(), 'E1_E2_Intron',
                                   event['Include1'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include1']))


    def p_alt3prime(self):
        with open(os.path.join(self.config.directory, 'p_alt3prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'p_alt3ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')

                            row = [event['C2'].range_str(), 'E2', event['A'].range_str(), 'E1_E2_Proximal',
                                   event['Include1'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1']))

                            row = [event['C2'].range_str(), 'E3', event['C1'].range_str(), 'E1_E3_Distal',
                                   event['Skip'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Skip']))

                            event['Include2'].junc['start'] += 1
                            event['Include2'].junc['end'] -= 1

                            row = [event['C1'].range_str(), 'E2', event['A'].range_str(), 'E3_E2_Intron',
                                   event['Include2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2']))

    def alt3and5prime(self):
        with open(os.path.join(self.config.directory, 'alt3and5prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'alt3and5ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_J1',
                                   event['J1'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['J1']))
                            row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_J2',
                                   event['J2'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['J2']))
                            row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_J1',
                                   event['J1'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['J1']))
                            row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_J2',
                                   event['J2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['J2']))

    def mutually_exclusive(self):
        with open(os.path.join(self.config.directory, 'mutually_exclusive.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'mutually_exclusive':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['C1'].range_str(), 'A1', event['A1'].range_str(), 'C1_A1',
                                   event['Include1'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1']))
                            row = [event['C1'].range_str(), 'A2', event['A2'].range_str(), 'C1_A2',
                                   event['SkipA1'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['SkipA1']))
                            row = [event['C2'].range_str(), 'A1', event['A1'].range_str(), 'C2_A1',
                                   event['Include2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2']))
                            row = [event['C2'].range_str(), 'A2', event['A2'].range_str(), 'C2_A2',
                                   event['SkipA2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['SkipA2']))

    def alternate_last_exon(self):
        with open(os.path.join(self.config.directory, 'alternate_last_exon.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'ale':
                            src_common = self.common_data(module, 's', node=event['Reference'])
                            if src_common[5]:
                                for junc in event['SkipA2']:
                                    row = [event['Reference'].range_str(), 'A', event['Proximal'].range_str(),
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    writer.writerow(src_common + row + self.quantifications(module, 's', junc))
                                for junc in event['SkipA1']:
                                    row = [event['Reference'].range_str(), 'A', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    writer.writerow(src_common + row + self.quantifications(module, 's', junc))

    def alternate_first_exon(self):
        with open(os.path.join(self.config.directory, 'alternate_first_exon.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'afe':
                            trg_common = self.common_data(module, 't', node=event['Reference'])
                            if trg_common[5]:
                                for junc in event['SkipA1']:
                                    row = [event['Reference'].range_str(), 'A', event['Proximal'].range_str(),
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                                for junc in event['SkipA2']:
                                    row = [event['Reference'].range_str(), 'A', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))

    def p_alternate_last_exon(self):
        with open(os.path.join(self.config.directory, 'p_alternate_last_exon.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'p_ale':
                            src_common = self.common_data(module, 's', node=event['Reference'])
                            if event['Proximal'].start == -1:
                                proxStr = "nan-{}".format(event['Proximal'].end)
                            else:
                                proxStr = "{}-nan".format(event['Proximal'].start)
                            if src_common[5]:
                                for junc in event['SkipA2']:
                                    row = [event['Reference'].range_str(), 'A', proxStr,
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    writer.writerow(src_common + row + self.quantifications(module, 's', junc))
                                for junc in event['SkipA1']:
                                    row = [event['Reference'].range_str(), 'A', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    writer.writerow(src_common + row + self.quantifications(module, 's', junc))

    def p_alternate_first_exon(self):
        with open(os.path.join(self.config.directory, 'p_alternate_first_exon.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'p_afe':
                            if event['Proximal'].start == -1:
                                proxStr = "nan-{}".format(event['Proximal'].end)
                            else:
                                proxStr = "{}-nan".format(event['Proximal'].start)
                            trg_common = self.common_data(module, 't', node=event['Reference'])
                            if trg_common[5]:
                                for junc in event['SkipA1']:
                                    row = [event['Reference'].range_str(), 'A', proxStr,
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                                for junc in event['SkipA2']:
                                    row = [event['Reference'].range_str(), 'A', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))

    def alternative_intron(self):
        with open(os.path.join(self.config.directory, 'alternative_intron.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'alternative_intron':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')

                            # put coordinates back to Jordi's offset numbers
                            event['Intron'].junc['start'] += 1
                            event['Intron'].junc['end'] -= 1
                            if any(':t:' in _l for _l in event['Intron'].lsvs) and not \
                               any(':s:' in _l for _l in event['Intron'].lsvs):
                                row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1_intron',
                                       event['Intron'].range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', event['Intron']))
                                row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1_spliced',
                                       semicolon((x.range_str() for x in event['Spliced']))]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', event['Spliced']))
                            else:
                                row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2_intron',
                                       event['Intron'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Intron']))
                                row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2_spliced',
                                       semicolon((x.range_str() for x in event['Spliced']))]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Spliced']))

    def multi_exon_spanning(self):
        with open(os.path.join(self.config.directory, 'multi_exon_spanning.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'multi_exon_spanning':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C1_C2',
                                   semicolon((x.range_str() for x in event['Skip']))]
                            writer.writerow(src_common + row + self.quantifications(module, 's'))
                            row = [event['C1'].range_str(), 'A1', event['As'][0].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C1_A',
                                   semicolon((x.range_str() for x in event['Include1']))]
                            writer.writerow(src_common + row + self.quantifications(module, 's')), event['C1']
                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C2_C1',
                                   semicolon((x.range_str() for x in event['Skip']))]
                            writer.writerow(trg_common + row + self.quantifications(module, 't'))
                            row = [event['C2'].range_str(), 'A_Last', '',
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'A_Last_C2',
                                   semicolon((x.range_str() for x in event['Include2']))]
                            writer.writerow(trg_common + row + self.quantifications(module, 't'))

    def tandem_cassette(self):
        with open(os.path.join(self.config.directory, 'tandem_cassette.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'tandem_cassette':
                            src_common = self.common_data(module, 's', node=event['C1'])
                            trg_common = self.common_data(module, 't', node=event['C2'])
                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C1_C2',
                                   semicolon((x.range_str() for x in event['Skip']))]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Skip'][0], event['C1']))
                            row = [event['C1'].range_str(), 'A1', event['As'][0].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C1_A',
                                   semicolon((x.range_str() for x in event['Include1']))]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1'][0], event['C1']))
                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C2_C1',
                                   semicolon((x.range_str() for x in event['Skip']))]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Skip'][0], event['C2']))
                            row = [event['C2'].range_str(), 'A_Last', event['As'][-1].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'A_Last_C2',
                                   semicolon((x.range_str() for x in event['Include2']))]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2'][0], event['C2']))

    def exitron(self):
        with open(os.path.join(self.config.directory, 'exitron.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'exitron':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['Exon'].range_str(), event['Junc'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's'))
                            row = [event['Exon'].range_str(), event['Junc'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't'))

    def constitutive(self):
        with open(os.path.join(self.config.directory, 'constitutive.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'constitutive':
                            common = self.common_data(module)
                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1',
                                   event['Junc'].range_str(), 'False', module.collapsed_event_name]
                            writer.writerow(common + row + self.quantifications(module, edge=event['Junc']))

                        elif event['event'] == 'constitutive_intron':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')

                            # put coordinates back to Jordi's offset numbers
                            event['Intron'].junc['start'] += 1
                            event['Intron'].junc['end'] -= 1
                            if any(':t:' in _l for _l in event['Intron'].lsvs) and not \
                               any(':s:' in _l for _l in event['Intron'].lsvs):
                                row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1_intron',
                                       event['Intron'].range_str(), 'True', module.collapsed_event_name]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', event['Intron']))
                            else:
                                row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2_intron',
                                       event['Intron'].range_str(), 'True', module.collapsed_event_name]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Intron']))


    def p_multi_gene_region(self):
        with open(os.path.join(self.config.directory, 'p_multi_gene_region.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                for event in events:
                    if event['event'] == 'p_multi_gene_region':
                        row = ["%s_Region%d" % (self.gene_id, event['idx']), self.gene_id, self.graph.gene_name,
                               self.graph.chromosome, self.graph.strand, event['ExonStart'].start,
                               event['ExonStart'].end, event['ExonEnd'].start,
                               event['ExonEnd'].end]
                        writer.writerow(row)

    def summary(self):
        """
        Write the summary style output file
        :param genes_modules: a list of (gene_id (str), gene_modules (obj)) tuples
        :return: NOTHING
        """
        with open(os.path.join(self.config.directory, 'summary.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                counts = OrderedDict()
                counts['cassette_exon'] = 0
                counts['tandem_cassette'] = 0
                counts['alt3ss'] = 0
                counts['alt5ss'] = 0
                counts['p_alt3ss'] = 0
                counts['p_alt5ss'] = 0
                counts['alt3and5ss'] = 0
                counts['mutually_exclusive'] = 0
                counts['alternative_intron'] = 0
                counts['ale'] = 0
                counts['afe'] = 0
                counts['p_ale'] = 0
                counts['p_afe'] = 0
                counts['orphan_junction'] = 0
                if self.config.keep_constitutive:
                    counts['constitutive'] = 0
                    counts['constitutive_intron'] = 0
                counts['multi_exon_spanning'] = 0
                counts['exitron'] = 0
                for event in events:
                    if event['event'] in counts:
                        counts[event['event']] += 1

                # we store collapsed event name on module, because we need it for constitutive
                module.collapsed_event_name = self._collapsed_event_name(counts)

                writer.writerow(["%s_%d" % (self.gene_id, module.idx),
                                 self.gene_id, self.graph.gene_name, self.graph.chromosome, self.graph.strand,
                                 semicolon(module.target_lsv_ids.union(module.source_lsv_ids))] +
                                [v if v else '' for v in counts.values()] + [str(_complex), str(_total_events),
                                                                             module.collapsed_event_name]
                                )

    def _collapsed_event_name(self, counts):
        """
        function to generate the collapsed event name from event counts
        """
        out = []
        for count in counts:
            if counts[count]:
                out.append("%sx%d" % (summaryVars2Headers[count], counts[count]))
        return '_'.join(out)

