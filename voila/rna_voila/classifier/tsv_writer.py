import csv
import os
from rna_voila.api import Matrix
from rna_voila import constants
from rna_voila.exceptions import GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile
from rna_voila.api.matrix_utils import generate_variances
from rna_voila.api import view_matrix
from collections import OrderedDict
from rna_voila.config import ClassifyConfig
import multiprocessing
from rna_voila.classifier.quantification_finder import MultiQuantWriter
import numpy as np
from rna_voila.vlsv import get_expected_psi, matrix_area
from itertools import combinations
from operator import itemgetter
import sys, json
import networkx as nx

summaryVars2Headers = {
    'cassette_exon': 'cassette',
    'tandem_cassette': 'tandem_cassette',
    'alt3ss': 'alt3',
    'alt5ss': 'alt5',
    'p_alt3ss': 'putative_alt3',
    'p_alt5ss': 'putative_alt5',
    'alt3and5ss': 'alt3_5',
    'mutually_exclusive': 'mxe',
    'alternative_intron': 'ir',
    'ale': 'ale',
    'afe': 'afe',
    'p_ale': 'putative_ale',
    'p_afe': 'putative_afe',
    'orphan_junction': 'orphan_junc',
    'constitutive': 'constitutive_junc',
    'constitutive_intron': 'persistent_ir',
    'multi_exon_spanning': 'multi_exon_spanning',
    'other_event': 'other'
}

reverse_strand_map = {
    'start': 'end',
    'end': 'start',
    'absolute_start': 'absolute_end',
    'absolute_end': 'absolute_start'
}

class BaseTsvWriter(MultiQuantWriter):
    """
    Output AS data from one gene
    """
    def __init__(self, graph, gene_id):
        """
        :param output_path: The folder where all output TSV files will be written under
        :param graph: the Graph object of the gene
        """
        self.graph = graph
        super().__init__()

        self.common_headers = ['module_id', 'gene_id', 'gene_name', 'seqid', 'strand']


        self.gene_id = gene_id
        self.strand = self.graph.strand if self.graph else None


        self.pid = multiprocessing.current_process().pid

        self.heatmap_cache = {}
        self.junction_cache = []
        self.training_junction_cache = {}

        self.individual_nonchange_cols = []
        self.individual_change_cols = []
        if self.config.include_change_cases:
            for voila_file in self.config.voila_files:
                with Matrix(voila_file) as m:
                    analysis_type = m.analysis_type
                    if analysis_type in (constants.ANALYSIS_DELTAPSI, constants.ANALYSIS_HETEROGEN):
                        self.individual_nonchange_cols.append(f"{m.group_names[0]}_{m.group_names[1]}_non_changing")
                        self.individual_change_cols.append(f"{m.group_names[0]}_{m.group_names[1]}_changing")


    @property
    def quantification_headers(self):
        return ['event_non_changing', 'event_changing'] + self.individual_nonchange_cols + self.individual_change_cols + list(self.quantifications_int.keys())

    @staticmethod
    def tsv_names():
        return []

    @classmethod
    def delete_tsvs(cls):
        for tsv_file in cls.tsv_names():
            config = ClassifyConfig()
            path = os.path.join(config.directory, tsv_file)
            if os.path.exists(path):
                os.remove(path)

    def common_data(self, module, parity=None, edge=None, node=None, event_name=None, event_ii=None):
        """
        Extract the certain cols from the TSV which are generally similar across all outputs,

        """

        if not edge and node is not None:
            lsvs = self.parity2lsv_node(module, parity, node)
        else:
            lsvs = self.parity2lsv(module, parity, edge)

        module_id = "%s_%d" % (self.gene_id, module.idx)
        if event_ii is None:
            event_id = "%s_%s" % (module_id, event_name)
        else:
            event_id = "%s_%s_%s" % (module_id, event_name, event_ii)

        out = [module_id,
               self.gene_id,
               self.graph.gene_name,
               self.graph.chromosome,
               self.graph.strand]

        out.append(self.semicolon(lsvs))
        out.append(event_id)
        out.append(str(module.is_complex))
        return out

    def start_headers(self, headers, filename):
        """
        Start a tsv file with the required headers, only if it does not yet exist

        """
        _has_dpsi_voila_files = False
        _has_het_voila_files = False
        for voila_file in self.config.voila_files:
            with Matrix(voila_file) as m:
                if m.analysis_type == constants.ANALYSIS_DELTAPSI:
                    _has_dpsi_voila_files = True
                if m.analysis_type == constants.ANALYSIS_HETEROGEN:
                    _has_het_voila_files = True

        if not os.path.exists(os.path.join(self.config.directory, filename)):
            with open(os.path.join(self.config.directory, filename), 'w', newline='') as csvfile:
                if not self.config.disable_metadata:
                    metadata = {
                        'voila_version': constants.VERSION,
                        'command': ' '.join(sys.argv),
                    }
                    if _has_dpsi_voila_files:
                        metadata['dpsi_changing_threshold'] = self.config.changing_between_group_dpsi
                        metadata['dpsi_probability_changing_threshold'] = self.config.probability_changing_threshold
                        metadata['dpsi_nonchanging_threshold'] = self.config.probability_non_changing_threshold

                    if _has_het_voila_files:
                        metadata['het_changing_threshold'] = self.config.changing_between_group_dpsi
                        metadata['het_pvalue_changing_threshold'] = self.config.changing_pvalue_threshold
                        metadata['het_nonchanging_threshold'] = self.config.non_changing_between_group_dpsi
                        metadata['het_pvalue_nonchanging_threshold'] = self.config.non_changing_pvalue_threshold
                        metadata['het_IQR_nonchanging_threshold'] = self.config.non_changing_within_group_iqr

                    metadata_string = '\n'.join([f"# {l}" for l in json.dumps(metadata, sort_keys=True,
                                                                              indent=4, separators=(',', ': ')).split(
                        '\n')]) + '\n'
                    csvfile.write(metadata_string)

                writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
                writer.writerow(headers)

    def _collapsed_event_name(self, counts, other_novel_count=0):
        """
        function to generate the collapsed event name from event counts
        """
        out = []
        for count in counts:
            # do not include 'other' in the collapsed event name.
            if count == "other_event":
                continue
            if counts[count]:
                out.append("%s^%d" % (summaryVars2Headers[count], counts[count]))
        if other_novel_count > 0:
            out.append("%s^%d" % (summaryVars2Headers['other_event'], other_novel_count))
        return '|'.join(out)


class TsvWriter(BaseTsvWriter):

    def __init__(self, graph, gene_id):
        super().__init__(graph, gene_id)

        self.common_headers.append('lsv_id')
        self.common_headers.append('event_id')
        self.common_headers.append('complex')

        # we could do some crazy thing to yield to all of the different output types at once (across each method)
        # (in order to save memory) But for now we just save modules in a list. Will ammend later if memory use
        # becomes an issue.
        if self.graph:
            self.modules = self.graph.modules()

            # self.as_types = {x.idx: x.as_types() for x in self.modules}
            self.as_types = {x.idx: x.as_types() for x in self.modules}
            self.mpe_regions = {x.idx: x.mpe_regions() for x in self.modules}


    @staticmethod
    def tsv_names():
        config = ClassifyConfig()
        if config.putative_multi_gene_regions:
            return ['p_multi_gene_region.tsv']
        names = []
        if 'events' in config.enabled_outputs:
            names += ['cassette.tsv', 'alt3prime.tsv', 'alt5prime.tsv', 'alt3and5prime.tsv',
                 'mutually_exclusive.tsv', 'alternate_last_exon.tsv', 'alternate_first_exon.tsv',
                 'alternative_intron.tsv', 'p_alt5prime.tsv', 'p_alt3prime.tsv', 'multi_exon_spanning.tsv',
                 'tandem_cassette.tsv', 'orphan_junction.tsv',
                 'p_alternate_last_exon.tsv', 'p_alternate_first_exon.tsv', 'other.tsv']
            if config.keep_constitutive:
                names.append('constitutive.tsv')
        if 'heatmap' in config.enabled_outputs:
            names += ['heatmap.tsv']
        if 'junctions' in config.enabled_outputs:
            names += ['junctions.tsv']
        if 'summary' in config.enabled_outputs:
            names += ['summary.tsv']
        if 'mpe' in config.enabled_outputs:
            names += ['mpe_primerable_regions.tsv']
        if 'training_junctions' in config.enabled_outputs:
            names += ['training_junctions.tsv']
        if 'training_paths' in config.enabled_outputs:
            names += ['training_paths.tsv']

        return names


    def start_all_headers(self):

        if self.config.putative_multi_gene_regions:
            headers = self.common_headers + ['Exon1Start', 'Exon1End', 'Exon2Start', 'Exon2End']
            self.start_headers(headers, 'p_multi_gene_region.tsv')

        else:
            relatively_common_headers = ['denovo', 'reference_exon_coord', 'spliced_with',
                                             'spliced_with_coord', 'junction_name',
                                             'junction_coord']
            if 'events' in self.config.enabled_outputs:

                headers = self.common_headers + relatively_common_headers + ['event_size'] + self.quantification_headers
                if self.config.cassettes_constitutive_column:
                    headers.append('event_constitutive')
                self.start_headers(headers, 'cassette.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + ['event_size'] + self.quantification_headers, 'alt3prime.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + ['event_size'] + self.quantification_headers, 'alt5prime.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + ['event_size'] + self.quantification_headers, 'p_alt5prime.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + ['event_size'] + self.quantification_headers, 'p_alt3prime.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + self.quantification_headers, 'alt3and5prime.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + self.quantification_headers, 'mutually_exclusive.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + self.quantification_headers, 'alternate_last_exon.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + self.quantification_headers, 'alternate_first_exon.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + self.quantification_headers, 'p_alternate_last_exon.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + self.quantification_headers, 'p_alternate_first_exon.tsv')
                self.start_headers(self.common_headers + relatively_common_headers + ['event_size'] + self.quantification_headers, 'alternative_intron.tsv')
                headers = self.common_headers + ['junction_name',
                                                 'junction_coord',
                                                 'denovo',
                                                 'intron',
                                                 'reference_exon_coord',
                                                 'spliced_with_coord',
                                                 'exons_skipped_coords',
                                                 'num_skipped_exons'
                                                 ] + self.quantification_headers
                self.start_headers(headers, 'multi_exon_spanning.tsv')
                headers = self.common_headers + ['denovo', 'reference_exon_coord', 'spliced_with',
                                                 'spliced_with_coord', 'exons_skipped_coords',
                                                 'num_skipped_exons',
                                                 'junction_name', 'junction_coord'] + ['event_size'] + self.quantification_headers
                self.start_headers(headers, 'tandem_cassette.tsv')
                headers = self.common_headers + ['denovo', 'exon1_coord', 'Exon2 coordinate', 'junction_coord'] + self.quantification_headers
                self.start_headers(headers, 'orphan_junction.tsv')
                headers = self.common_headers + ["other_junctions", "other_exons"]# + self.quantification_headers
                self.start_headers(headers, 'other.tsv')

                if self.config.keep_constitutive:
                    headers = self.common_headers + ['denovo', 'reference_exon_coord', 'spliced_with',
                                                     'spliced_with_coord', 'junction_name',
                                                     'junction_coord', 'Is Intron',
                                                     'module_event_combination'] + self.quantification_headers
                    self.start_headers(headers, 'constitutive.tsv')





            if 'summary' in self.config.enabled_outputs:
                headers = self.common_headers + ['cassette', 'tandem_cassette',
                                                 'alt3', 'alt5', 'putative_alt3', 'putative_alt5', 'alt3_5', 'mxe',
                                                 'ir', 'ale', 'afe',
                                                 'putative_ale', 'putative_afe', 'orphan_junc', 'other']
                if self.config.keep_constitutive:
                    headers.append('constitutive_junc')
                    headers.append('persistent_ir')

                headers.remove('complex')
                event_id_ii = headers.index("event_id")
                headers.pop(event_id_ii)

                headers += ['multi_exon_spanning', 'complex', 'denovo_juncs', 'denovo_introns',
                            'num_events', 'module_event_combination'
                            ]

                self.start_headers(headers, 'summary.tsv')

                headers.insert(event_id_ii, "event_id")

            if 'heatmap' in self.config.enabled_outputs:
                headers = self.common_headers + ['module_event_combination', 'denovo', 'junction_name',
                                                 'junction_coord'] + self.quantification_headers
                self.start_headers(headers, 'heatmap.tsv')

            if 'junctions' in self.config.enabled_outputs:
                headers = self.common_headers + ['module_event_combination', 'denovo', 'junction_name',
                                                 'junction_coord'] + self.quantification_headers
                self.start_headers(headers, 'junctions.tsv')

            if 'mpe' in self.config.enabled_outputs:
                headers = self.common_headers
                lids_col_ii = headers.index("lsv_id")
                # Only 1 LSV ID per row possible.
                headers[lids_col_ii].replace("(s)","")
                headers += ['module_event_combination',
                        'Type',
                        'Edge of the Module',
                        'Edge of Transcript',
                        'Reference Exon Coord',
                        'Reference Exon De Novo',
                        'Reference Exon Exitrons',
                        'Reference Exon Constant Region',
                        'Reference Exon Trimmed',
                        'Constitutive Direction',
                        'Constitutive Regions',
                        'Constitutive De Novo',
                        'Constitutive Exon or Intron']
                self.start_headers(headers, 'mpe_primerable_regions.tsv')

            if 'training_junctions' in self.config.enabled_outputs:
                headers = self.common_headers + ['junc_id', 'junc_start_coordinate', 'junc_end_coordinate',
                                                 'exon_1_id', 'exon_1_start_coordinate', 'exon_1_end_coordinate',
                                                 'exon_2_id', 'exon_2_start_coordinate', 'exon_2_end_coordinate',
                                                 'is_intron'] + self.quantification_headers
                self.start_headers(headers, 'training_junctions.tsv')

            if 'training_paths' in self.config.enabled_outputs:
                headers = ['gene_id', 'seqid', 'strand', 'path_id', 'path', 'edge_type', 'junc_ids', 'exon_length_mod_3']
                self.start_headers(headers, 'training_paths.tsv')



    def _determine_changing(self, all_quants):

        s = set(x[0] for x in all_quants)
        if 'True' in s:
            return 'True'
        elif 'False' in s:
            return 'False'
        return ''

    def _determine_non_changing(self, all_quants):

        s = set(x[1] for x in all_quants)
        if 'False' in s:
            return 'False'
        elif 'True' in s:
            return 'True'
        return ''

    def _trim_strand_case_range_str(self, e1, pos1, e2, pos2):
        """
        produce a range_str based on two nodes start or ends, and a strand
        :param e1: exon / junction to use for begin
        :param pos1: attribute (start or end) string to use for begin
        :param e2: exon / junction to use for end
        :param pos2: attribute (start or end) string to use for end
        :param plus_strand: true is plus strand, false if minus strand
        :return:
        """
        if self.strand == '+':
            return f"{getattr(e1, pos1)}-{getattr(e2, pos2)}"
        return f"{getattr(e2,  reverse_strand_map[pos2])}-{getattr(e1, reverse_strand_map[pos1])}"

    def _trim_strand_case_event_size(self, e1, pos1, e2, pos2):
        """
        Same as _trim_strand_case_range_str, but produce a numeric result of the resulting absolute difference
        Between the start and end of the specified nodes / positions
        """

        if self.strand == '+':
            return abs(getattr(e1, pos1)-getattr(e2, pos2))
        return abs(getattr(e2, reverse_strand_map[pos2])-getattr(e1, reverse_strand_map[pos1]))

    def _add_training_junction(self, common, edge, quants, module_idx):
        if 'training_junctions' in self.config.enabled_outputs or 'training_paths' in self.config.enabled_outputs:
            if not module_idx in self.training_junction_cache:
                self.training_junction_cache[module_idx] = []
            self.training_junction_cache[module_idx].append((common, edge, quants,))

    def training_paths(self):
        with open(os.path.join(self.config.directory, 'training_paths.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module_idx, junc_list in self.training_junction_cache.items():
                # for each module, build a subgraph based only on junctions from events in the module
                g = nx.DiGraph()
                for common, edge, quants in junc_list:
                    e1 = self.graph.start_node(edge)
                    e2 = self.graph.end_node(edge)

                    g.add_edge(edge.absolute_start, edge.absolute_end, ir=edge.ir, _type='j')
                    g.add_edge(e1.absolute_start, e1.absolute_end, _type='e')
                    g.add_edge(e2.absolute_start, e2.absolute_end, _type='e')


                sink_nodes = [node for node, outdegree in g.out_degree(g.nodes()) if outdegree == 0]
                source_nodes = [node for node, indegree in g.in_degree(g.nodes()) if indegree == 0]


                path_idx = 0
                for source in source_nodes:
                    for sink in sink_nodes:
                        for path in nx.all_simple_paths(g, source=source, target=sink):
                            path_idx += 1
                            # one row for one path

                            path_id = f"{self.gene_id}_{self.graph.chromosome}_{self.graph.strand}_{module_idx}_{path_idx}"
                            junc_coords = []
                            junc_types = []
                            junc_ids = []
                            exon_lengths = []

                            for n1, n2 in zip(path, path[1:]):
                                data = g.get_edge_data(n1, n2)
                                if data['_type'] == 'e':
                                    exon_lengths.append(n2-n1)
                                    continue
                                _type = 'I' if data['ir'] else 'J'
                                junc_coords.append(f"{n1}-{n2}")
                                junc_types.append(_type)
                                junc_ids.append(f"{self.gene_id}_{self.graph.chromosome}_{self.graph.strand}_{n1}_{n2}_{_type}")

                            if len(exon_lengths) <= 2:
                                exon_mod = 'N/A'
                            else:
                                exon_mod = str(sum(exon_lengths[1:-1]) % 3)

                            path_out = "_".join(junc_coords)
                            edge_type_out = "_".join(junc_types)
                            junc_ids_out = "_".join(junc_ids)

                            writer.writerow([self.gene_id, self.graph.chromosome, self.graph.strand,
                                             path_id, path_out, edge_type_out, junc_ids_out, exon_mod])





    def training_junctions(self):
        with open(os.path.join(self.config.directory, 'training_junctions.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for junc_list in self.training_junction_cache.values():
                for common, edge, quants in junc_list:
                    exon1 = self.graph.start_node(edge)
                    exon2 = self.graph.end_node(edge)
                    junc_id = f"{self.gene_id}_{self.graph.chromosome}_{self.graph.strand}_{edge.absolute_start}_{edge.absolute_end}_{'I' if edge.ir else 'J'}"
                    writer.writerow(common + [junc_id, edge.absolute_start, edge.absolute_end,
                                     '', exon1.absolute_start, exon1.absolute_end,
                                     '', exon2.absolute_start, exon2.absolute_end,
                                     edge.ir] + quants)

    def cassette(self):
        with open(os.path.join(self.config.directory, 'cassette.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'cassette_exon':



                        src_common = self.common_data(module,
                                                      's',
                                                      node=event['C1'],
                                                      event_name="CE",
                                                      event_ii=event_i)
                        trg_common = self.common_data(module,
                                                      't',
                                                      node=event['C2'],
                                                      event_name="CE",
                                                      event_ii=event_i)

                        # strand case per-event trimming of exons
                        # we dont have access to the strand case / graph here.
                        # instead we check strand by checking if the first exon has coordinates before the last exon
                        # if not, we take the forward trimming definition, set all start -> end and end -> start
                        # the then for each range, reverse the order of the first and second coordinate
                        # (if this needs to be redesigned, then, only design with the + strand in mind, and
                        # follow these rules to make the reverse for the minus strand

                        c1_range_str = self._trim_strand_case_range_str(event['C1'], 'start', event['Include1'], 'start')
                        a_range_str = self._trim_strand_case_range_str(event['Include1'], 'end', event['Include2'], 'start')
                        c2_range_str = self._trim_strand_case_range_str(event['Include2'], 'end', event['C2'], 'end')
                        event_size = self._trim_strand_case_event_size(event['Include1'], 'end', event['Include2'], 'start') + 1

                        row = [event['Skip'].de_novo,
                               c1_range_str,
                               'C2',
                               c2_range_str,
                               'C1_C2',
                               event['Skip'].range_str()]

                        # quantifications consist of a number of fixed columns that don't depend on each other
                        # plus three special columns junction_changing, event_changing, event_non_changing
                        # ONE CALL TO QUANTIFICATIONS GETS RESULTS FOR ALL VOILA FILES IN ONE JUNCTION BY DEFAULT
                        # junction_changing: this works well for quantifications() as is, because it works
                        # off a single junction, and all voila files. We add this case to the main quantification
                        # writer class

                        # event_non_changing and event_changing: to know this information, we need to have had
                        # calculated non_changing
                        # for all junctions in the event, first. This should only be done once, so we break out
                        # this calculation into a separate class which can take in a list of junctions at once
                        quant_identifiers = (
                            ('s', event['Skip']),
                            ('s', event['Include1']),
                            ('t', event['Skip']),
                            ('t', event['Include2']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Skip'], event['C1'])
                        writer.writerow(src_common + row + [event_size] + quants + event['constitutive'])
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self._add_training_junction(src_common, event['Skip'], quants, module.idx)
                        self.heatmap_add(module, src_common, quants,
                                         event['Skip'].absolute_end - event['Skip'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Include1'].de_novo,
                               c1_range_str,
                               'A',
                               a_range_str,
                               'C1_A',
                               event['Include1'].range_str()]

                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Include1'], event['C1'])
                        writer.writerow(src_common + row + [event_size] + quants + event['constitutive'])
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self._add_training_junction(src_common, event['Include1'], quants, module.idx)
                        self.heatmap_add(module, src_common, quants,
                                         event['Include1'].absolute_end - event['Include1'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Skip'].de_novo,
                               c2_range_str,
                               'C1',
                               c1_range_str,
                               'C2_C1',
                               event['Skip'].range_str()]

                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Skip'], event['C2'])
                        writer.writerow(trg_common + row + [event_size] + quants + event['constitutive'])
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self._add_training_junction(trg_common, event['Skip'], quants, module.idx)
                        self.heatmap_add(module, trg_common, quants,
                                         event['Skip'].absolute_end - event['Skip'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Include2'].de_novo,
                               c2_range_str,
                               'A',
                               a_range_str,
                               'C2_A',
                               event['Include2'].range_str()]

                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Include2'], event['C2'])
                        writer.writerow(trg_common + row + [event_size] + quants + event['constitutive'])
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self._add_training_junction(trg_common, event['Include2'], quants, module.idx)
                        self.heatmap_add(module, trg_common, quants,
                                         event['Include2'].absolute_end - event['Include2'].absolute_start,
                                         row[0], row[4], row[5])
                        event_i += 1


    def alt3prime(self):
        with open(os.path.join(self.config.directory, 'alt3prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'alt3ss':
                        src_common = self.common_data(module,
                                                      's',
                                                      node=event['E1'],
                                                      event_name="A3",
                                                      event_ii=event_i)
                        trg_common = self.common_data(module,
                                                      't',
                                                      node=event['E2'],
                                                      event_name="A3",
                                                      event_ii=event_i)

                        event_size = self._trim_strand_case_event_size(event['Proximal'], 'end', event['Distal'], 'end')

                        # preferentially use source LSV
                        if src_common[5]:

                            quant_identifiers = (
                                ('s', event['Proximal']),
                                ('s', event['Distal']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            e1_range_str = self._trim_strand_case_range_str(event['E1'], 'start', event['Proximal'], 'start')
                            e2_range_str = self._trim_strand_case_range_str(event['Proximal'], 'end', event['E2'], 'end')

                            row = [event['Proximal'].de_novo,
                                   e1_range_str,
                                   'E2',
                                   e2_range_str,
                                   'Proximal',
                                   event['Proximal'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', edge=event['Proximal'], node=event['E1'])
                            writer.writerow(src_common + row + [event_size] + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['Proximal'].absolute_end - event['Proximal'].absolute_start,
                                             row[0], row[4], row[5])

                            row = [event['Distal'].de_novo,
                                   e1_range_str,
                                   'E2',
                                   e2_range_str,
                                   'Distal',
                                   event['Distal'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', edge=event['Distal'], node=event['E1'])
                            writer.writerow(src_common + row + [event_size] + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['Distal'].absolute_end - event['Distal'].absolute_start,
                                             row[0], row[4], row[5])
                        elif trg_common[5]:

                            quant_identifiers = (
                                ('t', event['Proximal']),
                                ('t', event['Distal']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            e1_range_str = f"{event['E1'].start}-{event['Proximal'].start}"
                            e2_range_str = f"{event['Proximal'].end}-{event['E2'].end}"

                            row = [event['Proximal'].de_novo,
                                   e2_range_str,
                                   'E1',
                                   e1_range_str,
                                   'Proximal',
                                   event['Proximal'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', edge=event['Proximal'], node=event['E2'])
                            writer.writerow(trg_common + row + [event_size] + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['Proximal'].absolute_end - event['Proximal'].absolute_start,
                                             row[0], row[4], row[5])

                            row = [event['Distal'].de_novo,
                                   e2_range_str,
                                   'E1',
                                   e1_range_str,
                                   'Distal',
                                   event['Distal'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', edge=event['Distal'], node=event['E2'])
                            writer.writerow(trg_common + row + [event_size] + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['Distal'].absolute_end - event['Distal'].absolute_start,
                                             row[0], row[4], row[5])
                        event_i += 1


    def alt5prime(self):
        with open(os.path.join(self.config.directory, 'alt5prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'alt5ss':

                        src_common = self.common_data(module,
                                                      's',
                                                      node=event['E1'],
                                                      event_name="A5",
                                                      event_ii=event_i)
                        trg_common = self.common_data(module,
                                                      't',
                                                      node=event['E2'],
                                                      event_name="A5",
                                                      event_ii=event_i)

                        event_size = self._trim_strand_case_event_size(event['Proximal'], 'start', event['Distal'], 'start')

                        # preferentially use target LSV
                        if trg_common[5]:

                            quant_identifiers = (
                                ('t', event['Proximal']),
                                ('t', event['Distal']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            e1_range_str = self._trim_strand_case_range_str(event['E1'], 'start', event['Proximal'], 'start')
                            e2_range_str = self._trim_strand_case_range_str(event['Proximal'], 'end', event['E2'], 'end')

                            row = [event['Proximal'].de_novo,
                                   e2_range_str,
                                   'E1',
                                   e1_range_str,
                                   'Proximal',
                                   event['Proximal'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', edge=event['Proximal'], node=event['E2'])
                            writer.writerow(trg_common + row + [event_size] + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['Proximal'].absolute_end - event['Proximal'].absolute_start,
                                             row[0], row[4], row[5])

                            row = [event['Distal'].de_novo,
                                   e2_range_str,
                                   'E1',
                                   e1_range_str,
                                   'Distal',
                                   event['Distal'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', edge=event['Distal'], node=event['E2'])
                            writer.writerow(trg_common + row + [event_size] + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['Distal'].absolute_end - event['Distal'].absolute_start,
                                             row[0], row[4], row[5])
                        elif src_common[5]:

                            quant_identifiers = (
                                ('s', event['Proximal']),
                                ('s', event['Distal']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            e1_range_str = self._trim_strand_case_range_str(event['E1'], 'start', event['Proximal'],
                                                                            'start')

                            e2_range_str = self._trim_strand_case_range_str(event['Proximal'], 'end', event['E2'],
                                                                            'end')

                            row = [event['Proximal'].de_novo,
                                   e1_range_str,
                                   'E2',
                                   e2_range_str,
                                   'Proximal',
                                   event['Proximal'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', edge=event['Proximal'], node=event['E1'])
                            writer.writerow(src_common + row + [event_size] + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['Proximal'].absolute_end - event['Proximal'].absolute_start,
                                             row[0], row[4], row[5])

                            row = [event['Distal'].de_novo,
                                   e1_range_str,
                                   'E2',
                                   e2_range_str,
                                   'Distal',
                                   event['Distal'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', edge=event['Distal'], node=event['E1'])
                            writer.writerow(src_common + row + [event_size] + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['Distal'].absolute_end - event['Distal'].absolute_start,
                                             row[0], row[4], row[5])
                        event_i += 1


    def p_alt5prime(self):
        with open(os.path.join(self.config.directory, 'p_alt5prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'p_alt5ss':
                        src_common = self.common_data(module,
                                                      's',
                                                      node=event['C1'],
                                                      event_name="pA5",
                                                      event_ii=event_i)
                        trg_common = self.common_data(module,
                                                      't',
                                                      node=event['C2'],
                                                      event_name="pA5",
                                                      event_ii=event_i)

                        event_size = self._trim_strand_case_event_size(event['Include1'], 'absolute_end',
                                                                       event['Include1'], 'absolute_start') + 1 + \
                                     self._trim_strand_case_event_size(event['A'], 'end', event['A'], 'start') + 1

                        quant_identifiers = (
                            ('s', event['Skip']),
                            ('s', event['Include1']),
                            ('t', event['Skip']),
                            ('t', event['Include2']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        c1_range_str = event['C1'].range_str()
                        a_range_str = self._trim_strand_case_range_str(event['A'], 'start', event['Include2'], 'start')
                        c2_range_str = self._trim_strand_case_range_str(event['Include2'], 'end', event['C2'], 'end')

                        row = [event['Skip'].de_novo,
                               c1_range_str,
                               'E3',
                               c2_range_str,
                               'Distal',
                               event['Skip'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Skip'])
                        writer.writerow(src_common + row + [event_size] + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, src_common, quants,
                                         event['Skip'].absolute_end - event['Skip'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Include1'].de_novo,
                               c1_range_str,
                               'E2',
                               a_range_str,
                               'Intron',
                               event['Include1'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Include1'])
                        writer.writerow(src_common + row + [event_size] + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, src_common, quants,
                                         event['Include1'].absolute_end - event['Include1'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Skip'].de_novo,
                               c2_range_str,
                               'E1',
                               c1_range_str,
                               'Distal',
                               event['Skip'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Skip'])
                        writer.writerow(trg_common + row + [event_size] + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, trg_common, quants,
                                         event['Skip'].absolute_end - event['Skip'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Include2'].de_novo,
                               c2_range_str,
                               'E2',
                               a_range_str,
                               'Proximal',
                               event['Include2'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Include2'])
                        writer.writerow(trg_common + row + [event_size] + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, trg_common, quants,
                                         event['Include2'].absolute_end - event['Include2'].absolute_start,
                                         row[0], row[4], row[5])

                        event_i += 1


    def p_alt3prime(self):
        with open(os.path.join(self.config.directory, 'p_alt3prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'p_alt3ss':
                        # TODO: why am I using strand_case for putative alt3, but not putative alt 5?
                        src_common = self.common_data(module,
                                                      's',
                                                      node=module.strand_case(case_plus=event['C1'],
                                                                              case_minus=event['C2']),
                                                      event_name="pA3",
                                                      event_ii=event_i)
                        trg_common = self.common_data(module,
                                                      't',
                                                      node=module.strand_case(case_plus=event['C2'],
                                                                              case_minus=event['C1']),
                                                      event_name="pA3",
                                                      event_ii=event_i)

                        event_size = self._trim_strand_case_event_size(event['Include2'], 'absolute_end',
                                                                       event['Include2'], 'absolute_start') + 1 + \
                                     self._trim_strand_case_event_size(event['A'], 'end', event['A'], 'start') + 1

                        quant_identifiers = (
                            ('s', event['Skip']),
                            ('s', event['Include1']),
                            ('t', event['Include2']),
                            ('t', event['Skip']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        c1_range_str = self._trim_strand_case_range_str(event['C1'], 'start', event['Include1'], 'start')
                        a_range_str = self._trim_strand_case_range_str(event['Include1'], 'end', event['A'], 'end')
                        c2_range_str = event['C2'].range_str()

                        row = [event['Skip'].de_novo,
                               c2_range_str,
                               'E3',
                               c1_range_str,
                               'Distal',
                               event['Skip'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Skip'])
                        writer.writerow(src_common + row + [event_size] + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, src_common, quants,
                                         event['Skip'].absolute_end - event['Skip'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Include1'].de_novo,
                               c2_range_str,
                               'E2',
                               a_range_str,
                               'Proximal',
                               event['Include1'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Include1'])
                        writer.writerow(src_common + row + [event_size] + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, src_common, quants,
                                         event['Include1'].absolute_end - event['Include1'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Include2'].de_novo,
                               c1_range_str,
                               'E2',
                               a_range_str,
                               'Intron',
                               event['Include2'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Include2'])
                        writer.writerow(trg_common + row + [event_size] + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, trg_common, quants,
                                         event['Include2'].absolute_end - event['Include2'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Skip'].de_novo,
                               c1_range_str,
                               'E1',
                               c2_range_str,
                               'Distal',
                               event['Skip'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Skip'])
                        writer.writerow(trg_common + row + [event_size] + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, trg_common, quants,
                                         event['Skip'].absolute_end - event['Skip'].absolute_start,
                                         row[0], row[4], row[5])

                        event_i += 1

    def alt3and5prime(self):
        with open(os.path.join(self.config.directory, 'alt3and5prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'alt3and5ss':
                        src_common = self.common_data(module, 's',
                                                      edge=event["J1"],
                                                      node=event["E1"],
                                                      event_name="A3A5",
                                                      event_ii=event_i)
                        trg_common = self.common_data(module, 't',
                                                      edge=event["J1"],
                                                      node=event["E2"],
                                                      event_name="A3A5",
                                                      event_ii=event_i)

                        quant_identifiers = (
                            ('s', event['J1']),
                            ('s', event['J2']),
                            ('t', event['J1']),
                            ('t', event['J2']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        e1_range_str = self._trim_strand_case_range_str(event['E1'], 'start', event['J2'], 'start')
                        e2_range_str = self._trim_strand_case_range_str(event['J1'], 'end', event['E2'], 'end')

                        row = [event['J1'].de_novo,
                               e1_range_str,
                               'E2',
                               e2_range_str,
                               'E1_E2_J1',
                               event['J1'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['J1'])
                        writer.writerow(src_common + row + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, src_common, quants,
                                         event['J1'].absolute_end - event['J1'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['J2'].de_novo,
                               e1_range_str,
                               'E2',
                               e2_range_str,
                               'E1_E2_J2',
                               event['J2'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['J2'])
                        writer.writerow(src_common + row + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, src_common, quants,
                                         event['J2'].absolute_end - event['J2'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['J1'].de_novo,
                               e2_range_str,
                               'E1',
                               e1_range_str,
                               'E2_E1_J1',
                               event['J1'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['J1'])
                        writer.writerow(trg_common + row + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, trg_common, quants,
                                         event['J1'].absolute_end - event['J1'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['J2'].de_novo,
                               e2_range_str,
                               'E1',
                               e1_range_str,
                               'E2_E1_J2',
                               event['J2'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['J2'])
                        writer.writerow(trg_common + row + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, trg_common, quants,
                                         event['J2'].absolute_end - event['J2'].absolute_start,
                                         row[0], row[4], row[5])
                        event_i += 1


    def mutually_exclusive(self):
        with open(os.path.join(self.config.directory, 'mutually_exclusive.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'mutually_exclusive':

                        quant_identifiers = (
                            ('s', event['Include1']),
                            ('s', event['SkipA1']),
                            ('t', event['Include2']),
                            ('t', event['SkipA2']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        c1_range_str = self._trim_strand_case_range_str(event['C1'], 'start', event['Include1'], 'start')
                        a1_range_str = self._trim_strand_case_range_str(event['Include1'], 'end', event['SkipA2'], 'start')
                        a2_range_str = self._trim_strand_case_range_str(event['SkipA1'], 'end', event['Include2'], 'start')
                        c2_range_str = self._trim_strand_case_range_str(event['Include2'], 'end', event['C2'], 'end')

                        row = [event['Include1'].de_novo,
                               c1_range_str,
                               'A1',
                               a1_range_str,
                               'C1_A1',
                               event['Include1'].range_str()]
                        src_common = self.common_data(module,
                                                      's',
                                                      node=event['C1'],
                                                      edge=event['Include1'],
                                                      event_name='mxe',
                                                      event_ii=event_i)
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + \
                                 self.quantifications(module, 's', edge=event['Include1'],node=event['C1'])
                        writer.writerow(src_common + row + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, src_common, quants,
                                         event['Include1'].absolute_end - event['Include1'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['SkipA1'].de_novo,
                               c1_range_str,
                               'A2',
                               a2_range_str,
                               'C1_A2',
                               event['SkipA1'].range_str()]
                        src_common = self.common_data(module,
                                                      's',
                                                      node=event['C1'],
                                                      edge=event['SkipA1'],
                                                      event_name='mxe',
                                                      event_ii=event_i)
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + \
                                 self.quantifications(module, 's', edge=event['SkipA1'],node=event['C1'])
                        writer.writerow(src_common + row + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, src_common, quants,
                                         event['SkipA1'].absolute_end - event['SkipA1'].absolute_start,
                                         row[0], row[4], row[5])

                        row = [event['Include2'].de_novo,
                               c2_range_str,
                               'A1',
                               a1_range_str,
                               'C2_A1',
                               event['Include2'].range_str()]
                        trg_common = self.common_data(module,
                                                      't',
                                                      node=event['C2'],
                                                      edge=event['Include2'],
                                                      event_name='mxe',
                                                      event_ii=event_i)
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + \
                                 self.quantifications(module, 't', edge=event['Include2'], node=event['C2'])
                        writer.writerow(trg_common + row + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, trg_common, quants,
                                        event['Include2'].absolute_end - event['Include2'].absolute_start,
                                        row[0], row[4], row[5])

                        row = [event['SkipA2'].de_novo,
                               c2_range_str,
                               'A2',
                               a2_range_str,
                               'C2_A2',
                               event['SkipA2'].range_str()]
                        trg_common = self.common_data(module,
                                                      't',
                                                      node=event['C2'],
                                                      edge=event['SkipA2'],
                                                      event_name='mxe',
                                                      event_ii=event_i)
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + \
                                 self.quantifications(module, 't', edge=event['SkipA2'], node=event['C2'])
                        writer.writerow(trg_common + row + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        self.heatmap_add(module, trg_common, quants,
                                         event['SkipA2'].absolute_end - event['SkipA2'].absolute_start,
                                         row[0], row[4], row[5])
                        event_i += 1


    def alternate_last_exon(self):
        with open(os.path.join(self.config.directory, 'alternate_last_exon.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'ale':
                        src_common = self.common_data(module,
                                                      's',
                                                      node=event['Reference'],
                                                      event_name='ale',
                                                      event_ii=event_i)
                        # only ever write ALEs that are quantified by LSV

                        if src_common[5]:

                            quant_identifiers = (
                                ('s', event['SkipA2']),
                                ('s', event['SkipA1']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            reference_range_str = self._trim_strand_case_range_str(event['Reference'],
                                                                                   'start', event['SkipA2'], 'start')
                            proximal_range_str = self._trim_strand_case_range_str(event['SkipA2'],
                                                                                  'end', event['Proximal'], 'end')

                            row = [event['SkipA2'].de_novo,
                                   reference_range_str,
                                   'A',
                                   proximal_range_str,
                                   'Proximal']
                            if event['SkipA2'].ir:
                                row.append('{}-{}'.format(event['SkipA2'].absolute_start, event['SkipA2'].absolute_end))
                            else:
                                row.append(event['SkipA2'].range_str())
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['SkipA2'], event['Reference'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['SkipA2'].absolute_end - event['SkipA2'].absolute_start,
                                             row[0], row[4], row[5])

                            reference_range_str = self._trim_strand_case_range_str(event['Reference'],
                                                                                   'start', event['SkipA1'], 'start')
                            distal_range_str = self._trim_strand_case_range_str(event['SkipA1'],
                                                                                'end', event['Distal'], 'end')

                            row = [event['SkipA1'].de_novo,
                                   reference_range_str,
                                   'A',
                                   distal_range_str,
                                   'Distal',
                                   event['SkipA1'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing +\
                                     self.quantifications(module, 's', event['SkipA1'], event['Reference'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['SkipA1'].absolute_end - event['SkipA1'].absolute_start,
                                             row[0], row[4], row[5])
                            event_i += 1


    def alternate_first_exon(self):
        with open(os.path.join(self.config.directory, 'alternate_first_exon.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'afe':
                        all_event_quants = []
                        trg_common = self.common_data(module,
                                                      't',
                                                      node=event['Reference'],
                                                      event_ii=event_i,
                                                      event_name='afe')

                        # only ever write AFEs that are quantified by LSV
                        if trg_common[5]:

                            quant_identifiers = (
                                ('t', event['SkipA1']),
                                ('t', event['SkipA2']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            reference_range_str = self._trim_strand_case_range_str(event['SkipA1'],
                                                                                   'end', event['Reference'], 'end')
                            proximal_range_str = self._trim_strand_case_range_str(event['Proximal'],
                                                                                  'start', event['SkipA1'], 'start')

                            row = [event['SkipA1'].de_novo,
                                   reference_range_str,
                                   'A',
                                   proximal_range_str,
                                   'Proximal']
                            if event['SkipA1'].ir:
                                row.append('{}-{}'.format(event['SkipA1'].absolute_start, event['SkipA1'].absolute_end))
                            else:
                                row.append(event['SkipA1'].range_str())
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['SkipA1'], event['Reference'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['SkipA1'].absolute_end - event['SkipA1'].absolute_start,
                                             row[0], row[4], row[5])

                            reference_range_str = self._trim_strand_case_range_str(event['SkipA2'],
                                                                                   'end', event['Reference'], 'end')
                            distal_range_str = self._trim_strand_case_range_str(event['Distal'],
                                                                                'start', event['SkipA2'], 'start')

                            row = [event['SkipA2'].de_novo,
                                   reference_range_str,
                                   'A',
                                   distal_range_str,
                                   'Distal',
                                   event['SkipA2'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + \
                                     self.quantifications(module, 't', event['SkipA2'], event['Reference'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['SkipA2'].absolute_end - event['SkipA2'].absolute_start,
                                             row[0], row[4], row[5])
                            event_i += 1


    def p_alternate_last_exon(self):
        with open(os.path.join(self.config.directory, 'p_alternate_last_exon.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'p_ale':
                        src_common = self.common_data(module,
                                                      's',
                                                      node=event['Reference'],
                                                      event_ii=event_i,
                                                      event_name="pALE")


                        if src_common[5]:

                            quant_identifiers = (
                                ('s', event['SkipA2']),
                                ('s', event['SkipA1']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            reference_range_str = self._trim_strand_case_range_str(event['Reference'],
                                                                                   'start', event['SkipA2'], 'start')

                            row = [event['SkipA2'].de_novo,
                                   reference_range_str,
                                   'A', event['Proximal'].range_str(),
                                   'Proximal',
                                   event['SkipA2'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['SkipA2'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['SkipA2'].absolute_end - event['SkipA2'].absolute_start,
                                             row[0], row[4], row[5])

                            reference_range_str = self._trim_strand_case_range_str(event['Reference'],
                                                                                   'start', event['SkipA1'], 'start')

                            row = [event['SkipA1'].de_novo,
                                   reference_range_str,
                                   'A',
                                   event['Distal'].range_str(),
                                   'Distal',
                                   event['SkipA1'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['SkipA1'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['SkipA1'].absolute_end - event['SkipA1'].absolute_start,
                                             row[0], row[4], row[5])
                            event_i += 1


    def p_alternate_first_exon(self):
        with open(os.path.join(self.config.directory, 'p_alternate_first_exon.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'p_afe':

                        trg_common = self.common_data(module,
                                                      't',
                                                      node=event['Reference'],
                                                      event_ii=event_i,
                                                      event_name="pAFE")

                        if trg_common[5]:

                            quant_identifiers = (
                                ('t', event['SkipA1']),
                                ('t', event['SkipA2']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            reference_range_str = self._trim_strand_case_range_str(event['SkipA1'], 'end',
                                                                                   event['Reference'], 'end')

                            row = [event['SkipA1'].de_novo,
                                   reference_range_str,
                                   'A', event['Proximal'].range_str(),
                                   'Proximal',
                                   event['SkipA1'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['SkipA1'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['SkipA1'].absolute_end - event['SkipA1'].absolute_start,
                                             row[0], row[4], row[5])

                            reference_range_str = self._trim_strand_case_range_str(event['SkipA2'],
                                                                                   'end', event['Reference'], 'end')

                            row = [event['SkipA2'].de_novo,
                                   reference_range_str,
                                   'A',
                                   event['Distal'].range_str(),
                                   'Distal',
                                   event['SkipA2'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['SkipA2'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['SkipA2'].absolute_end - event['SkipA2'].absolute_start,
                                             row[0], row[4], row[5])
                            event_i += 1


    def alternative_intron(self):
        with open(os.path.join(self.config.directory, 'alternative_intron.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'alternative_intron':

                        # put coordinates back to Jordi's offset numbers

                        c1_node = module.strand_case(case_plus=event['C1'],
                                                     case_minus=event['C2'])
                        c2_node = module.strand_case(case_plus=event['C2'],
                                                     case_minus=event['C1'])
                        src_common = self.common_data(module,
                                                      's',
                                                      event['Intron'],
                                                      node=c1_node,
                                                      event_ii=event_i,
                                                      event_name="AI")
                        trg_common = self.common_data(module,
                                                      't',
                                                      event['Intron'],
                                                      node=c2_node,
                                                      event_ii=event_i,
                                                      event_name="AI")

                        event_size = abs(event['Intron'].absolute_end - event['Intron'].absolute_start) + 1

                        # I think this fails when there is a source LSV in the middle exon
                        # if any(':t:' in _l for _l in event['Intron'].lsvs) and not \
                        #    any(':s:' in _l for _l in event['Intron'].lsvs):
                        # The quantifications for the alternative intron must be from EITHER Target or Source point of view
                        #   To check if we should use the target point of view:
                        junc_in_trg_lsv_count = 0

                        # ensure the junction and intron are both quantified by the same target LSV
                        if len(set(trg_common[5].split(";")) & set(event['Spliced'].lsvs.keys())) == 1:
                            # if yes, add one to the junc count
                            junc_in_trg_lsv_count += 1
                        # if every junction passed above test, we know the Target point of view quantified the intron
                        # as well as all possible 'spliced' junctions
                        if junc_in_trg_lsv_count:

                            quant_identifiers = (
                                ('t', event['Intron']),
                                ('t', event['Spliced']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            row = [event['Intron'].de_novo,
                                   c2_node.range_str(),
                                   'C1',
                                   c1_node.range_str(),
                                   'C2_C1_intron',
                                   event['Intron'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module,
                                                          't',
                                                          edge=event['Intron'],
                                                          node=c2_node)
                            writer.writerow(trg_common + row + [event_size] + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['Intron'].absolute_end - event['Intron'].absolute_start,
                                             row[0], row[4], row[5])

                            row = [event['Spliced'].de_novo,
                                   c1_node.range_str(),
                                   'C1',
                                   c1_node.range_str(),
                                   'C2_C1_spliced',
                                   event['Spliced'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module,
                                                          't',
                                                          edge=event['Spliced'],
                                                          node=c2_node)
                            writer.writerow(trg_common + row + [event_size] + quants)
                            self.junction_cache.append((module, trg_common, quants,
                                                        event['Spliced'].de_novo, row[4], event['Spliced'].range_str()))
                            self.heatmap_add(module, trg_common, quants,
                                             event['Spliced'].absolute_end - event['Spliced'].absolute_start,
                                             event['Spliced'].de_novo, row[4], event['Spliced'].range_str())

                        # Else the intron and 'spliced' junctions are quantified from source point of view...
                        else:

                            quant_identifiers = (
                                ('s', event['Intron']),
                                ('s', event['Spliced']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            row = [event['Intron'].de_novo,
                                   c1_node.range_str(),
                                   'C2',
                                   c2_node.range_str(),
                                   'C1_C2_intron',
                                   event['Intron'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', edge=event['Intron'], node=c1_node)
                            writer.writerow(src_common + row + [event_size] + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            self.heatmap_add(module, src_common, quants,
                                             event['Intron'].absolute_end - event['Intron'].absolute_start,
                                             row[0], row[4], row[5])

                            row = [event['Spliced'].de_novo,
                                   c1_node.range_str(),
                                   'C2',
                                   c2_node.range_str(),
                                   'C1_C2_spliced',
                                   event['Spliced'].range_str()]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module,
                                                          's',
                                                          edge=event['Spliced'],
                                                          node=c1_node)
                            writer.writerow(src_common + row + [event_size] + quants)
                            self.junction_cache.append((module, src_common, quants,
                                                        event['Spliced'].de_novo, row[4], event['Spliced'].range_str()))
                            self.heatmap_add(module, src_common, quants,
                                             event['Spliced'].absolute_end - event['Spliced'].absolute_start,
                                             event['Spliced'].de_novo, row[4], event['Spliced'].range_str())

                        event_i += 1


    def multi_exon_spanning(self):
        with open(os.path.join(self.config.directory, 'multi_exon_spanning.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'multi_exon_spanning':

                        quant_identifiers = (
                            ('s', event['Skip']),
                            ('s', event['Include1']),
                            ('t', event['Skip']),
                            ('t', event['Include2'])
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        c1_range_str = self._trim_strand_case_range_str(event['C1'], 'start', event['Skip'], 'start')
                        c2_range_str = self._trim_strand_case_range_str(event['Skip'], 'end', event['C2'], 'end')
                        # C1_A
                        # how do you trim?
                        # did I use strand case correctly here?
                        if event['Include1']:
                            c1_a_range_str = module.strand_case(event['As'][0], event['As'][-1]).range_str()
                        else:
                            c1_a_range_str = ""
                        # C2_A_Last
                        # how do you trim?
                        # did I use strand case correctly here?
                        if event['Include2']:
                            c2_alast_range_str = module.strand_case(event['As'][-1], event['As'][0]).range_str()
                        else:
                            c2_alast_range_str = ""

                        # Source LSV side
                        row = ['Distal',
                               event['Skip'].range_str(),  # junction coord
                               event['Skip'].de_novo, # de novo?
                               event['Skip'].ir,
                               c1_range_str, # reference exon
                               c2_range_str, # exon spliced with
                               self.semicolon((x.range_str() for x in event['As'])), # exons spanned
                               len(event['As'])] # num exons spanned
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + \
                                 self.quantifications(module, 's', event['Skip'], event['C1'])

                        common = self.common_data(module,
                                                  's',
                                                  node=event['C1'],
                                                  edge=event['Skip'],
                                                  event_ii=event_i,
                                                  event_name="MES")

                        writer.writerow(common + row + quants)
                        self.junction_cache.append((module, common, quants, row[2], 'C1_C2', row[1]))

                        if event['Include1']:
                            a_1 = self._trim_strand_case_range_str(event['Skip'], 'end', event['C2'], 'end')
                        else:
                            a_1 = ""
                        row = ['Proximal',
                               event['Include1'].range_str() if event['Include1'] else "",  # junction coord
                               event['Include1'].de_novo if event['Include1'] else "",  # de novo?
                               event['Include1'].ir if event['Include1'] else "",
                               c1_range_str,  # reference exon
                               c1_a_range_str,  # exon spliced with
                               self.semicolon((x.range_str() for x in event['As'])),  # exons spanned
                               len(event['As'])]  # num exons spanned
                        quants = [
                                     event_non_changing, event_changing
                                 ] + self.quantifications(module, 's', event['Include1'], event['C1'])
                        writer.writerow(common + row + quants)
                        self.junction_cache.append((module, common, quants, row[2], 'C1_A', row[1]))

                        # Target LSV side
                        row = ['Distal',
                               event['Skip'].range_str(),  # junction coord
                               event['Skip'].de_novo,  # de novo?
                               event['Skip'].ir,
                               c2_range_str,  # reference exon
                               c1_range_str,  # exon spliced with
                               self.semicolon((x.range_str() for x in event['As'])),  # exons spanned
                               len(event['As'])]  # num exons spanned
                        quants = [
                                     event_non_changing, event_changing
                                 ] + self.quantifications(module, 't', event['Skip'], event['C2'])
                        common = self.common_data(module,
                                                  't',
                                                  node=event['C2'],
                                                  edge=event['Skip'],
                                                  event_ii=event_i,
                                                  event_name="MES")
                        writer.writerow(common + row + quants)
                        self.junction_cache.append((module, common, quants, row[2], 'C2_C1', row[1]))

                        row = ['Proximal',
                               event['Include2'].range_str() if event['Include2'] else "",  # junction coord
                               event['Include2'].de_novo  if event['Include2'] else "",  # de novo?
                               event['Include2'].ir if event['Include2'] else "",
                               c1_range_str,  # reference exon
                               c2_alast_range_str,  # exon spliced with
                               self.semicolon((x.range_str() for x in event['As'])),  # exons spanned
                               len(event['As'])]  # num exons spanned
                        quants = [
                                     event_non_changing, event_changing
                                 ] + self.quantifications(module, 't', event['Include2'], event['C2'])

                        writer.writerow(common + row + quants)
                        self.junction_cache.append((module, common, quants, row[2], 'A_C2', row[1]))


                        event_i += 1


    def tandem_cassette(self):
        with open(os.path.join(self.config.directory, 'tandem_cassette.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'tandem_cassette':

                        quant_identifiers = (
                            ('s', event['Skip']),
                            ('s', event['Include1']),
                            ('t', event['Skip']),
                            ('t', event['Include2']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        c1_range_str = self._trim_strand_case_range_str(event['C1'], 'start', event['Skip'], 'start')
                        c2_range_str = self._trim_strand_case_range_str(event['Skip'], 'end', event['C2'], 'end')

                        event_sizes = []
                        tandem_exons_range_strs = []
                        tandem_exons_range_strs.append(self._trim_strand_case_range_str(event['Include1'], 'end', event['Includes'][0], 'start'))
                        event_sizes.append(self._trim_strand_case_event_size(event['Include1'], 'end', event['Includes'][0], 'start') + 1)
                        for include_junc_i in range(len(event['Includes'])-1):
                            tandem_exons_range_strs.append(self._trim_strand_case_range_str(event['Includes'][include_junc_i], 'end', event['Includes'][include_junc_i+1], 'start'))
                            event_sizes.append(self._trim_strand_case_event_size(event['Includes'][include_junc_i], 'end', event['Includes'][include_junc_i+1], 'start') + 1)
                        tandem_exons_range_strs.append(self._trim_strand_case_range_str(event['Includes'][-1], 'end', event['Include2'], 'start'))
                        event_sizes.append(self._trim_strand_case_event_size(event['Includes'][-1], 'end', event['Include2'], 'start') + 1)
                        tandem_exons_range_strs = self.semicolon(tandem_exons_range_strs)
                        event_size = sum(event_sizes)

                        row = [event['Skip'].de_novo,
                               c1_range_str,
                               'C2',
                               c2_range_str,
                               tandem_exons_range_strs,
                               len(event['Tandem_Exons']),
                               'C1_C2',
                               event['Skip'].range_str()]
                        common = self.common_data(module,
                                                  's',
                                                  node=event['C1'],
                                                  edge=event['Skip'],
                                                  event_ii=event_i,
                                                  event_name="TCE")
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Skip'], event['C1'])
                        writer.writerow(common + row + [event_size] + quants)
                        self.junction_cache.append((module, common, quants, row[0], row[6], row[7]))
                        self.heatmap_add(module, common, quants,
                                         event['Skip'].absolute_end - event['Skip'].absolute_start,
                                         row[0], row[6], row[7])

                        row = [event['Include1'].de_novo,
                               c1_range_str,
                               'A1',
                               event['Tandem_Exons'][0].range_str(),
                               tandem_exons_range_strs,
                               len(event['Tandem_Exons']),
                               'C1_A',
                               event['Include1'].range_str()]
                        common = self.common_data(module,
                                                  's',
                                                  node=event['C1'],
                                                  edge=event['Include1'],
                                                  event_ii=event_i,
                                                  event_name="TCE")
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Include1'], event['C1'])
                        writer.writerow(common + row + [event_size] + quants)
                        self.junction_cache.append((module, common, quants, row[0], row[6], row[7]))
                        self.heatmap_add(module, common, quants,
                                         event['Include1'].absolute_end - event['Include1'].absolute_start,
                                         row[0], row[6], row[7])

                        row = [event['Skip'].de_novo,
                               event['C2'].range_str(),
                               'C1',
                               c1_range_str,
                               tandem_exons_range_strs,
                               len(event['Tandem_Exons']),
                               'C2_C1',
                               event['Skip'].range_str()]
                        common = self.common_data(module,
                                                  't',
                                                  node=event['C2'],
                                                  edge=event['Skip'],
                                                  event_ii=event_i,
                                                  event_name="TCE")
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Skip'], event['C2'])
                        writer.writerow(common + row + [event_size] + quants)
                        self.junction_cache.append((module, common, quants, row[0], row[6], row[7]))
                        self.heatmap_add(module, common, quants,
                                         event['Skip'].absolute_end - event['Skip'].absolute_start,
                                         row[0], row[6], row[7])

                        row = [event['Include2'].de_novo,
                               c2_range_str,
                               'A_Last',
                               event['Tandem_Exons'][-1].range_str(),
                               tandem_exons_range_strs,
                               len(event['Tandem_Exons']),
                               'C2_A_Last',
                               event['Include2'].range_str()]
                        common = self.common_data(module,
                                                  't',
                                                  node=event['C2'],
                                                  edge=event['Include2'],
                                                  event_ii=event_i,
                                                  event_name="TCE")
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Include2'], event['C2'])
                        writer.writerow(common + row + [event_size] + quants)
                        self.junction_cache.append((module, common, quants, row[0], row[6], row[7]))
                        self.heatmap_add(module, common, quants,
                                         event['Include2'].absolute_end - event['Include2'].absolute_start,
                                         row[0], row[6], row[7])
                        event_i += 1


    def exitron(self):
        with open(os.path.join(self.config.directory, 'exitron.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'exitron':

                        quant_identifiers = (
                            ('s', event['Junc']),
                            ('t', event['Junc']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        src_common = self.common_data(module,
                                                      's',
                                                      event_ii=event_i,
                                                      event_name="EX")
                        trg_common = self.common_data(module,
                                                      't',
                                                      event_ii=event_i,
                                                      event_name="EX")
                        row = [event['Junc'].de_novo,
                               event['Exon'].range_str(),
                               event['Junc'].range_str()]
                        # just in case exitrons are ever quantified, somehow, *try* to get
                        # the exitron's junction quantification (won't exist for now... which is desired)
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', edge=event['Junc'])
                        writer.writerow(src_common + row + quants)
                        self.junction_cache.append((module, src_common, quants, row[0], 'exitron', row[2]))
                        row = [event['Junc'].de_novo,
                               event['Exon'].range_str(),
                               event['Junc'].range_str()]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', edge=event['Junc'])
                        writer.writerow(trg_common + row + quants)
                        self.junction_cache.append((module, trg_common, quants, row[0], 'exitron', row[2]))
                        event_i += 1

    def orphan_junction(self):
        with open(os.path.join(self.config.directory, 'orphan_junction.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'orphan_junction':

                        quant_identifiers = (
                            ('s', event['Junc']),
                            ('t', event['Junc']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        src_common = self.common_data(module,
                                                      's',
                                                      event_ii=event_i,
                                                      event_name="OR")
                        trg_common = self.common_data(module,
                                                      't',
                                                      event_ii=event_i,
                                                      event_name="OR")

                        row = [event['Junc'].de_novo,
                              event['A1'].range_str(),
                              event['A2'].range_str(),
                              event['Junc'].range_str()]

                        if src_common[5]:

                            # just in case exitrons are ever quantified, somehow, *try* to get
                            # the exitron's junction quantification (won't exist for now... which is desired)
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', edge=event['Junc'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], 'Orphan', row[3]))
                            self.heatmap_add(module, src_common, quants,
                                             event['Junc'].absolute_end - event['Junc'].absolute_start,
                                             row[0], 'Orphan', row[3])

                            event_i += 1

                        elif trg_common[5]:

                            # just in case exitrons are ever quantified, somehow, *try* to get
                            # the exitron's junction quantification (won't exist for now... which is desired)
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't',
                                                                                                 edge=event['Junc'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], 'Orphan', row[3]))
                            self.heatmap_add(module, trg_common, quants,
                                             event['Junc'].absolute_end - event['Junc'].absolute_start,
                                             row[0], 'Orphan', row[3])

                            event_i += 1

    def constitutive(self):
        with open(os.path.join(self.config.directory, 'constitutive.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'constitutive':

                        common = self.common_data(module,
                                                  event_ii=event_i,
                                                  event_name="CJ")

                        quant_identifiers = (
                            ('s', event['Junc']),
                            ('t', event['Junc']),
                        )
                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        c1_range_str = self._trim_strand_case_range_str(event['C1'], 'start', event['Junc'], 'start')
                        c2_range_str = self._trim_strand_case_range_str(event['Junc'], 'end', event['C2'], 'end')

                        row = [event['Junc'].de_novo,
                               c2_range_str,
                               'C1', c1_range_str,
                               'C2_C1',
                               event['Junc'].range_str(),
                               'False',
                               module.collapsed_event_name]
                        quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, edge=event['Junc'])
                        writer.writerow(common + row + quants)
                        self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))

                    elif event['event'] == 'constitutive_intron':
                        src_common = self.common_data(module,
                                                      's',
                                                      event_ii=event_i,
                                                      event_name="CI")
                        trg_common = self.common_data(module,
                                                      't',
                                                      event_ii=event_i,
                                                      event_name="CI")

                        if any(':t:' in _l for _l in event['Intron'].lsvs) and not \
                           any(':s:' in _l for _l in event['Intron'].lsvs):

                            quant_identifiers = (
                                ('t', event['Intron']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            row = [event['Intron'].de_novo,
                                   event['C2'].range_str(),
                                   'C1',
                                   event['C1'].range_str(),
                                   'C2_C1_intron',
                                   event['Intron'].range_str(),
                                   'True',
                                   module.collapsed_event_name]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 't', event['Intron'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                        else:
                            quant_identifiers = (
                                ('s', event['Intron']),
                            )
                            event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                            event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                            row = [event['Intron'].de_novo,
                                   event['C1'].range_str(),
                                   'C2',
                                   event['C2'].range_str(),
                                   'C1_C2_intron',
                                   event['Intron'].range_str(),
                                   'True',
                                   module.collapsed_event_name]
                            quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + self.quantifications(module, 's', event['Intron'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                    event_i += 1

    def other_event(self):
        with open(os.path.join(self.config.directory, 'other.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if _complex and self.config.only_binary:
                    continue
                event_i = 1
                for event in events:
                    if event['event'] == 'other_event':
                        all_event_quants = []
                        quant_identifiers = []

                        juncs = event['edges']
                        exons = event['nodes']
                        lsvs = event['lsvs']

                        # gets all LSVs
                        row_common = self.common_data(module,
                                                      parity=None,
                                                      edge=None,
                                                      node=None,
                                                      event_ii=None,
                                                      event_name='other')
                        # update lsvs with relevant LSVs
                        row_common[5]  = ";".join(lsvs)
                        juncs_l = [j.range_str() for j in juncs]
                        juncs_str = ";".join(juncs_l)
                        exons_str = ";".join([e.range_str() for e in exons])
                        row = [juncs_str, exons_str]

                        writer.writerow(row_common + row)
                        # still want 1 row per quantification in junctions.tsv and heatmap.tsv
                        junc_i = 0
                        # for n1, n2 in combinations(exons, 2):
                        lsvs_seen = []
                        seen_junc = {
                            "s": None,
                            "t": None
                        }

                        for junc in juncs:
                            junc_i += 1
                            this_node = junc.node

                            for lsv_type in ["s", "t"]:
                                all_event_quants.append(self.quantifications(module,
                                                              lsv_type,
                                                              edge=junc,
                                                              node=this_node))

                                quant_identifiers.append((lsv_type, junc,))

                        event_non_changing, cases_non_changing = self.event_non_changing(module, quant_identifiers)
                        event_changing, cases_changing = self.event_changing(module, quant_identifiers)

                        for junc in juncs:
                            junc_i += 1
                            this_node = junc.node
                            # TODO: handle na type?
                            for lsv_type in ["s", "t"]:
                                j_quants = [event_non_changing, event_changing] + cases_non_changing + cases_changing + all_event_quants.pop(0)
                                j_common = self.common_data(module,
                                              parity=lsv_type,
                                              edge=junc,
                                              node=this_node,
                                              event_ii=None,
                                              event_name='other')
                                # de_novo, junction_name, coordinates
                                both_infos = [
                                    junc.de_novo,
                                    "J%s" % junc_i,  # junction name
                                    junc.range_str()
                                ]
                                # module, common_data, quantifications, de_novo, junction_name, coordinates
                                junction_cache_row = (module,
                                     j_common,
                                     j_quants,
                                     both_infos[0],
                                     both_infos[1],
                                     both_infos[2]
                                     )
                                # module, common, quants, junc_len, denovo, juncname, junccoord
                                heatmap_row = (module,
                                    j_common,
                                    j_quants,
                                    junc.absolute_end - junc.absolute_start,
                                    both_infos[0],
                                    both_infos[1],
                                    both_infos[2]
                                )
                                if seen_junc[lsv_type] is None:
                                    seen_junc[lsv_type] = (junction_cache_row, heatmap_row)
                                # if LSV here
                                if len(j_common[5]) > 0:
                                    lsvs_seen.append(lsv_type)
                                    seen_junc[lsv_type] = (junction_cache_row, heatmap_row)
                            # only add 1 row per junction
                            # (or 2, if the junction was quantifable by source & target LSV...
                            if len(lsvs_seen) == 0:
                                lsvs_seen = ["s"]
                            for lsv_type in lsvs_seen:
                                self.junction_cache.append(
                                    seen_junc[lsv_type][0]
                                )
                                self.heatmap_add(
                                    seen_junc[lsv_type][1][0],
                                    seen_junc[lsv_type][1][1],
                                    seen_junc[lsv_type][1][2],
                                    seen_junc[lsv_type][1][3],
                                    seen_junc[lsv_type][1][4],
                                    seen_junc[lsv_type][1][5],
                                    seen_junc[lsv_type][1][6]
                                )
                    event_i += 1


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

    def heatmap_add(self, module, common, quants, junc_len, denovo, juncname, junccoord):
        """
        Conditionally add a row to heatmap cache by comparing it to what exists there already
            **The final heatmap should have the shortest quantifiable junction from the module.**
            If none of the junctions in the module were quantifiable, heatmap has the shortest junction
        """

        if self.config.heatmap_selection == 'max_abs_dpsi':

            try:
                # i + 2, because quants contains event_changing and event_non_changing which are
                # not considered for offsets by the rest of the quantifications class.
                max_abs_dpsi = max((abs(float(quants[i+2])) if quants[i+2] else 0.0 for i in self.dpsi_quant_idxs))

            except ValueError:
                # this problem only happens with semicolon separated values, which we are planning to get rid of
                # which is why it is just try/except for now
                max_abs_dpsi = 0.0

            if not module.idx in self.heatmap_cache:
                self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, max_abs_dpsi)
            else:
                # if existing junc lacks an LSV (thus no quantification)
                if not self.heatmap_cache[module.idx][1][5]:
                    # if new junc does have LSV, update
                    if common[5]:
                        self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, max_abs_dpsi)
                    # if there is a tie, retain the shorter junction
                    elif round(self.heatmap_cache[module.idx][7], 3) == round(max_abs_dpsi, 3):
                        if junc_len < self.heatmap_cache[module.idx][3]:
                            self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, max_abs_dpsi)
                    # else, if new max dpsi larger, update
                    elif self.heatmap_cache[module.idx][7] < max_abs_dpsi:
                        self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, max_abs_dpsi)
                # else existing junc has LSV
                else:
                    # only if new junc has LSV
                    if common[5]:
                        # and if new max dpsi ties, and junction is shorter
                        if round(self.heatmap_cache[module.idx][7], 3) == round(max_abs_dpsi, 3):
                            if junc_len < self.heatmap_cache[module.idx][3]:
                                self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, max_abs_dpsi)
                        # or if new max dpsi is larger
                        elif self.heatmap_cache[module.idx][7] < max_abs_dpsi:
                            self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, max_abs_dpsi)

        else:
            if not module.idx in self.heatmap_cache:
                self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, None)
            else:
                # if existing junc lacks an LSV (thus no quantification)
                if not self.heatmap_cache[module.idx][1][5]:
                    # if new junc does have LSV, update
                    if common[5]:
                        self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, None)
                    # else, if new junc shorter, update
                    elif self.heatmap_cache[module.idx][3] > junc_len:
                        self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, None)
                # else existing junc has LSV
                else:
                    # only if new junc has LSV
                    if common[5]:
                        # and if new junc is shorter
                        if self.heatmap_cache[module.idx][3] > junc_len:
                            self.heatmap_cache[module.idx] = (module, common, quants, junc_len, denovo, juncname, junccoord, None)

    def junctions(self):
        """
        Write a file with a listing of all junctions
        :return:
        """
        with open(os.path.join(self.config.directory, 'junctions.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module, common_data, quantifications, de_novo, junction_name, coordinates in self.junction_cache:
                writer.writerow(common_data + [module.collapsed_event_name, de_novo, junction_name, coordinates] + quantifications)

    def mpe(self):
        """
        Write a file with a listing of all primerable regions
        :return:
        """
        with open(os.path.join(self.config.directory, 'mpe_primerable_regions.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events = self.mpe_regions[module.idx]
                event_i = 1
                for event in events:
                    n_assigned_lids = 0
                    lsvid = ""
                    edge_type = ""
                    if event['event'] == "mpe_source":

                        common = self.common_data(module,
                                                  's',
                                                  node=event['reference_exon'],
                                                  event_ii=event_i,
                                                  event_name="MPE")

                        eventtype = "Single Source" # SingleSource
                        constitutive_direction = "Upstream"
                        if event['edge_of_transcript']:
                            edge_type = "first_exon"
                    else:
                        common = self.common_data(module,
                                                  't',
                                                  node=event['reference_exon'],
                                                  event_ii=event_i,
                                                  event_name="MPE")

                        eventtype = "Single Target" # SingleTarget
                        constitutive_direction = "Downstream"
                        if event['edge_of_transcript']:
                            edge_type = "last_exon"
                    lsvids = common[-1].split(";")
                    if len(lsvids) > 1:
                        raise ValueError("Multiple LSV IDs (%s) unexpected ... %s" % (lsvids, event))
                    if event['at_module_edge']:
                        isfirst = "True"
                    else:
                        isfirst = "False"
                    constitutive_regions = event['constitutive_regions']
                    constitutive_coords = ";".join([region.range_str() for region in constitutive_regions])
                    constitutive_denovo = ";".join([str(region.is_de_novo()) for region in constitutive_regions])
                    constitutive_types = ";".join([region.short_name for region in constitutive_regions])
                    ref_exon = event['reference_exon']
                    ref_exon_coord = ref_exon.range_str()
                    ref_exon_exitrons = ";".join(ref_exon.get_exitrons())
                    const_reg = ref_exon.get_constant_region()
                    if const_reg == ref_exon_coord:
                        was_trimmed = "False"
                    else:
                        was_trimmed = "True"
                    if "events" in self.config.enabled_outputs or "summary" in self.config.enabled_outputs:
                        collapsed_event_name = module.collapsed_event_name
                    else:
                        collapsed_event_name = "ND"
                    row = common # ModID, GeneID, GeneName, Chr, Strand, LSV(s)
                    row += [collapsed_event_name, eventtype, isfirst, edge_type]
                    row += [ref_exon_coord, ref_exon.is_de_novo(), ref_exon_exitrons]
                    row += [const_reg, was_trimmed, constitutive_direction]
                    row += [constitutive_coords, constitutive_denovo, constitutive_types]
                    writer.writerow(row)
                    event_i += 1

    def heatmap(self):
        """
        Write the easily excel-able file
        This is similar to the other non-summary event files and also part of the summary file
        It lists rows like the non-summary files, but only chooses the shortest junction
        """
        with open(os.path.join(self.config.directory, 'heatmap.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module, common_data, quantifications, _, de_novo, junction_name, coordinates, _ in self.heatmap_cache.values():
                writer.writerow(common_data + [module.collapsed_event_name, de_novo, junction_name, coordinates] + quantifications)

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
                if _complex and self.config.only_binary:
                    continue
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
                counts['other_event'] = 0
                if self.config.keep_constitutive:
                    counts['constitutive'] = 0
                    counts['constitutive_intron'] = 0
                counts['multi_exon_spanning'] = 0
                other_novel_count = 0
                for event in events:
                    if event['event'] in counts:
                        counts[event['event']] += 1
                        if event['event'] == 'other_event' and event.get('novel', False) is True:
                            other_novel_count += 1
                de_novo_junctions = 0
                de_novo_introns = 0
                for edge in module.get_all_edges(ir=True):
                    if edge.de_novo:
                        if edge.ir:
                            de_novo_introns += 1
                        else:
                            de_novo_junctions += 1

                # we store collapsed event name on module, because we need it for constitutive
                module.collapsed_event_name = self._collapsed_event_name(counts, other_novel_count)

                writer.writerow(
                    ["%s_%d" % (self.gene_id, module.idx),
                     self.gene_id,
                     self.graph.gene_name,
                     self.graph.chromosome,
                     self.graph.strand,
                     self.semicolon(module.target_lsv_ids.union(module.source_lsv_ids))] +
                                [v if v else '' for v in counts.values()] +
                                [str(_complex),
                                 str(de_novo_junctions),
                                 str(de_novo_introns),
                                 str(_total_events),
                                 module.collapsed_event_name
                                 ]
                )



