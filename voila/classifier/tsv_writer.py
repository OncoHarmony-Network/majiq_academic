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
from voila.classifier.quantification_finder import QuantificationWriter
import numpy as np
from voila.vlsv import get_expected_psi, matrix_area
from itertools import combinations
from operator import itemgetter


summaryVars2Headers = {
    'cassette_exon': 'Cassette',
    'tandem_cassette': 'Tandem Cassette',
    'alt3ss': 'Alt 3',
    'alt5ss': 'Alt 5',
    'p_alt3ss': 'P_Alt 3',
    'p_alt5ss': 'P_Alt 5',
    'alt3and5ss': 'Alt 3 and Alt 5',
    'mutually_exclusive': 'MXE',
    'alternative_intron': 'Alternative Intron',
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

class BaseTsvWriter(QuantificationWriter):
    """
    Output AS data from one gene
    """
    def __init__(self, graph, gene_id):
        """
        :param output_path: The folder where all output TSV files will be written under
        :param graph: the Graph object of the gene
        """
        super().__init__()

        self.common_headers = ['Module ID', 'Gene ID', 'Gene Name', 'Chr', 'Strand']

        self.graph = graph
        self.gene_id = gene_id


        self.pid = multiprocessing.current_process().pid

        self.heatmap_cache = {}
        self.junction_cache = []



    @property
    def quantification_headers(self):
        return list(self.quantifications_int.keys())

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

    def common_data(self, module, parity=None, edge=None, node=None):
        """
        Extract the certain cols from the TSV which are generally similar across all outputs,

        """
        lsvs = self.parity2lsv(module, parity, edge, node)

        out = ["%s_%d" % (self.gene_id, module.idx), self.gene_id, self.graph.gene_name,
               self.graph.chromosome, self.graph.strand]

        if self.config.output_complex:
            out.append(str(module.is_complex))
        out.append(self.semicolon(lsvs))
        return out


    def start_headers(self, headers, filename):
        """
        Start a tsv file with the required headers, only if it does not yet exist

        """
        if not os.path.exists(os.path.join(self.config.directory, filename)):
            with open(os.path.join(self.config.directory, filename), 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
                writer.writerow(headers)

    def _collapsed_event_name(self, counts):
        """
        function to generate the collapsed event name from event counts
        """
        out = []
        for count in counts:
            if counts[count]:
                out.append("%sx%d" % (summaryVars2Headers[count], counts[count]))
        return '_'.join(out)


class TsvWriter(BaseTsvWriter):

    def __init__(self, graph, gene_id):
        super().__init__(graph, gene_id)

        if self.config.output_complex:
            self.common_headers.append('Complex')
        self.common_headers.append('LSV ID(s)')

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
                 'tandem_cassette.tsv', 'exitron.tsv', 'p_alternate_last_exon.tsv', 'p_alternate_first_exon.tsv']
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

        return names


    def start_all_headers(self):

        if self.config.putative_multi_gene_regions:
            headers = self.common_headers + ['Exon1Start', 'Exon1End', 'Exon2Start', 'Exon2End']
            self.start_headers(headers, 'p_multi_gene_region.tsv')

        else:
            headers = self.common_headers + ['De Novo', 'Reference Exon Coordinate', 'Exon Spliced With',
                                             'Exon Spliced With Coordinate', 'Junction Name',
                                             'Junction Coordinate'] + self.quantification_headers
            if 'events' in self.config.enabled_outputs:

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
                headers = self.common_headers + ['Junction Coordinate',
                                                 'De Novo',
                                                 'Reference Exon Coordinate',
                                                 'Exon Spliced With Coordinate',
                                                 'Coordinates of Exons Spanned',
                                                 'Number of Exons Spanned'
                                                 ] + self.quantification_headers
                self.start_headers(headers, 'multi_exon_spanning.tsv')
                headers = self.common_headers + ['De Novo', 'Reference Exon Coordinate', 'Exon Spliced With',
                                                 'Exon Spliced With Coordinate', 'Tandem Exon Coordinates',
                                                 'Num_Tandem_Exons',
                                                 'Junction Name', 'Junction Coordinate'] + self.quantification_headers
                self.start_headers(headers, 'tandem_cassette.tsv')
                headers = self.common_headers + ['De Novo', 'Exon coordinate', 'Junction Coordinate'] + self.quantification_headers
                self.start_headers(headers, 'exitron.tsv')

                if self.config.keep_constitutive:
                    headers = self.common_headers + ['De Novo', 'Reference Exon Coordinate', 'Exon Spliced With',
                                                     'Exon Spliced With Coordinate', 'Junction Name',
                                                     'Junction Coordinate', 'Is Intron',
                                                     'Collapsed Event Name'] + self.quantification_headers
                    self.start_headers(headers, 'constitutive.tsv')





            if 'summary' in self.config.enabled_outputs:
                headers = self.common_headers + ["Cassette", "Tandem Cassette",
                                                 "Alt 3", "Alt 5", "P_Alt 3", "P_Alt 5", "Alt 3 and Alt 5", "MXE",
                                                 "Alternative Intron", "ALE", "AFE",
                                                 "P_ALE", "P_AFE", "Orphan Junction"]
                if self.config.keep_constitutive:
                    headers.append("Constitutive Junction")
                    headers.append("Constitutive Intron")

                if self.config.output_complex:
                    headers.remove("Complex")

                headers += ["Multi Exon Spanning", "Exitron", "Complex", 'De Novo Junctions', 'De Novo Introns',
                            "Number of Events", "Collapsed Event Name"
                            ]
                self.start_headers(headers, 'summary.tsv')

            if 'heatmap' in self.config.enabled_outputs:
                headers = self.common_headers + ['Collapsed Event Name'] + self.quantification_headers
                self.start_headers(headers, 'heatmap.tsv')

            if 'junctions' in self.config.enabled_outputs:
                headers = self.common_headers + ['Collapsed Event Name', 'Junction Name',
                                                 'Junction Coordinate', 'De Novo'] + self.quantification_headers
                self.start_headers(headers, 'junctions.tsv')

            if 'mpe' in self.config.enabled_outputs:
                headers = self.common_headers
                lids_col_ii = headers.index("LSV ID(s)")
                # Only 1 LSV ID per row possible.
                headers[lids_col_ii].replace("(s)","")
                headers += ['Collapsed Event Name',
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

                            row = [event['Skip'].de_novo, event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2',
                                   event['Skip'].range_str()]
                            quants = self.quantifications(module, 's', event['Skip'], event['C1'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            row = [event['Include1'].de_novo, event['C1'].range_str(), 'A', event['A'].range_str(), 'C1_A',
                                   event['Include1'].range_str()]
                            quants = self.quantifications(module, 's', event['Include1'], event['C1'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            row = [event['Skip'].de_novo, event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1',
                                   event['Skip'].range_str()]
                            quants = self.quantifications(module, 't', event['Skip'], event['C2'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                            row = [event['Include2'].de_novo, event['C2'].range_str(), 'A', event['A'].range_str(), 'C2_A',
                                   event['Include2'].range_str()]
                            quants = self.quantifications(module, 't', event['Include2'], event['C2'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                            if True:
                                if trg_common[5]:
                                    self.heatmap_add(module, trg_common, self.quantifications(module, 't', event['Include2'], event['C2']),
                                                     event['Include2'].end - event['Include2'].start)
                                elif src_common[5]:
                                    self.heatmap_add(module, src_common,
                                                     self.quantifications(module, 's', event['Include1'], event['C1']),
                                                     event['Include1'].end - event['Include1'].start)

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
                                row = [event['Proximal'].de_novo, event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
                                       event['Proximal'].range_str()]
                                quants = self.quantifications(module, 's', event['Proximal'])
                                writer.writerow(src_common + row + quants)
                                self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                                row = [event['Distal'].de_novo, event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
                                       event['Distal'].range_str()]
                                quants = self.quantifications(module, 's', event['Distal'])
                                writer.writerow(src_common + row + quants)
                                self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            elif trg_common[5]:
                                row = [event['Proximal'].de_novo, event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Proximal',
                                       event['Proximal'].range_str()]
                                quants = self.quantifications(module, 't', event['Proximal'])
                                writer.writerow(trg_common + row + quants)
                                self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                                row = [event['Distal'].de_novo, event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Distal',
                                       event['Distal'].range_str()]
                                quants = self.quantifications(module, 't', event['Distal'])
                                writer.writerow(trg_common + row + quants)
                                self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                            if True:
                                if src_common[5]:
                                    self.heatmap_add(module, src_common, self.quantifications(module, 's', event['Proximal'], event['E1']),
                                                     event['Proximal'].end - event['Proximal'].start)
                                elif trg_common[5]:
                                    self.heatmap_add(module, trg_common, self.quantifications(module, 't', event['Proximal'], event['E2']),
                                                     event['Proximal'].end - event['Proximal'].start)


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
                            if trg_common[5]:
                                row = [event['Proximal'].de_novo, event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Proximal',
                                       event['Proximal'].range_str()]
                                quants = self.quantifications(module, 't', event['Proximal'])
                                writer.writerow(trg_common + row + quants)
                                self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                                row = [event['Distal'].de_novo, event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Distal',
                                       event['Distal'].range_str()]
                                quants = self.quantifications(module, 't', event['Distal'])
                                writer.writerow(trg_common + row + quants)
                                self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            elif src_common[5]:
                                row = [event['Proximal'].de_novo, event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
                                       event['Proximal'].range_str()]
                                quants = self.quantifications(module, 's', event['Proximal'])
                                writer.writerow(src_common + row + quants)
                                self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                                row = [event['Distal'].de_novo, event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
                                       event['Distal'].range_str()]
                                quants = self.quantifications(module, 's', event['Distal'])
                                writer.writerow(src_common + row + quants)
                                self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            if True:
                                if trg_common[5]:
                                    self.heatmap_add(module, trg_common, self.quantifications(module, 't', event['Proximal'], event['E2']),
                                                     event['Proximal'].end - event['Proximal'].start)
                                elif src_common[5]:
                                    self.heatmap_add(module, src_common, self.quantifications(module, 's', event['Proximal'], event['E1']),
                                                     event['Proximal'].end - event['Proximal'].start)



    def p_alt5prime(self):
        with open(os.path.join(self.config.directory, 'p_alt5prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'p_alt5ss':
                            src_common = self.common_data(module, 's', node=event['C1'])
                            trg_common = self.common_data(module, 't', node=event['C2'])

                            row = [event['Skip'].de_novo, event['C1'].range_str(), 'E3', event['C2'].range_str(), 'E1_E3_Distal',
                                   event['Skip'].range_str()]
                            quants = self.quantifications(module, 's', event['Skip'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            event['Include1'].junc['start'] += 1
                            event['Include1'].junc['end'] -= 1

                            row = [event['Include1'].de_novo, event['C1'].range_str(), 'E2', event['A'].range_str(), 'E1_E2_Intron',
                                   event['Include1'].range_str()]
                            quants = self.quantifications(module, 's', event['Include1'])
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1']))
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            row = [event['Skip'].de_novo, event['C2'].range_str(), 'E1', event['C1'].range_str(), 'E3_E1_Distal',
                                   event['Skip'].range_str()]
                            quants = self.quantifications(module, 't', event['Skip'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                            row = [event['Include2'].de_novo, event['C2'].range_str(), 'E2', event['A'].range_str(), 'E3_E2_Proximal',
                                   event['Include2'].range_str()]
                            quants = self.quantifications(module, 't', event['Include2'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                            if True:
                                if trg_common[5]:
                                    self.heatmap_add(module, trg_common, self.quantifications(module, 't', event['Include2'], event['C2']),
                                                 event['Include2'].end - event['Include2'].start)
                                else:
                                    self.heatmap_add(module, src_common,
                                                     self.quantifications(module, 's', event['Include1'], event['C1']),
                                                     event['Include1'].end - event['Include1'].start)

                            event['Include1'].junc['start'] -= 1
                            event['Include1'].junc['end'] += 1


    def p_alt3prime(self):
        with open(os.path.join(self.config.directory, 'p_alt3prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'p_alt3ss':
                            src_common = self.common_data(module, 's', node=event['C1'])
                            trg_common = self.common_data(module, 't', node=event['C2'])

                            event['Include2'].junc['start'] += 1
                            event['Include2'].junc['end'] -= 1

                            row = [event['Skip'].de_novo, event['C2'].range_str(), 'E3', event['C1'].range_str(), 'E1_E3_Distal',
                                   event['Skip'].range_str()]
                            quants = self.quantifications(module, 's', event['Skip'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            row = [event['Include1'].de_novo, event['C2'].range_str(), 'E2', event['A'].range_str(), 'E1_E2_Proximal',
                                   event['Include1'].range_str()]
                            quants = self.quantifications(module, 's', event['Include1'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            row = [event['Include2'].de_novo, event['C1'].range_str(), 'E2', event['A'].range_str(), 'E3_E2_Intron',
                                   event['Include2'].range_str()]
                            quants = self.quantifications(module, 't', event['Include2'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                            row = [event['Skip'].de_novo, event['C1'].range_str(), 'E1', event['C2'].range_str(), 'E3_E1_Distal',
                                   event['Skip'].range_str()]
                            quants = self.quantifications(module, 't', event['Skip'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                            if True:
                                if src_common[5]:
                                    self.heatmap_add(module, src_common, self.quantifications(module, 's', event['Include1'], event['C1']),
                                                 event['Include1'].end - event['Include1'].start)
                                else:
                                    self.heatmap_add(module, trg_common,
                                                     self.quantifications(module, 't', event['Include2'], event['C2']),
                                                     event['Include2'].end - event['Include2'].start)

                            event['Include2'].junc['start'] -= 1
                            event['Include2'].junc['end'] += 1

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
                            row = [event['J1'].de_novo, event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_J1',
                                   event['J1'].range_str()]
                            quants = self.quantifications(module, 's', event['J1'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            row = [event['J2'].de_novo, event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_J2',
                                   event['J2'].range_str()]
                            quants = self.quantifications(module, 's', event['J2'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                            row = [event['J1'].de_novo, event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_J1',
                                   event['J1'].range_str()]
                            quants = self.quantifications(module, 't', event['J1'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            row = [event['J2'].de_novo, event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_J2',
                                   event['J2'].range_str()]
                            quants = self.quantifications(module, 't', event['J2'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                            if True:
                                if trg_common[5]:
                                    self.heatmap_add(module, trg_common, self.quantifications(module, 't', event['J2']),
                                                 event['J2'].end - event['J2'].start)
                                else:
                                    self.heatmap_add(module, src_common, self.quantifications(module, 's', event['J1']),
                                                     event['J1'].end - event['J1'].start)

    def mutually_exclusive(self):
        with open(os.path.join(self.config.directory, 'mutually_exclusive.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'mutually_exclusive':

                            row = [event['Include1'].de_novo, event['C1'].range_str(), 'A1', event['A1'].range_str(), 'C1_A1',
                                   event['Include1'].range_str()]
                            common = self.common_data(module, 's', node=event['C1'], edge=event['Include1'])
                            quants = self.quantifications(module, 's', event['Include1'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))
                            row = [event['SkipA1'].de_novo, event['C1'].range_str(), 'A2', event['A2'].range_str(), 'C1_A2',
                                   event['SkipA1'].range_str()]
                            common = self.common_data(module, 's', node=event['C1'], edge=event['SkipA1'])
                            quants = self.quantifications(module, 's', event['SkipA1'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))
                            row = [event['Include2'].de_novo, event['C2'].range_str(), 'A1', event['A1'].range_str(), 'C2_A1',
                                   event['Include2'].range_str()]
                            common = self.common_data(module, 't', node=event['C2'], edge=event['Include2'])
                            quants = self.quantifications(module, 't', event['Include2'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))
                            row = [event['SkipA2'].de_novo, event['C2'].range_str(), 'A2', event['A2'].range_str(), 'C2_A2',
                                   event['SkipA2'].range_str()]
                            common = self.common_data(module, 't', node=event['C2'], edge=event['SkipA2'])
                            quants = self.quantifications(module, 't', event['SkipA2'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))

                            if True:
                                common = self.common_data(module, 't', node=event['C2'], edge=event['Include2'])
                                self.heatmap_add(module, common, self.quantifications(module, 't', event['Include2']),
                                                 event['Include2'].end - event['Include2'].start)

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
                                    row = [junc.de_novo, event['Reference'].range_str(), 'A_Proximal', event['Proximal'].range_str(),
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    quants = self.quantifications(module, 's', junc)
                                    writer.writerow(src_common + row + quants)
                                    self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                                for junc in event['SkipA1']:
                                    row = [junc.de_novo, event['Reference'].range_str(), 'A_Distal', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    quants = self.quantifications(module, 's', junc)
                                    writer.writerow(src_common + row + quants)
                                    self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                                if True:
                                    if event['SkipA2']:
                                        self.heatmap_add(module, src_common,
                                                         self.quantifications(module, 's', event['SkipA2'][0], event['Reference']),
                                                         event['SkipA2'][0].end - event['SkipA2'][0].start)
                                    elif event['SkipA1']:
                                        self.heatmap_add(module, src_common,
                                                         self.quantifications(module, 's', event['SkipA1'][0], event['Reference']),
                                                         event['SkipA1'][0].end - event['SkipA1'][0].start)

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
                                    row = [junc.de_novo, event['Reference'].range_str(), 'A_Proximal', event['Proximal'].range_str(),
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    quants = self.quantifications(module, 't', junc)
                                    writer.writerow(trg_common + row + quants)
                                    self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                                for junc in event['SkipA2']:
                                    row = [junc.de_novo, event['Reference'].range_str(), 'A_Distal', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    quants = self.quantifications(module, 't', junc)
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                                    self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                                if True:
                                    if event['SkipA1']:
                                        self.heatmap_add(module, trg_common,
                                                         self.quantifications(module, 't', event['SkipA1'][0], event['Reference']),
                                                         event['SkipA1'][0].end - event['SkipA1'][0].start)
                                    elif event['SkipA2']:
                                        self.heatmap_add(module, trg_common,
                                                         self.quantifications(module, 't', event['SkipA2'][0], event['Reference']),
                                                         event['SkipA2'][0].end - event['SkipA2'][0].start)



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
                                    row = [junc.de_novo, event['Reference'].range_str(), 'A', proxStr,
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    quants = self.quantifications(module, 's', junc)
                                    writer.writerow(src_common + row + quants)
                                    self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                                for junc in event['SkipA1']:
                                    row = [junc.de_novo, event['Reference'].range_str(), 'A', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    quants = self.quantifications(module, 's', junc)
                                    writer.writerow(src_common + row + quants)
                                    self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                                if True:
                                    if event['SkipA2']:
                                        self.heatmap_add(module, src_common,
                                                         self.quantifications(module, 's', event['SkipA2'][0], event['Reference']),
                                                         event['SkipA2'][0].end - event['SkipA2'][0].start)
                                    elif event['SkipA1']:
                                        self.heatmap_add(module, src_common,
                                                         self.quantifications(module, 's', event['SkipA1'][0], event['Reference']),
                                                         event['SkipA1'][0].end - event['SkipA1'][0].start)

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
                                    row = [junc.de_novo, event['Reference'].range_str(), 'A', proxStr,
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    quants = self.quantifications(module, 't', junc)
                                    writer.writerow(trg_common + row + quants)
                                    self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                                for junc in event['SkipA2']:
                                    row = [junc.de_novo, event['Reference'].range_str(), 'A', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    quants = self.quantifications(module, 't', junc)
                                    writer.writerow(trg_common + row + quants)
                                    self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))

                                if True:
                                    if event['SkipA1']:
                                        self.heatmap_add(module, trg_common,
                                                         self.quantifications(module, 't', event['SkipA1'][0], event['Reference']),
                                                         event['SkipA1'][0].end - event['SkipA1'][0].start)
                                    elif event['SkipA2']:
                                        self.heatmap_add(module, trg_common,
                                                         self.quantifications(module, 't', event['SkipA2'][0], event['Reference']),
                                                         event['SkipA2'][0].end - event['SkipA2'][0].start)

    def alternative_intron(self):
        with open(os.path.join(self.config.directory, 'alternative_intron.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'alternative_intron':

                            # put coordinates back to Jordi's offset numbers
                            event['Intron'].junc['start'] += 1
                            event['Intron'].junc['end'] -= 1

                            src_common = self.common_data(module, 's', event['Intron'])
                            trg_common = self.common_data(module, 't', event['Intron'])


                            if any(':t:' in _l for _l in event['Intron'].lsvs) and not \
                               any(':s:' in _l for _l in event['Intron'].lsvs):
                                row = [event['Intron'].de_novo, event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1_intron',
                                       event['Intron'].range_str()]
                                quants = self.quantifications(module, 't', event['Intron'])
                                writer.writerow(trg_common + row + quants)
                                self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                                row = [self.semicolon((x.de_novo for x in event['Spliced'])),
                                       event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1_spliced',
                                       self.semicolon((x.range_str() for x in event['Spliced']))]
                                quants = self.quantifications(module, 't', event['Spliced'])
                                writer.writerow(trg_common + row + quants)
                                self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))



                            else:
                                row = [event['Intron'].de_novo,
                                       event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2_intron',
                                       event['Intron'].range_str()]
                                quants = self.quantifications(module, 's', event['Intron'])
                                writer.writerow(src_common + row + quants)
                                self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))
                                row = [self.semicolon((x.de_novo for x in event['Spliced'])),
                                       event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2_spliced',
                                       self.semicolon((x.range_str() for x in event['Spliced']))]
                                quants = self.quantifications(module, 's', event['Spliced'])
                                writer.writerow(src_common + row + quants)
                                self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            if True:
                                if trg_common[5]:
                                    self.heatmap_add(module, trg_common,
                                                     self.quantifications(module, 't', event['Intron']),
                                                     event['Intron'].end - event['Intron'].start)

                            if True:
                                if src_common[5]:
                                    self.heatmap_add(module, src_common,
                                                     self.quantifications(module, 's', event['Intron']),
                                                     event['Intron'].end - event['Intron'].start)

                            event['Intron'].junc['start'] -= 1
                            event['Intron'].junc['end'] += 1


    def multi_exon_spanning(self):
        with open(os.path.join(self.config.directory, 'multi_exon_spanning.tsv.%s' % self.pid), 'a',
                  newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'multi_exon_spanning':
                            #src_common = self.common_data(module, 's')
                            # Source LSV side
                            row = [self.semicolon((x.range_str() for x in event['Skip'])),  # junction coord
                                    self.semicolon((x.de_novo for x in event['Skip'])), # de novo?
                                    event['C1'].range_str(), # reference exon
                                    event['C2'].range_str(), # exon spliced with
                                    self.semicolon((x.range_str() for x in event['As'])), # exons spanned
                                    len(event['As'])] # num exons spanned
                            quants = self.quantifications(module, 's', event['Skip'], event['C1'])
                            common = self.common_data(module, 's', node=event['C1'], edge=event['Skip'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[1], '', row[0]))
                            # Target LSV side
                            row = [self.semicolon((x.range_str() for x in event['Skip'])),  # junction coord
                                   self.semicolon((x.de_novo for x in event['Skip'])),  # de novo?
                                   event['C2'].range_str(),  # reference exon
                                   event['C1'].range_str(),  # exon spliced with
                                   self.semicolon((x.range_str() for x in event['As'])),  # exons spanned
                                   len(event['As'])]  # num exons spanned
                            quants = self.quantifications(module, 't', event['Skip'], event['C2'])
                            common = self.common_data(module, 't', node=event['C2'], edge=event['Skip'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[1], '', row[0]))


    def tandem_cassette(self):
        with open(os.path.join(self.config.directory, 'tandem_cassette.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'tandem_cassette':
                            # if event['Skip'][0].range_str() == "65976202-65980224":
                            #     print(event['Skip'][0].range_str())
                            #     exit(1)
                            row = [event['Skip'].de_novo,
                                   event['C1'].range_str(),
                                   'C2',
                                   event['C2'].range_str(),
                                   self.semicolon((x.range_str() for x in event['Tandem_Exons'])),
                                   len(event['Tandem_Exons']),
                                   'C1_C2',
                                   event['Skip'].range_str()]
                            common = self.common_data(module, 's', node=event['C1'], edge=event['Skip'])
                            quants = self.quantifications(module, 's', event['Skip'], event['C1'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))
                            row = [event['Include1'].de_novo,
                                   event['C1'].range_str(),
                                   'A1',
                                   event['Tandem_Exons'][0].range_str(),
                                   self.semicolon((x.range_str() for x in event['Tandem_Exons'])),
                                   len(event['Tandem_Exons']),
                                   'C1_A',
                                   event['Include1'].range_str()]
                            common = self.common_data(module, 's', node=event['C1'], edge=event['Include1'])
                            quants = self.quantifications(module, 's', event['Include1'], event['C1'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))
                            row = [event['Skip'].de_novo,
                                   event['C2'].range_str(),
                                   'C1',
                                   event['C1'].range_str(),
                                   self.semicolon((x.range_str() for x in event['Tandem_Exons'])),
                                   len(event['Tandem_Exons']),
                                   'C2_C1',
                                   event['Skip'].range_str()]
                            common = self.common_data(module, 't', node=event['C2'], edge=event['Skip'])
                            quants = self.quantifications(module, 't', event['Skip'], event['C2'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))
                            row = [event['Include2'].de_novo,
                                   event['C2'].range_str(),
                                   'A_Last',
                                   event['Tandem_Exons'][-1].range_str(),
                                   self.semicolon((x.range_str() for x in event['Tandem_Exons'])),
                                   len(event['Tandem_Exons']),
                                   'C2_A_Last',
                                   event['Include2'].range_str()]
                            common = self.common_data(module, 't', node=event['C2'], edge=event['Include2'])
                            quants = self.quantifications(module, 't', event['Include2'], event['C2'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))

                            if True:
                                common = self.common_data(module, 't', node=event['C2'], edge=event['Include2'])
                                self.heatmap_add(module, common, self.quantifications(module, 't', event['Include2'], event['C2']),
                                                 event['Include2'].end - event['Include2'].start)

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
                            row = [event['Junc'].de_novo, event['Exon'].range_str(), event['Junc'].range_str()]
                            # just in case exitrons are ever quantified, somehow, *try* to get
                            # the exitron's junction quantification (won't exist for now... which is desired)
                            quants = self.quantifications(module, 's', edge=event['Junc'])
                            writer.writerow(src_common + row + quants)
                            self.junction_cache.append((module, src_common, quants, row[0], 'Exitron', row[2]))
                            row = [event['Junc'].de_novo, event['Exon'].range_str(), event['Junc'].range_str()]
                            quants = self.quantifications(module, 't', edge=event['Junc'])
                            writer.writerow(trg_common + row + quants)
                            self.junction_cache.append((module, trg_common, quants, row[0], 'Exitron', row[2]))

    def constitutive(self):
        with open(os.path.join(self.config.directory, 'constitutive.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _total_events = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'constitutive':
                            common = self.common_data(module)
                            row = [event['Junc'].de_novo,
                                   event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1',
                                   event['Junc'].range_str(), 'False', module.collapsed_event_name]
                            quants = self.quantifications(module, edge=event['Junc'])
                            writer.writerow(common + row + quants)
                            self.junction_cache.append((module, common, quants, row[0], row[4], row[5]))

                        elif event['event'] == 'constitutive_intron':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')

                            # put coordinates back to Jordi's offset numbers
                            event['Intron'].junc['start'] += 1
                            event['Intron'].junc['end'] -= 1
                            if any(':t:' in _l for _l in event['Intron'].lsvs) and not \
                               any(':s:' in _l for _l in event['Intron'].lsvs):
                                row = [event['Intron'].de_novo,
                                       event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1_intron',
                                       event['Intron'].range_str(), 'True', module.collapsed_event_name]
                                quants = self.quantifications(module, 't', event['Intron'])
                                writer.writerow(trg_common + row + quants)
                                self.junction_cache.append((module, trg_common, quants, row[0], row[4], row[5]))
                            else:
                                row = [event['Intron'].de_novo,
                                       event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2_intron',
                                       event['Intron'].range_str(), 'True', module.collapsed_event_name]
                                quants = self.quantifications(module, 's', event['Intron'])
                                writer.writerow(src_common + row + quants)
                                self.junction_cache.append((module, src_common, quants, row[0], row[4], row[5]))

                            event['Intron'].junc['start'] -= 1
                            event['Intron'].junc['end'] += 1


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

    def heatmap_add(self, module, common, quants, junc_len):
        """
        Conditionally add a row toe heatmap cache by comparing it to what exists there already
        """

        if not module.idx in self.heatmap_cache:
            self.heatmap_cache[module.idx] = (module, common, quants, junc_len)
        else:
            if self.heatmap_cache[module.idx][3] > junc_len:
                self.heatmap_cache[module.idx] = (module, common, quants, junc_len)

    def junctions(self):
        """
        Write a file with a listing of all junctions
        :return:
        """
        with open(os.path.join(self.config.directory, 'junctions.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module, common_data, quantifications, de_novo, junction_name, coordinates in self.junction_cache:
                writer.writerow(common_data + [module.collapsed_event_name, junction_name, coordinates, de_novo] + quantifications)

    def mpe(self):
        """
        Write a file with a listing of all primerable regions
        :return:
        """
        with open(os.path.join(self.config.directory, 'mpe_primerable_regions.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events = self.mpe_regions[module.idx]
                for event in events:
                    n_assigned_lids = 0
                    lsvid = ""
                    edge_type = ""
                    if event['event'] == "mpe_source":
                        common = self.common_data(module, 's', event['reference_exon'])
                        eventtype = "Single Source" # SingleSource
                        constitutive_direction = "Upstream"
                        if event['edge_of_transcript']:
                            edge_type = "first_exon"
                    else:
                        common = self.common_data(module, 't', event['reference_exon'])
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

    def heatmap(self):
        """
        Write the easily excel-able file
        This is similar to the other non-summary event files and also part of the summary file
        It lists rows like the non-summary files, but only chooses the shortest junction
        """
        with open(os.path.join(self.config.directory, 'heatmap.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module, common_data, quantifications, junc_len in self.heatmap_cache.values():
                writer.writerow(common_data + [module.collapsed_event_name] + quantifications)

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
                de_novo_junctions = 0
                de_novo_introns = 0
                for edge in module.get_all_edges(ir=True):
                    if edge.de_novo:
                        if edge.ir:
                            de_novo_introns += 1
                        else:
                            de_novo_junctions += 1

                # we store collapsed event name on module, because we need it for constitutive
                module.collapsed_event_name = self._collapsed_event_name(counts)

                writer.writerow(["%s_%d" % (self.gene_id, module.idx),
                                 self.gene_id, self.graph.gene_name, self.graph.chromosome, self.graph.strand,
                                 self.semicolon(module.target_lsv_ids.union(module.source_lsv_ids))] +
                                [v if v else '' for v in counts.values()] +
                                [str(_complex), str(de_novo_junctions), str(de_novo_introns),
                                 str(_total_events), module.collapsed_event_name]
                                )



