import csv
import os
from voila.api import Matrix
from voila import constants
from voila.exceptions import GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile
from voila.api.matrix_utils import generate_variances
from collections import OrderedDict
from voila.config import ClassifyConfig
import multiprocessing

def semicolon(value_list):
    return ';'.join(str(x) for x in value_list)


class TsvWriter:
    """
    Output AS data from one gene
    """

    def __init__(self, graph, gene_id, quantifications=('psi', 'var',)):
        """

        :param output_path: The folder where all output TSV files will be written under
        :param graph: the Graph object of the gene
        """
        
        self.common_headers = ['Module ID', 'LSV ID(s)', 'Gene ID', 'Gene Name', 'Chr', 'Strand']
        self.graph = graph
        self.gene_id = gene_id
        self.quantifications_enabled = quantifications
        self.config = ClassifyConfig()
        self.pid = multiprocessing.current_process().pid

        # we could do some crazy thing to yield to all of the different output types at once (across each method)
        # (in order to save memory) But for now we just save modules in a list. Will ammend later if memory use
        # becomes an issue.
        if self.graph:
            self.modules = self.graph.modules()

            #self.as_types = {x.idx: x.as_types() for x in self.modules}
            self.as_types = {x.idx: x.as_types() for x in self.modules}

    @property
    def quantification_headers(self):
        headers = []
        for input_file in self.config.voila_files:
            trunc_name = str(input_file).split('/')[-1].split('.')[0]
            if 'psi' in self.quantifications_enabled:
                headers.append('%s_E(PSI)' % trunc_name)
            if 'var' in self.quantifications_enabled:
                headers.append('%s_Var(E(PSI))' % trunc_name)
        return headers

    @staticmethod
    def tsv_names():
        return ('summary.tsv', 'cassette.tsv', 'alt3prime.tsv', 'alt5prime.tsv', 'alt3and5prime.tsv',
                'mutually_exclusive.tsv', 'alternate_last_exon.tsv', 'alternate_first_exon.tsv',
                'intron_retention.tsv', 'p_alt5prime.tsv', 'p_alt3prime.tsv', 'multi_exon_spanning.tsv',
                'tandem_cassette.tsv', 'exitron.tsv', 'p_multi_gene_region.tsv')

    @staticmethod
    def delete_tsvs():
        for tsv_file in TsvWriter.tsv_names():
            config = ClassifyConfig()
            path = os.path.join(config.directory, tsv_file)
            if os.path.exists(path):
                os.remove(path)

    @staticmethod
    def parity2lsv(module, parity):
        if parity == 's':
            lsvs = module.source_lsv_ids
        elif parity == 't':
            lsvs = module.target_lsv_ids
        else:
            lsvs = module.target_lsv_ids.union(module.source_lsv_ids)
        return lsvs

    def common_data(self, module, parity=None, edge=None):
        """
        Extract the certain cols from the CSV which are generally similar across all outputs,

        """
        lsvs = self.parity2lsv(module, parity)


        return ["%s_%d" % (self.gene_id, module.idx), semicolon(lsvs), self.gene_id, self.graph.gene_name,
                self.graph.chromosome, self.graph.strand]

    def quantifications(self, module, parity=None, edge=None):

        lsvs = self.parity2lsv(module, parity)

        with Matrix(self.config.voila_file) as m:
            analysis_type = m.analysis_type
        quantification_fields = []
        for i, voila_file in enumerate(self.config.voila_files):
            try:
                with Matrix(voila_file) as m:
                    means = []
                    vars = []
                    for lsv_id in lsvs:
                        if analysis_type == constants.ANALYSIS_PSI:
                            lsv = m.psi(lsv_id)
                        else:
                            lsv = m.delta_psi(lsv_id)

                        if edge:
                            # loop through junctions to find one matching range of edge
                            for j, junc in enumerate(lsv.get('junctions')):
                                if junc[0] == edge.start and junc[1] == edge.end:
                                    means.append(lsv.get('means')[j])
                                    vars.append(generate_variances([lsv.get('bins')[i]])[0])
                                    break
                            else:
                                # junction not quantified by majiq
                                pass
                        else:
                            means += list(lsv.get('means'))
                            vars += list(generate_variances(lsv.get('bins')))
                if 'psi' in self.quantifications_enabled:
                    quantification_fields.append(semicolon(means))
                if 'var' in self.quantifications_enabled:
                    quantification_fields.append(semicolon(vars))
            except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                if 'psi' in self.quantifications_enabled:
                    quantification_fields.append('')
                if 'var' in self.quantifications_enabled:
                    quantification_fields.append('')

        return quantification_fields


    def start_headers(self, headers, filename):
        """
        Start a tsv file with the required headers, only if it does not yet exist

        """
        if not os.path.exists(os.path.join(self.config.directory, filename)):
            with open(os.path.join(self.config.directory, filename), 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
                writer.writerow(headers)

    def start_all_headers(self):
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
        self.start_headers(headers, 'intron_retention.tsv')
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Tandem Exon Coordinates', 'Num_Tandem_Exons',
                                         'Junction Name', 'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'multi_exon_spanning.tsv')
        self.start_headers(headers, 'tandem_cassette.tsv')
        headers = self.common_headers + ['Exon coordinate', 'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'exitron.tsv')
        headers = ['Module', 'Gene ID', 'Gene Name', 'LSV ID(s)', 'Cassette', 'Alt 3',
                   'Alt 5', 'P_Alt 3', 'P_Alt 5', 'Alt 3 and Alt 5', 'MXE', 'ALE',
                   'AFE', 'P_ALE', 'P_AFE', 'Orphan Junction', 'Multi Exon Spanning',
                   'Tandem Cassette', 'Intron Retention', 'Exitron', 'Complex', 'Multi-Event']
        self.start_headers(headers, 'summary.tsv')
        headers = ['Gene ID_Region', 'Gene ID', 'Gene Name', 'Chr', 'Strand', 'First Exon Start coord',
                   'First Exon End coord', 'Last Exon Start coord', "Last Exon End coord"]
        self.start_headers(headers, 'p_multi_gene_region.tsv')


    def cassette(self):
        with open(os.path.join(self.config.directory, 'cassette.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'cassette_exon':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')

                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2',
                                   event['Skip'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Skip']))

                            row = [event['C1'].range_str(), 'A', event['A'].range_str(), 'C1_A',
                                   event['Include1'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1']))

                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1',
                                   event['Skip'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Skip']))

                            row = [event['C2'].range_str(), 'A', event['A'].range_str(), 'C2_A',
                                   event['Include2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2']))

    def alt3prime(self):
        with open(os.path.join(self.config.directory, 'alt3prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'alt3ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            if src_common[1]:
                                row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
                                       event['Proximal'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Proximal']))
                                row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
                                       event['Distal'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Distal']))
                            if trg_common[1]:
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
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'alt5ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            if src_common[1]:
                                row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
                                       event['Proximal'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Proximal']))
                                row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
                                       event['Distal'].range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', event['Distal']))
                            if trg_common[1]:
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
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'p_alt5ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2',
                                   event['Skip'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Skip']))

                            row = [event['C1'].range_str(), 'A', event['A'].range_str(), 'C1_A',
                                   event['Include1'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1']))

                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1',
                                   event['Skip'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Skip']))

                            row = [event['C2'].range_str(), 'A', event['A'].range_str(), 'C2_A',
                                   event['Include2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2']))

    def p_alt3prime(self):
        with open(os.path.join(self.config.directory, 'p_alt3prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'p_alt3ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2',
                                   event['Skip'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Skip']))

                            row = [event['C1'].range_str(), 'A', event['A'].range_str(), 'C1_A',
                                   event['Include1'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1']))

                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1',
                                   event['Skip'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Skip']))

                            row = [event['C2'].range_str(), 'A', event['A'].range_str(), 'C2_A',
                                   event['Include2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2']))

    def alt3and5prime(self):
        with open(os.path.join(self.config.directory, 'alt3and5prime.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
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
        with open(os.path.join(self.config.directory, 'mutually_exclusive.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
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
                            row = [event['C2'].range_str(), 'A2', event['A2'].range_str(), 'C2_A2',
                                   event['Include2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2']))
                            row = [event['C2'].range_str(), 'A1', event['A1'].range_str(), 'C2_A1',
                                   event['SkipA2'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['SkipA2']))

    # def p_alternate_first_exon(self):
    #     headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
    #                                      'Exon Spliced With Coordinate', 'Junction Name',
    #                                      'Junction Coordinate'] + self.quantification_headers
    #     self.start_headers(headers, 'p_alternate_first_exon.tsv.%s' % self.pid)
    #     with open(os.path.join(self.config.directory, 'p_alternate_first_exon.tsv.%s' % self.pid), 'a', newline='') as csvfile:
    #         writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
    #         for module in self.modules:
    #             events, _complex, _multi_event = self.as_types[module.idx]
    #             if not _complex or self.config.output_complex:
    #                 for event in events:
    #                     if event['event'] == 'p_afe':
    #                         src_common = self.common_data(module, 's')
    #                         trg_common = self.common_data(module, 't')
    #                         for junc in event['SkipA2']:
    #                             row = [event['C1'].range_str(), 'A1', event['A1'].range_str(), 'C1_A1',
    #                                    junc.range_str()]
    #                             writer.writerow(src_common + row + self.quantifications(module, 's', junc))
    #                         for junc in event['SkipA1']:
    #                             row = [event['C1'].range_str(), 'A2', event['A2'].range_str(), 'C1_A2',
    #                                    junc.range_str()]
    #                             writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
    #                     elif event['event'] == 'p_ale':
    #                         trg_common = self.common_data(module, 't')
    #                         row = ['N/A', 'A1', event['A1'].range_str(), 'C1_A1',
    #                                'N/A']
    #                         writer.writerow(trg_common + row + self.quantifications(module, 't', event['A1']))
    #
    # def p_alternate_last_exon(self):
    #     headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
    #                                      'Exon Spliced With Coordinate', 'Junction Name',
    #                                      'Junction Coordinate'] + self.quantification_headers
    #     self.start_headers(headers, 'p_alternate_last_exon.tsv.%s' % self.pid)
    #     with open(os.path.join(self.config.directory, 'p_alternate_last_exon.tsv.%s' % self.pid), 'a', newline='') as csvfile:
    #         writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
    #         for module in self.modules:
    #             events, _complex, _multi_event = self.as_types[module.idx]
    #             if not _complex or self.config.output_complex:
    #                 for event in events:
    #                     if event['event'] == 'p_ale':
    #                         src_common = self.common_data(module, 's')
    #                         trg_common = self.common_data(module, 't')
    #                         for junc in event['SkipA2']:
    #                             row = [event['C1'].range_str(), 'A1', event['A1'].range_str(), 'C1_A1',
    #                                    junc.range_str()]
    #                             writer.writerow(src_common + row + self.quantifications(module, 's', junc))
    #                         for junc in event['SkipA1']:
    #                             row = [event['C1'].range_str(), 'A2', event['A2'].range_str(), 'C1_A2',
    #                                    junc.range_str()]
    #                             writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
    #                     elif event['event'] == 'p_ale':
    #                         trg_common = self.common_data(module, 't')
    #                         row = ['N/A', 'A1', event['A1'].range_str(), 'C1_A1',
    #                                'N/A']
    #                         writer.writerow(trg_common + row + self.quantifications(module, 't', event['A1']))

    def alternate_last_exon(self):
        with open(os.path.join(self.config.directory, 'alternate_last_exon.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'ale':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            if src_common[1]:
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
                            if trg_common[1]:
                                for junc in event['SkipA2']:
                                    row = [event['Reference'].range_str(), 'A', event['Proximal'].range_str(),
                                           'C_A_Proximal',
                                           junc.range_str()]
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                                for junc in event['SkipA1']:
                                    row = [event['Reference'].range_str(), 'A', event['Distal'].range_str(),
                                           'C_A_Distal',
                                           junc.range_str()]
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                        elif event['event'] == 'p_ale':
                            trg_common = self.common_data(module, 't')
                            row = ['N/A', 'A1', event['A1'].range_str(), 'C1_A1',
                                   'N/A']
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['A1']))

    def alternate_first_exon(self):
        with open(os.path.join(self.config.directory, 'alternate_first_exon.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'afe':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            if src_common[1]:
                                for junc in event['SkipA2']:
                                    row = [event['Reference'].range_str(), 'A', event['Proximal'].range_str(), 'C_A_Proximal',
                                           junc.range_str()]
                                    writer.writerow(src_common + row + self.quantifications(module, 's', junc))
                                for junc in event['SkipA1']:
                                    row = [event['Reference'].range_str(), 'A', event['Distal'].range_str(), 'C_A_Distal',
                                           junc.range_str()]
                                    writer.writerow(src_common + row + self.quantifications(module, 's', junc))
                            if trg_common[1]:
                                for junc in event['SkipA2']:
                                    row = [event['Reference'].range_str(), 'A', event['Proximal'].range_str(), 'C_A_Proximal',
                                           junc.range_str()]
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                                for junc in event['SkipA1']:
                                    row = [event['Reference'].range_str(), 'A', event['Distal'].range_str(), 'C_A_Distal',
                                           junc.range_str()]
                                    writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                        elif event['event'] == 'p_afe':
                            trg_common = self.common_data(module, 't')
                            row = ['N/A', 'A1', event['A1'].range_str(), 'C1_A1',
                                   'N/A']
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['A1']))


    def intron_retention(self):
        with open(os.path.join(self.config.directory, 'intron_retention.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'intron_retention':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')

                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2_intron',
                                   event['Intron'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Intron']))
                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2_spliced',
                                   event['Intron'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Intron']))
                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1_intron',
                                   event['Intron'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Intron']))
                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1_spliced',
                                   event['Intron'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Intron']))

    def multi_exon_spanning(self):
        with open(os.path.join(self.config.directory, 'multi_exon_spanning.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
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
                            writer.writerow(src_common + row + self.quantifications(module, 's'))
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
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'tandem_cassette':

                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C1_C2',
                                   semicolon((x.range_str() for x in event['Skip']))]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Skip'][0]))
                            row = [event['C1'].range_str(), 'A1', event['As'][0].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C1_A',
                                   semicolon((x.range_str() for x in event['Include1']))]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Include1'][0]))
                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'C2_C1',
                                   semicolon((x.range_str() for x in event['Skip']))]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Skip'][0]))
                            row = [event['C2'].range_str(), 'A_Last', '',
                                   semicolon((x.range_str() for x in event['As'])), len(event['As']), 'A_Last_C2',
                                   semicolon((x.range_str() for x in event['Include2']))]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Include2'][0]))

    def exitron(self):
        with open(os.path.join(self.config.directory, 'exitron.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
                if not _complex or self.config.output_complex:
                    for event in events:
                        if event['event'] == 'exitron':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['Exon'].range_str(), event['Junc'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's'))
                            row = [event['Exon'].range_str(), event['Junc'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't'))

    def p_multi_gene_region(self):
        with open(os.path.join(self.config.directory, 'p_multi_gene_region.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex, _multi_event = self.as_types[module.idx]
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
                events, _complex, _multi_event = self.as_types[module.idx]
                counts = OrderedDict()
                counts['cassette_exon'] = 0
                counts['alt3ss'] = 0
                counts['alt5ss'] = 0
                counts['p_alt3ss'] = 0
                counts['p_alt5ss'] = 0
                counts['alt3and5ss'] = 0
                counts['mutually_exclusive'] = 0
                counts['ale'] = 0
                counts['afe'] = 0
                counts['p_ale'] = 0
                counts['p_afe'] = 0
                counts['orphan_junction'] = 0
                counts['multi_exon_spanning'] = 0
                counts['tandem_cassette'] = 0
                counts['intron_retention'] = 0
                counts['exitron'] = 0
                for event in events:
                    if event['event'] in counts:
                        counts[event['event']] += 1



                writer.writerow(["%s_%d" % (self.gene_id, module.idx),
                                 self.gene_id, self.graph.gene_name,
                                 semicolon(module.target_lsv_ids.union(module.source_lsv_ids))] +
                                [v if v else '' for v in counts.values()] + [str(_complex), str(_multi_event)]
                                )
