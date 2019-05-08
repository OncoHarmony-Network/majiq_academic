import csv
import os
from voila.api import Matrix
from voila import constants
from voila.exceptions import GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile
from voila.api.matrix_utils import generate_variances
from collections import OrderedDict
from voila.config import ClassifyConfig

def semicolon(value_list):
    return ';'.join(str(x) for x in value_list)

SHOW_COMPLEX_IN_ALL = False


class TsvWriter:
    """
    Output AS data from one gene
    """

    def __init__(self, graph, gene_id, quantifications=('psi', 'var',)):
        """

        :param output_path: The folder where all output TSV files will be written under
        :param graph: the Graph object of the gene
        """
        
        self.common_headers = ['Module ID', 'LSV ID(s)', 'Gene Name', 'Gene ID', 'Chr', 'Strand']
        self.graph = graph
        self.gene_id = gene_id
        self.quantifications_enabled = quantifications
        self.config = ClassifyConfig()

        # we could do some crazy thing to yield to all of the different output types at once (across each method)
        # (in order to save memory) But for now we just save modules in a list. Will ammend later if memory use
        # becomes an issue.
        self.modules = [x for x in self.graph.modules()]

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
    def delete_tsvs():
        for tsv_file in ('summary.tsv', 'cassette.tsv', 'alt3prime.tsv', 'alt5prime.tsv', 'alt3and5prime.tsv',
                         'mutually_exclusive.tsv', 'alternate_last_exon.tsv', 'alternate_fist_exon.tsv',
                         'intron_retention.tsv'):
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

    def cassette(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                                       'Exon Spliced With Coordinate', 'Junction Name',
                                                       'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'cassette.tsv')
        with open(os.path.join(self.config.directory, 'cassette.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:

                events, _complex = self.as_types[module.idx]
                if not _complex or SHOW_COMPLEX_IN_ALL:
                    for event in events:
                        if event['event'] == 'cassette_exon':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            for junc in event['Skip']:
                                row = [event['C1'].range_str(), 'C2', event['C2'].range_str(), 'C1_C2',
                                       junc.range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', junc))
                            for junc in event['Include1']:
                                row = [event['C1'].range_str(), 'A', event['A'].range_str(), 'C1_A',
                                       junc.range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', junc))
                            for junc in event['Skip']:
                                row = [event['C2'].range_str(), 'C1', event['C1'].range_str(), 'C2_C1',
                                       junc.range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                            for junc in event['Include2']:
                                row = [event['C2'].range_str(), 'A', event['A'].range_str(), 'C2_A',
                                       junc.range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', junc))

    def alt3prime(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Junction Name',
                                         'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'alt3prime.tsv')
        with open(os.path.join(self.config.directory, 'alt3prime.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex = self.as_types[module.idx]
                if not _complex or SHOW_COMPLEX_IN_ALL:
                    for event in events:
                        if event['event'] == 'alt3ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
                                   event['Proximal'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Proximal']))
                            row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
                                   event['Distal'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Distal']))
                            row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Proximal',
                                   event['Proximal'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Proximal']))
                            row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Distal',
                                   event['Distal'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Distal']))
                            
    def alt5prime(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Junction Name',
                                         'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'alt5prime.tsv')
        with open(os.path.join(self.config.directory, 'alt5prime.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex = self.as_types[module.idx]
                if not _complex or SHOW_COMPLEX_IN_ALL:
                    for event in events:
                        if event['event'] == 'alt5ss':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
                                   event['Proximal'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Proximal']))
                            row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
                                   event['Distal'].range_str()]
                            writer.writerow(src_common + row + self.quantifications(module, 's', event['Distal']))
                            row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Proximal',
                                   event['Proximal'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Proximal']))
                            row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Distal',
                                   event['Distal'].range_str()]
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['Distal']))

    def alt3and5prime(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Junction Name',
                                         'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'alt3and5prime.tsv')
        with open(os.path.join(self.config.directory, 'alt3and5prime.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            # for module in self.modules:
            #     events, _complex = self.as_types[module.idx]
            #     if not _complex or SHOW_COMPLEX_IN_ALL:
            #         for event in events:
            #             if event['event'] == 'alt5ss':
            #                 src_common = self.common_data(module, 's')
            #                 trg_common = self.common_data(module, 't')
            #                 row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Proximal',
            #                        event['Proximal'].range_str()]
            #                 writer.writerow(src_common[0] + row + src_common[1])
            #                 row = [event['E1'].range_str(), 'E2', event['E2'].range_str(), 'E1_E2_Distal',
            #                        event['Distal'].range_str()]
            #                 writer.writerow(src_common[0] + row + src_common[1])
            #                 row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Proximal',
            #                        event['Proximal'].range_str()]
            #                 writer.writerow(trg_common[0] + row + trg_common[1])
            #                 row = [event['E2'].range_str(), 'E1', event['E1'].range_str(), 'E2_E1_Distal',
            #                        event['Distal'].range_str()]
            #                 writer.writerow(trg_common[0] + row + trg_common[1])

    def mutually_exclusive(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Junction Name',
                                         'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'mutually_exclusive.tsv')
        with open(os.path.join(self.config.directory, 'mutually_exclusive.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex = self.as_types[module.idx]
                if not _complex or SHOW_COMPLEX_IN_ALL:
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

    def alternate_last_exon(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Junction Name',
                                         'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'alternate_last_exon.tsv')
        with open(os.path.join(self.config.directory, 'alternate_last_exon.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex = self.as_types[module.idx]
                if not _complex or SHOW_COMPLEX_IN_ALL:
                    for event in events:
                        if event['event'] == 'ale':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            for junc in event['SkipA2']:
                                row = [event['C1'].range_str(), 'A1', event['A1'].range_str(), 'C1_A1',
                                       junc.range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', junc))
                            for junc in event['SkipA1']:
                                row = [event['C1'].range_str(), 'A2', event['A2'].range_str(), 'C1_A2',
                                       junc.range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                        elif event['event'] == 'p_ale':
                            trg_common = self.common_data(module, 't')
                            row = ['N/A', 'A1', event['A1'].range_str(), 'C1_A1',
                                   'N/A']
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['A1']))

    def alternate_first_exon(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Junction Name',
                                         'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'alternate_first_exon.tsv')
        with open(os.path.join(self.config.directory, 'alternate_first_exon.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex = self.as_types[module.idx]
                if not _complex or SHOW_COMPLEX_IN_ALL:
                    for event in events:
                        if event['event'] == 'afe':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            for junc in event['SkipA2']:
                                row = [event['C1'].range_str(), 'A1', event['A1'].range_str(), 'C1_A1',
                                       junc.range_str()]
                                writer.writerow(src_common + row + self.quantifications(module, 's', junc))
                            for junc in event['SkipA1']:
                                row = [event['C1'].range_str(), 'A2', event['A2'].range_str(), 'C1_A2',
                                       junc.range_str()]
                                writer.writerow(trg_common + row + self.quantifications(module, 't', junc))
                        elif event['event'] == 'p_ale':
                            trg_common = self.common_data(module, 't')
                            row = ['N/A', 'A1', event['A1'].range_str(), 'C1_A1',
                                   'N/A']
                            writer.writerow(trg_common + row + self.quantifications(module, 't', event['A1']))


    def intron_retention(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Junction Name',
                                         'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'intron_retention.tsv')
        with open(os.path.join(self.config.directory, 'intron_retention.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex = self.as_types[module.idx]
                if not _complex or SHOW_COMPLEX_IN_ALL:
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

    def multi_exon_skipping(self):
        headers = self.common_headers + ['Reference Exon Coordinate', 'Exon Spliced With',
                                         'Exon Spliced With Coordinate', 'Tandem Exon Coordinates',
                                         'Junction Name', 'Junction Coordinate'] + self.quantification_headers
        self.start_headers(headers, 'multi_exon_skipping.tsv')
        with open(os.path.join(self.config.directory, 'multi_exon_skipping.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:
                events, _complex = self.as_types[module.idx]
                if not _complex or SHOW_COMPLEX_IN_ALL:
                    for event in events:
                        if event['event'] == 'multi_exon_skipping':
                            src_common = self.common_data(module, 's')
                            trg_common = self.common_data(module, 't')
                            row = [event['C1'].range_str(), 'C2', event['C2'].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), 'C1_C2',
                                   semicolon((x.range_str() for x in event['Skip']))]
                            writer.writerow(src_common + row + self.quantifications(module, 's'))
                            row = [event['C1'].range_str(), 'A1', event['As'][0].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), 'C1_A',
                                   semicolon((x.range_str() for x in event['Include1']))]
                            writer.writerow(src_common + row + self.quantifications(module, 's'))
                            row = [event['C2'].range_str(), 'C1', event['C1'].range_str(),
                                   semicolon((x.range_str() for x in event['As'])), 'C2_C1',
                                   semicolon((x.range_str() for x in event['Skip']))]
                            writer.writerow(trg_common + row + self.quantifications(module, 't'))
                            row = [event['C2'].range_str(), 'A<N>',
                                   semicolon((x.range_str() for x in event['As'])), 'A<N>_C2',
                                   semicolon((x.range_str() for x in event['Includes']))]
                            writer.writerow(trg_common + row + self.quantifications(module, 't'))

    def summary(self):
        """
        Write the summary style output file
        :param genes_modules: a list of (gene_id (str), gene_modules (obj)) tuples
        :return: NOTHING
        """
        headers = ['Module', 'LSV ID(s)', 'Cassette', 'Alt 3',
                             'Alt 5', 'Alt 3 and Alt 5', 'MXE', 'ALE',
                             'AFE', 'P_ALE', 'P_AFE', 'Multi Exon Skipping', 'Intron Retention', 'Complex']
        self.start_headers(headers, 'summary.tsv')

        with open(os.path.join(self.config.directory, 'summary.tsv'), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')


            for module in self.modules:
                events, _complex = self.as_types[module.idx]
                counts = OrderedDict()
                counts['cassette_exon'] = 0
                counts['alt3ss'] = 0
                counts['alt5ss'] = 0
                counts['alt3ss+alt5ss'] = 0
                counts['mutually_exclusive'] = 0
                counts['ale'] = 0
                counts['afe'] = 0
                counts['p_ale'] = 0
                counts['p_afe'] = 0
                counts['multi_exon_skipping'] = 0
                counts['intron_retention'] = 0
                for event in events:
                    if event['event'] in counts:
                        counts[event['event']] += 1



                writer.writerow(["%s_%d" % (self.gene_id, module.idx),
                                 semicolon(module.target_lsv_ids.union(module.source_lsv_ids))] +
                                [v if v else '' for v in counts.values()] + [str(_complex)]
                                )
