import csv
import itertools
import os
from multiprocessing import Manager
from multiprocessing import Process, Pool
from multiprocessing.queues import JoinableQueue
from os.path import join

import h5py
import numpy as np

from voila import constants, vlsv
from voila.api import Voila, SpliceGraph
from voila.api.view_matrix import ViewDeltaPsi, ViewPsi
from voila.constants import JUNCTION_TYPE_RNASEQ
from voila.hdf5 import HDF5
from voila.utils import utils_voila
from voila.utils.run_voila_utils import get_output_html
from voila.utils.voila_log import voila_log
from voila.vlsv import VoilaLsv

__author__ = 'abarrera'


class VoilaInput(HDF5):
    """Standard input interface by experiment used by Voila"""

    def __init__(self, lsvs=(), metainfo=None):
        super(VoilaInput, self).__init__()
        print('VoilaInput has been deprecated.  Use Voila instead.')
        self.lsvs = lsvs
        self.metainfo = metainfo

    def get_lsvs(self):
        return self.lsvs

    def add_lsv(self, l):
        self.lsvs.append(l)

    def samples_metainfo(self):
        """Sample information generator, one by one."""
        for sample_info in self.metainfo:
            yield sample_info

    def encode_metainfo(self, h):
        h = h.create_group('/metainfo')
        metainfo = self.metainfo.copy()

        experiments1 = h.create_group('experiments1')
        for index, item in enumerate(metainfo['experiments1']):
            experiments1.attrs[str(index)] = item

        del metainfo['experiments1']

        if 'group2' in metainfo and 'experiments2' in metainfo:
            experiments2 = h.create_group('experiments2')
            for index, item in enumerate(metainfo['experiments2']):
                experiments2.attrs[str(index)] = item

            del metainfo['experiments2']

        for key in metainfo:
            h.attrs[key] = metainfo[key]

    def decode_metainfo(self, h):
        self.metainfo = {}
        for key in h.attrs:
            self.metainfo[key] = h.attrs[key]

        for key in h:
            self.metainfo[key] = [h[key].attrs[attr] for attr in h[key].attrs]

    def exclude(self):
        return ['metainfo', 'lsvs']

    def to_hdf5(self, h, use_id=True):
        # metainfo
        self.encode_metainfo(h)

        # lsvs
        for lsv in self.lsvs:
            lsv.to_hdf5(h, use_id)

        super(VoilaInput, self).to_hdf5(h, use_id)

    def from_hdf5(self, h):
        # metainfo
        self.decode_metainfo(h['metainfo'])

        # lsvs
        self.lsvs = [VoilaLsv.easy_from_hdf5(h['lsvs'][lsv_id]) for lsv_id in h['lsvs']]

        return super(VoilaInput, self).from_hdf5(h)

    @classmethod
    def metainfo_to_hdf5(cls, h, genome, group1, experiments1, group2=None, experiments2=None):
        metainfo = {'group1': group1, 'experiments1': experiments1, 'genome': genome}
        if group2 and experiments2:
            metainfo['experiments2'] = experiments2
            metainfo['group2'] = group2
        vi = cls()
        vi.metainfo = metainfo
        vi.encode_metainfo(h['/'])

    @classmethod
    def from_hdf5_file(cls, hdf5_filename, lsv_types=None, lsv_ids=None, gene_names=None):
        """
        Create VoilaInput object from HDF5 file.  This will process each of the VoilaLsvs in their own thread using the
        Producer Consumer design pattern.
        :param hdf5_filename: HDF5 filename string
        :return: VoilaInput object
        """

        def worker():
            with h5py.File(hdf5_filename, 'r') as h:
                while True:
                    id = queue.get()
                    lsv = VoilaLsv.easy_from_hdf5(h['lsvs'][id])

                    if not lsv_types or lsv.lsv_type in lsv_types:
                        if not gene_names or lsv.name in gene_names:
                            if not lsv_ids or lsv.lsv_id in lsv_ids:
                                manage_dict[id] = lsv
                    queue.task_done()

        def producer():
            with h5py.File(hdf5_filename, 'r') as h:
                for id in h['lsvs']:
                    queue.put(id)

        log = voila_log()
        if not os.path.isfile(hdf5_filename):
            log.error('unable to load file: {0}'.format(hdf5_filename))
            raise IOError('Voila input file does not exist.')

        log.info('Loading {0}.'.format(hdf5_filename))

        voila_input = VoilaInput()

        queue = JoinableQueue()

        manage_dict = Manager().dict()

        producer_proc = Process(target=producer)
        producer_proc.daemon = True
        producer_proc.start()

        pool = Pool(constants.PROCESS_COUNT, worker)

        producer_proc.join()
        queue.join()

        pool.close()
        queue.close()

        with h5py.File(hdf5_filename, 'r') as h:
            voila_input.decode_metainfo(h['metainfo'])

        voila_input.lsvs = manage_dict.values()

        return voila_input


def load_dpairs(pairwise_dir, majiq_output):
    """
    Load pairwise files from MAJIQ analysis.

    :param str pairwise_dir: directory containing pairwise comparisons produced by MAJIQ.
    :param majiq_output: parsed data from old_majiq.
    :return: list of deltapsi lsvs
    :return: name of condition 1
    :return: name of condition 2
    """
    meta_exps = majiq_output['meta_exps']
    lmajiq_pairs = [[None for i in range(len(meta_exps[1]))] for j in range(len(meta_exps[0]))]

    lsv_names = majiq_output['genes_dict'].keys()

    group1_name = meta_exps[0][0]['group']
    group2_name = meta_exps[1][0]['group']

    for idx1 in range(len(meta_exps[0])):
        for idx2 in range(len(meta_exps[1])):
            pairwise_file = "%s/%s_%d_%s_%d.deltapsi.pickle" % (
                pairwise_dir, group1_name, idx1 + 1, group2_name, idx2 + 1)
            try:
                lmajiq_pairs[idx1][idx2] = utils_voila.get_lsv_delta_exp_data(
                    pairwise_file,
                    show_all=True,
                    gene_name_list=lsv_names
                )
            except IOError:
                pass
    return lmajiq_pairs, group1_name, group2_name


def filter_exons(exons):
    for exon in exons:
        if exon.start == -1:
            yield 'nan', exon.end
        elif exon.end == -1:
            yield exon.start, 'nan'
        else:
            yield exon.start, exon.end


def delta_psi_tab_output(args, voila_links):
    def semicolon_join(value_list):
        return ';'.join(str(x) for x in value_list)

    log = voila_log()
    log.info("Creating Tab-delimited output file")
    output_html = get_output_html(args, args.voila_file)
    tsv_file = join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

    with ViewDeltaPsi(args.voila_file, 'r') as m, SpliceGraph(args.splice_graph) as sg:
        metadata = m.metadata
        experiment = m.experiment_names[0][0]

        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction']
        group1 = metadata['group_names'][0]
        group2 = metadata['group_names'][1]
        fieldnames = fieldnames[:3] + ['E(dPSI) per LSV junction',
                                       'P(|dPSI|>=%.2f) per LSV junction' % args.threshold,
                                       '%s E(PSI)' % group1, '%s E(PSI)' % group2]

        fieldnames += ['LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                       'strand', 'Junctions coords', 'Exons coords', 'IR coords']

        if voila_links:
            fieldnames.append('Voila link')

        with open(tsv_file, 'w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for gene_id in m.get_gene_ids(args):
                gene = sg.gene(gene_id).get

                for lsv_id in m.get_lsv_ids(args, gene_id):

                    lsv = m.delta_psi(lsv_id)
                    lsv_junctions = tuple(gene.lsv_junctions(lsv_id))
                    lsv_exons = tuple(gene.lsv_exons(lsv_id, lsv_junctions))
                    group_means = tuple(lsv.group_means)
                    excl_incl = tuple(lsv.excl_incl)

                    row = {
                        '#Gene Name': gene.name,
                        'Gene ID': gene_id,
                        'LSV ID': lsv_id,
                        'LSV Type': lsv.lsv_type,
                        'A5SS': lsv.prime5,
                        'A3SS': lsv.prime3,
                        'ES': lsv.exon_skipping,
                        'Num. Junctions': lsv.junction_count,
                        'Num. Exons': lsv.exon_count,
                        'chr': gene.chromosome,
                        'strand': gene.strand,
                        'De Novo Junctions': semicolon_join(
                            int(next(junc.get_junction_types([experiment])) == JUNCTION_TYPE_RNASEQ) for junc in
                            lsv_junctions
                        ),
                        'Junctions coords': semicolon_join(
                            '{0}-{1}'.format(junc.start, junc.end) for junc in lsv_junctions
                        ),
                        'Exons coords': semicolon_join(
                            '{0}-{1}'.format(start, end) for start, end in filter_exons(lsv_exons)
                        ),
                        'IR coords': semicolon_join(
                            '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                        ),
                        'E(dPSI) per LSV junction': semicolon_join(
                            excl_incl[i][1] - excl_incl[i][0] for i in
                            range(np.size(lsv.bins, 0))
                        ),
                        'P(|dPSI|>=%.2f) per LSV junction' % args.threshold: semicolon_join(
                            vlsv.matrix_area(np.array(bin), args.threshold, collapsed_mat=True).sum() for bin in
                            lsv.bins
                        ),
                        '%s E(PSI)' % group1: semicolon_join(
                            '%.3f' % i for i in group_means[0]
                        ),
                        '%s E(PSI)' % group2: semicolon_join(
                            '%.3f' % i for i in group_means[1]
                        )
                    }

                    if voila_links:
                        summary_path = voila_links[gene_id]
                        if not os.path.isabs(summary_path):
                            summary_path = join(os.getcwd(), args.output, summary_path)
                        row['Voila link'] = "file://{0}".format(summary_path)

                    log.debug('Write TSV row for {0}'.format(lsv_id))

                    writer.writerow(row)

    log.info("Delimited output file successfully created in: %s" % tsv_file)


def psi_tab_output(args, voila_links):
    def semicolon_join(value_list):
        return ';'.join(str(x) for x in value_list)

    log = voila_log()
    log.info("Creating Tab-delimited output file")

    output_html = get_output_html(args, args.voila_file)
    tsv_file = join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

    with ViewPsi(args.voila_file, 'r') as m, SpliceGraph(args.splice_graph) as sg:
        experiment = m.experiment_names[0][0]
        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction',
                      'LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords']
        if voila_links:
            fieldnames.append('Voila link')

        with open(tsv_file, 'w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for gene_id in m.get_gene_ids(args):
                gene = sg.gene(gene_id).get

                for lsv_id in m.get_lsv_ids(args, gene_id):
                    lsv = m.psi(lsv_id)
                    lsv_junctions = tuple(gene.lsv_junctions(lsv_id))
                    lsv_exons = tuple(gene.lsv_exons(lsv_id, lsv_junctions))

                    row = {
                        '#Gene Name': gene.name,
                        'Gene ID': gene_id,
                        'LSV ID': lsv_id,
                        'LSV Type': lsv.lsv_type,
                        'A5SS': lsv.prime5,
                        'A3SS': lsv.prime3,
                        'ES': lsv.exon_skipping,
                        'Num. Junctions': lsv.junction_count,
                        'Num. Exons': lsv.exon_count,
                        'chr': gene.chromosome,
                        'strand': gene.strand,
                        'De Novo Junctions': semicolon_join(
                            int(next(junc.get_junction_types([experiment])) == JUNCTION_TYPE_RNASEQ) for junc in
                            lsv_junctions
                        ),
                        'Junctions coords': semicolon_join(
                            '{0}-{1}'.format(junc.start, junc.end) for junc in lsv_junctions
                        ),
                        'Exons coords': semicolon_join(
                            '{0}-{1}'.format(start, end) for start, end in filter_exons(lsv_exons)
                        ),
                        'IR coords': semicolon_join(
                            '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                        ),
                        'E(PSI) per LSV junction': semicolon_join(lsv.means),
                        'Var(E(PSI)) per LSV junction': semicolon_join(lsv.variances)
                    }

                    if voila_links:
                        summary_path = voila_links[gene_id]
                        if not os.path.isabs(summary_path):
                            summary_path = join(os.getcwd(), args.output, summary_path)
                        row['Voila link'] = "file://{0}".format(summary_path)

                    log.debug('Write TSV row for {0}'.format(lsv_id))

                    writer.writerow(row)

    log.info("Delimited output file successfully created in: %s" % tsv_file)


def generic_feature_format_txt_files(args, out_gff3=False):
    """
    Create GFF3 files for each LSV.
    :param majiq_output: majiq data
    :param args: parsed input data
    :param out_gff3: output as a GFF3 file
    :return: None
    """

    log = voila_log()
    output_dir = args.output

    if out_gff3:
        log.info("Create GFF files for LSVs")
    else:
        log.info("Create GTF files for LSVs")

    header = "##gff-version 3"

    odir = join(output_dir, "static/doc/lsvs")
    utils_voila.create_if_not_exists(odir)

    with Voila(args.voila_file, 'r') as v:
        for lsv in v.get_voila_lsvs(args):
            lsv_file_basename = "%s/%s" % (odir, lsv.lsv_id)
            try:
                lsv_gff3_str = lsv.to_gff3()
                utils_voila.gff2gtf(lsv_gff3_str.split('\n'), "%s.gtf" % lsv_file_basename)

                # not accessible from command line
                if out_gff3:
                    gff_file = "%s.gff3" % (lsv_file_basename)
                    with open(gff_file, 'w') as ofile:
                        ofile.write(header + "\n")
                        ofile.write(lsv_gff3_str + "\n")

            except UnboundLocalError as e:
                log.warning("problem generating GTF file for %s" % lsv.id)
                log.error(e)


def get_lsv_info_fieldnames():
    return ['Gene Name', 'Gene ID', 'LSV ID']


def get_lsv_info(lsv):
    return {
        'Gene Name': lsv.name,
        'Gene ID': lsv.lsv_id.split(':')[0],
        'LSV ID': lsv.lsv_id,
    }


def get_lsv_extra_info_fieldnames():
    return ['LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr', 'strand',
            'Junctions coords', 'Exons coords', 'Exons Alternative Start', 'Exons Alternative End', 'IR coords']


def get_lsv_extra_info(lsv):
    return {
        'LSV Type': lsv.lsv_type,
        'A5SS': lsv.categories['prime5'],
        'A3SS': lsv.categories['prime3'],
        'ES': lsv.categories['ES'],
        'Num. Junctions': lsv.categories['njuncs'],
        'Num. Exons': lsv.categories['nexons'],
        'chr': lsv.chromosome,
        'strand': lsv.strand,
        'De Novo Junctions': semicolon_join(
            int(junc.junction_type == JUNCTION_TYPE_RNASEQ) for junc in lsv.junctions
        ),
        'Junctions coords': semicolon_join(
            '{0}-{1}'.format(junc.start, junc.end) for junc in lsv.junctions
        ),
        'Exons coords': semicolon_join(
            '{0}-{1}'.format(e.start, e.end) for e in lsv.exons
        ),
        'Exons Alternative Start': semicolon_join(
            '|'.join(str(a) for a in e.alt_starts) for e in lsv.exons if e.alt_starts
        ),
        'Exons Alternative End': semicolon_join(
            '|'.join(str(a) for a in e.alt_ends) for e in lsv.exons if e.alt_ends
        ),
        'IR coords': semicolon_join(
            '{0}-{1}'.format(e.start, e.end) for e in lsv.exons if e.intron_retention
        )
    }


def semicolon_join(value_list):
    return ';'.join(str(x) for x in value_list)


def het_tab_output(args):
    voila_log().info('Creating HET TSV file...')

    output_html = get_output_html(args, args.voila_file)
    tsv_file = join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

    with Voila(args.voila_file, 'r') as v:

        metainfo = v.get_metainfo()

        fieldnames = get_lsv_info_fieldnames() + list(metainfo['stat_names']) + get_lsv_extra_info_fieldnames()

        lsv_fieldnames = ['Junction ID'] + list(itertools.chain.from_iterable(metainfo['experiment_names']))

        with open(tsv_file, 'w') as tsvfile:
            tsv = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
            tsv.writeheader()

            for lsv in v.get_voila_lsvs(args):
                row = get_lsv_info(lsv)
                row.update(get_lsv_extra_info(lsv))

                for stat_name, junction_stat in zip(metainfo['stat_names'], lsv.het.junction_stats):
                    row[stat_name] = semicolon_join(junction_stat)

                tsv.writerow(row)

                rows = {}

        for lsv in v.get_voila_lsvs(args):

            for group, experiment_names in zip(lsv.het.groups, metainfo['experiment_names']):
                for experiment_index, experiment_name in enumerate(experiment_names):
                    for junction_index, junction_id in enumerate(lsv.junction_ids()):

                        psi = group.get_psi(experiment_index=experiment_index, junction_index=junction_index)

                        try:
                            rows[junction_id][experiment_name] = psi
                        except KeyError:
                            rows[junction_id] = {experiment_name: psi}

            with open(join(args.output, lsv.lsv_id.replace(':', '_') + '.tsv'), 'w') as tsvfile:
                tsv = csv.DictWriter(tsvfile, fieldnames=lsv_fieldnames, delimiter='\t')
                tsv.writeheader()
                for junction_id, row_dict in rows.items():
                    row = {'Junction ID': junction_id}
                    row.update({column: value for column, value in row_dict.items()})
                    tsv.writerow(row)
