import csv
import os
from multiprocessing import Manager
from multiprocessing import Process, Pool
from multiprocessing.queues import JoinableQueue
from os.path import join

import h5py
import numpy as np

import vlsv
from voila import constants
from voila.constants import JUNCTION_TYPE_RNASEQ
from voila.hdf5 import HDF5
from voila.producer_consumer import ProducerConsumer
from voila.utils import utils_voila
from voila.utils.run_voila_utils import get_output_html
from voila.utils.voila_log import voila_log
from voila.vlsv import VoilaLsv

__author__ = 'abarrera'


class Voila(ProducerConsumer):
    VERSION = '/voila_file_version'

    def __init__(self, voila_file_name, mode):
        """
        Parse or edit the voila (quantifier output) file.
        :param voila_file_name: location of voila file
        :param mode: file mode passed to h5py
        """
        super(Voila, self).__init__()
        self.mode = mode
        self.file_name = voila_file_name
        self.hdf5 = None
        self.lsv_ids = None
        self.lsv_types = None
        self.gene_names = None
        self.file_version = None

    def __enter__(self):
        self.hdf5 = h5py.File(self.file_name, self.mode)

        if self.VERSION not in self.hdf5:
            if self.mode == constants.FILE_MODE.write:
                self.hdf5[self.VERSION] = constants.VOILA_FILE_VERSION

        try:
            self.file_version = self.hdf5[self.VERSION].value
        except KeyError:
            pass

        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def _producer(self):
        with h5py.File(self.file_name, 'r') as h:
            for lsv_id in h['lsvs']:
                self.queue.put(lsv_id)

    def _worker(self):
        with h5py.File(self.file_name, 'r') as h:
            while True:
                lsv_id = self.queue.get()
                lsv = VoilaLsv.easy_from_hdf5(h['lsvs'][lsv_id])

                if not self.lsv_types or lsv.lsv_type in self.lsv_types:
                    if not self.gene_names or lsv.name in self.gene_names:
                        if not self.lsv_ids or lsv.lsv_id in self.lsv_ids:
                            self.dict(lsv_id, lsv)

                self.queue.task_done()

    def close(self):
        try:
            self.hdf5.close()
        except ValueError:
            pass

    def add_lsv(self, voilaLsv):
        """
        Add VoilaLsv to Voila file.
        :param voilaLsv: VoilaLsv object
        :return: None
        """
        voilaLsv.to_hdf5(self.hdf5)

    def add_metainfo(self, genome, group1, experiments1, group2=None, experiments2=None):
        """
        Add metainfo to Voila file.
        :param genome: genome where the genes are found
        :param group1: first group name (used in psi and deltapsi)
        :param experiments1: first list of experiment names (used in psi and deltapsi)
        :param group2: second group (only deltapsi)
        :param experiments2: second list of experiment names (only deltapsi)
        :return: None
        """

        h = self.hdf5.create_group('/metainfo')

        h.attrs['genome'] = genome

        h.attrs['group1'] = group1
        experiments1_grp = h.create_group('experiments1')
        for index, item in enumerate(experiments1):
            experiments1_grp.attrs[str(index)] = item

        if group2 and experiments2:
            h.attrs['group2'] = group2
            experiments2_grp = h.create_group('experiments2')
            for index, item in enumerate(experiments2):
                experiments2_grp.attrs[str(index)] = item

    def get_metainfo(self):
        """
        Get metainfo from voila file.
        :return: dict
        """
        voila_log().info('Getting Voila Metainfo from {0} ...'.format(self.file_name))
        metainfo = {}
        h = self.hdf5['/metainfo']

        for key in h.attrs:
            metainfo[key] = h.attrs[key]

        for key in h:
            metainfo[key] = [h[key].attrs[attr] for attr in h[key].attrs]

        return metainfo

    def get_voila_lsv(self, lsv_id):
        """
        Get LSV by LSV id.
        :param lsv_id:
        :return: VoilaLsv
        """
        return VoilaLsv.easy_from_hdf5(self.hdf5['lsvs'][lsv_id])

    def get_voila_lsvs(self, lsv_types=None, lsv_ids=None, gene_names=None):
        """
        Get list of LSVs from voila file.
        :param lsv_types: search for lsv types
        :param lsv_ids: search for lsv ids
        :param gene_names: search for gene names
        :return: list
        """
        voila_log().info('Getting Voila LSVs from {0} ...'.format(self.file_name))

        self.lsv_types = lsv_types
        self.lsv_ids = lsv_ids
        self.gene_names = gene_names

        self.run()

        voila_lsvs = self.get_values()
        self.manager_shutdown()

        return voila_lsvs

    def check_version(self):
        if self.file_version != constants.VOILA_FILE_VERSION:
            voila_log().warning('Voila file version isn\'t current.  This will probably cause significant '
                                'issues with the voila output.  It would be best to run quantifier again with the '
                                'current version of MAJIQ.')


class VoilaInput(HDF5):
    """Standard input interface by experiment used by Voila"""

    def __init__(self, lsvs=(), metainfo=None):
        super(VoilaInput, self).__init__()
        print 'VoilaInput has been deprecated.  Use Voila instead.'
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


def tab_output(args, majiq_output):
    def semicolon_join(value_list):
        return ';'.join(str(x) for x in value_list)

    output_dir = args.output
    output_html = get_output_html(args, args.voila_file)
    log = voila_log()
    type_summary = args.type_analysis
    group1 = None
    group2 = None
    tsv_file = join(output_dir, output_html.rsplit('.html', 1)[0] + '.tsv')

    log.info("Creating Tab-delimited output file in %s..." % tsv_file)

    fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction']

    if 'delta' in type_summary:
        group1 = majiq_output['meta_exps']['group1']
        group2 = majiq_output['meta_exps']['group2']
        fieldnames = fieldnames[:3] + ['E(dPSI) per LSV junction',
                                       'P(|E(dPSI)|>=%.2f) per LSV junction' % args.threshold,
                                       '%s E(PSI)' % group1, '%s E(PSI)' % group2]

    fieldnames += ['LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                   'strand', 'Junctions coords', 'Exons coords', 'Exons Alternative Start', 'Exons Alternative End',
                   'IR coords']

    if 'voila_links' in majiq_output:
        fieldnames.append('Voila link')

    with open(tsv_file, 'w') as tsv:
        writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for gene in majiq_output['genes_dict']:
            for lsv in majiq_output['genes_dict'][gene]:

                if constants.ANALYSIS_DELTAPSI == type_summary:
                    lsv = lsv['lsv']

                row = {
                    '#Gene Name': lsv.name,
                    'Gene ID': gene,
                    'LSV ID': lsv.lsv_id,
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

                if constants.ANALYSIS_DELTAPSI == type_summary:
                    row.update({
                        'E(dPSI) per LSV junction': semicolon_join(
                            lsv.excl_incl[i][1] - lsv.excl_incl[i][0] for i in range(len(lsv.bins))
                        ),
                        'P(|dPSI|>=%.2f) per LSV junction' % args.threshold: semicolon_join(
                            vlsv.matrix_area(np.array(bin), args.threshold, collapsed_mat=True).sum() for bin in
                            lsv.bins
                        ),
                        '%s E(PSI)' % group1: semicolon_join(
                            '%.3f' % i for i in lsv.means_psi1
                        ),
                        '%s E(PSI)' % group2: semicolon_join(
                            '%.3f' % i for i in lsv.means_psi2
                        )
                    })

                if constants.ANALYSIS_PSI == type_summary:
                    row.update({
                        'E(PSI) per LSV junction': semicolon_join(lsv.means),
                        'Var(E(PSI)) per LSV junction': semicolon_join(lsv.variances)
                    })

                if 'voila_links' in majiq_output:
                    summary_path = majiq_output['voila_links'][lsv.name]
                    if not os.path.isabs(summary_path):
                        summary_path = join(os.getcwd(), output_dir, summary_path)
                    row['Voila link'] = constants.URL_COMPOSITE % (summary_path, lsv.name)

                writer.writerow(row)

    log.info("Delimited output file successfully created in: %s" % tsv_file)


def generic_feature_format_txt_files(args, majiq_output, out_gff3=False):
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
        log.info("Saving LSVs files in gff format ...")
    else:
        log.info("Saving LSVs files in gtf format ...")

    if 'genes_dict' not in majiq_output or not len(majiq_output['genes_dict']):
        log.warning("No gene information provided. Genes files are needed to calculate the gtf/gff files.")
        return

    header = "##gff-version 3"

    odir = join(output_dir, "static/doc/lsvs")
    utils_voila.create_if_not_exists(odir)

    for gkey, gvalue in majiq_output['genes_dict'].iteritems():
        for lsv_dict in gvalue:
            lsv = lsv_dict
            if type(lsv_dict) == dict:
                lsv = lsv_dict['lsv']
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

            except UnboundLocalError, e:
                log.warning("problem generating GTF file for %s" % lsv.id)
                log.error(e.message)
