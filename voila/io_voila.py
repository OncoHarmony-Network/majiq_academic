import csv
import os
from collections import defaultdict
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
    def __init__(self, voila_file_name, mode):
        super(Voila, self).__init__()
        self.mode = mode
        self.file_name = voila_file_name
        self.hdf5 = None
        self.lsv_ids = None
        self.lsv_types = None
        self.gene_names = None

    def __enter__(self):
        self.hdf5 = h5py.File(self.file_name, self.mode)
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def _producer(self):
        with h5py.File(self.file_name, 'r') as h:
            for id in h['lsvs']:
                self.queue.put(id)

    def _worker(self):
        with h5py.File(self.file_name, 'r') as h:
            while True:
                id = self.queue.get()
                lsv = VoilaLsv.easy_from_hdf5(h['lsvs'][id])

                if not self.lsv_types or lsv.lsv_type in self.lsv_types:
                    if not self.gene_names or lsv.name in self.gene_names:
                        if not self.lsv_ids or lsv.lsv_id in self.lsv_ids:
                            self.dict(id, lsv)

                self.queue.task_done()

    def close(self):
        try:
            self.hdf5.close()
        except ValueError:
            pass

    def add_lsv(self, voilaLsv):
        voilaLsv.to_hdf5(self.hdf5)

    def add_metainfo(self, genome, group1, experiments1, group2=None, experiments2=None):
        metainfo = {'group1': group1, 'experiments1': experiments1, 'genome': genome}

        if group2 and experiments2:
            metainfo.update({
                'experiments2': experiments2,
                'group2': group2
            })

        h = self.hdf5.create_group('/metainfo')
        metainfo = metainfo.copy()

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

    def get_metainfo(self):
        voila_log().info('Getting Voila Metainfo from {0} ...'.format(self.file_name))
        metainfo = {}
        h = self.hdf5['/metainfo']

        for key in h.attrs:
            metainfo[key] = h.attrs[key]

        for key in h:
            metainfo[key] = [h[key].attrs[attr] for attr in h[key].attrs]

        return metainfo

    def get_voila_lsv(self, lsv_id):
        return VoilaLsv.easy_from_hdf5(self.hdf5['lsvs'][lsv_id])

    def get_voila_lsvs(self, lsv_types=None, lsv_ids=None, gene_names=None):
        voila_log().info('Getting Voila LSVs from {0} ...'.format(self.file_name))
        self.lsv_types = lsv_types
        self.lsv_ids = lsv_ids
        self.gene_names = gene_names
        self.run()
        voila_lsvs = self.get_values()
        self.manager_shutdown()
        return voila_lsvs


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


def write_tab_output(args, majiq_output):
    """
    Create tab-delimited output file summarizing all the LSVs detected and quantified with MAJIQ.

    :param args: parsed input data
    :return:
    """

    if args.type_analysis == constants.COND_TABLE:
        cond_table_tab_output(args)
    else:
        tab_output(args, majiq_output)


def cond_table_tab_output(input_parsed):
    majiq_output = input_parsed.majiq_output
    lsvs = majiq_output['lsvs']
    sample_names = majiq_output['sample_names']
    log = voila_log()

    log.info('Creating cond-table TSV...')

    tsv_file = os.path.join(input_parsed.output_dir, input_parsed.output_html.split('.html')[0] + '.tsv')
    log.info(tsv_file)

    with open(tsv_file, 'w') as csvfile:
        fieldnames = ['Gene', 'LSV ID', '#Disagreeing', '#Changing samples', 'Junction'] + sample_names
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()

        for lsv in lsvs:
            row = {
                'Gene': lsvs[lsv]['gene'],
                'LSV ID': lsv,
                '#Disagreeing': lsvs[lsv]['ndisagree'],
                '#Changing samples': lsvs[lsv]['nchangs'],
                'Junction': lsvs[lsv]['njunc']
            }

            for index, sample_name in enumerate(sample_names):
                sample = round(lsvs[lsv]['expecs'][index], 3)
                if sample > -1:
                    row[sample_name] = sample

            writer.writerow(row)


def tab_output(args, majiq_output):
    def semicolon_join(value_list):
        return ';'.join(str(x) for x in value_list)

    output_dir = args.output
    output_html = get_output_html(args, args.majiq_quantifier)
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
                    'De Novo Junctions': any(
                        junc.junction_type == JUNCTION_TYPE_RNASEQ for junc in lsv.junctions
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
                        'P(|E(dPSI)|>=%.2f) per LSV junction' % args.threshold: semicolon_join(
                            vlsv.matrix_area(np.array(bin), args.threshold, collapsed_mat=True).sum() for bin in
                            lsv.bins
                        ),
                        '%s E(PSI)' % group1: semicolon_join(
                            '%.3f' % vlsv.get_expected_psi(np.array(lsv.psi1[i])) for i in range(len(lsv.bins))
                        ),
                        '%s E(PSI)' % group2: semicolon_join(
                            '%.3f' % vlsv.get_expected_psi(np.array(lsv.psi2[i])) for i in range(len(lsv.bins))
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


# def load_dpsi_tab(tab_files_list, sample_names, thres_change=None, filter_genes=None, filter_lsvs=None,
#                   pairwise_dir=None, outdir=None):

def load_dpsi_tab(args):
    pairwise_dir = args.pairwise_dir
    outdir = args.output
    tab_files_list = args.sample_files
    filter_genes = args.gene_names
    filter_lsvs = args.lsv_ids
    sample_names = args.sample_names
    thres_change = args.threshold_change

    """Load LSV delta psi information from tab-delimited file."""
    lsvs_dict = defaultdict(lambda: defaultdict(lambda: None))
    if pairwise_dir is None:
        pairwise_dir = os.getcwd()

    root_path = None
    path_prefix = '/'.join(['..'] * len(outdir.strip('./').split('/'))) + '/'
    if len(outdir.strip('./')) == 0:
        path_prefix = './'
    # 2-step process:
    #   1. create a data structure finding the most changing junction,
    #   2. select the expected psi from the most changing junction more frequent in the set
    for idx, tab_file in enumerate(tab_files_list):
        with open(tab_file, 'r') as tabf:
            for line in tabf:
                if line.startswith("#"):
                    continue

                fields = line.split()

                if root_path is None:
                    pr, linkk = os.path.split(fields[-1])
                    linkk = linkk.split('#')[0]
                    while not os.path.exists(pairwise_dir + '/' + linkk) and len(pr) > 0:
                        pr, aux = os.path.split(pr)
                        linkk = aux + '/' + linkk

                    if len(pr) == 0:
                        raise Exception('Couldn\'t determine links to delta psi summaries')
                    root_path = pr

                if filter_genes:
                    if fields[0] not in filter_genes and fields[1] not in filter_genes:
                        continue

                if filter_lsvs:
                    if fields[2].upper() not in filter_lsvs:
                        continue

                expecs = [float(aa) for aa in fields[3].split(";")]

                if lsvs_dict[fields[2]]['expecs'] is None:
                    lsvs_dict[fields[2]]['expecs'] = [[]] * len(sample_names)
                    lsvs_dict[fields[2]]['expecs_marks'] = [None] * len(sample_names)
                    lsvs_dict[fields[2]]['links'] = [None] * len(sample_names)
                    lsvs_dict[fields[2]]['njunc'] = [-1] * len(sample_names)

                idx_max = np.argmax([abs(ee) for ee in expecs])

                lsvs_dict[fields[2]]['expecs'][idx] = expecs
                lsvs_dict[fields[2]]['njunc'][idx] = idx_max
                lsvs_dict[fields[2]]['links'][idx] = path_prefix + pairwise_dir + fields[-1].split(root_path)[1]
                lsvs_dict[fields[2]]['gene'] = fields[0]

    for lsv_idx in lsvs_dict.keys():
        if np.max([abs(bb) for ff in lsvs_dict[lsv_idx]['expecs'] for bb in ff]) < thres_change:
            del lsvs_dict[lsv_idx]  # Remove LSVs not passing the changing threshold
            continue

        idx_most_freq = np.argmax(np.bincount(np.array(lsvs_dict[lsv_idx]['njunc'])[
                                                  (np.array(lsvs_dict[lsv_idx]['njunc']) > -1) & np.array(
                                                      [np.any(np.array([abs(fff) for fff in expec]) > thres_change) for
                                                       expec in lsvs_dict[lsv_idx]['expecs']])]))

        # Mark adjusted most changing junction
        lsvs_dict[lsv_idx]['expecs_marks'] = ~np.array(idx_most_freq == lsvs_dict[lsv_idx]['njunc'])

        for idx_exp, expec in enumerate(lsvs_dict[lsv_idx]['expecs']):
            if len(expec) > 0:
                lsvs_dict[lsv_idx]['expecs'][idx_exp] = expec[idx_most_freq]
            else:
                lsvs_dict[lsv_idx]['expecs'][idx_exp] = -1

        lsvs_dict[lsv_idx]['nchangs'] = np.count_nonzero(
            [abs(ee) > thres_change for ee in lsvs_dict[lsv_idx]['expecs'] if ee > -1])

        lsvs_dict[lsv_idx]['njunc'] = idx_most_freq

        exist_expecs = np.array(lsvs_dict[lsv_idx]['expecs'])[(np.array(lsvs_dict[lsv_idx]['expecs']) > -1) & (
            np.array([abs(xx) for xx in lsvs_dict[lsv_idx]['expecs']]) > thres_change)]

        lsvs_dict[lsv_idx]['ndisagree'] = len(exist_expecs) - max(
            (np.count_nonzero(exist_expecs > 0), np.count_nonzero(exist_expecs <= 0)))

        lsvs_dict[lsv_idx]['nagree'] = len(exist_expecs) - min(
            (np.count_nonzero(exist_expecs > 0), np.count_nonzero(exist_expecs <= 0)))

    return lsvs_dict


def generic_feature_format_txt_files(input_parsed, out_gff3=False):
    """
    Create GFF3 files for each LSV.
    :param input_parsed: parsed input data
    :param out_gff3: output as a GFF3 file
    :return: None
    """

    log = voila_log()
    majiq_output = input_parsed.majiq_output
    output_dir = input_parsed.output_dir

    if input_parsed.type_summary == constants.COND_TABLE:
        log.info('Skipping generating gff3 format.')
        return

    log.info("Saving LSVs files in gff3 format ...")
    if 'genes_dict' not in majiq_output or len(majiq_output['genes_dict']) < 1:
        log.warning("No gene information provided. Genes files are needed to calculate the gff3 files.")
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
                lsv_gff3_str = lsv.get_gff3()
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

    log.info("GTF files for LSVs saved in %s" % odir)
