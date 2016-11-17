import csv
import os
from collections import defaultdict
from multiprocessing import Manager, Process, Pool
from multiprocessing.queues import JoinableQueue

import h5py
import numpy as np

import vlsv
from voila import constants
from voila.hdf5 import HDF5
from voila.utils import utils_voila
from voila.utils.voilaLog import voilaLog
from voila.vlsv import VoilaLsv

__author__ = 'abarrera'


class Voila(object):
    def __init__(self, voila_file_name, mode):
        self.mode = mode
        self.file_name = voila_file_name
        self.hdf5 = None

    def __enter__(self):
        self.hdf5 = h5py.File(self.file_name, self.mode)
        return self

    def __exit__(self, type, value, traceback):
        self.hdf5.close()

    def add_lsv(self, voilaLsv):
        voilaLsv.to_hdf5(self.hdf5)

    def add_metainfo(self, genome, group1, experiments1, group2=None, experiments2=None):
        metainfo = {'group1': group1, 'experiments1': experiments1, 'genome': genome}
        if group2 and experiments2:
            metainfo['experiments2'] = experiments2
            metainfo['group2'] = group2
        vi = VoilaInput()
        vi.metainfo = metainfo
        vi.encode_metainfo(self.hdf5['/'])


class VoilaInput(HDF5):
    """Standard input interface by experiment used by Voila"""

    def __init__(self, lsvs=(), metainfo=None):
        super(VoilaInput, self).__init__()
        self.lsvs = lsvs
        self.metainfo = metainfo

    def get_lsvs(self):
        return self.lsvs

    def add_lsv(self, l):
        self.lsvs.append(l)

    def get_metainfo(self):
        """Retrieve experiment metainfo (primarily samples information, expandable)"""
        return self.metainfo

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
    def from_hdf5_file(cls, hdf5_filename):
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
                    manage_dict[id] = VoilaLsv.easy_from_hdf5(h['lsvs'][id])
                    queue.task_done()

        def producer():
            with h5py.File(hdf5_filename, 'r') as h:
                for id in h['lsvs']:
                    queue.put(id)

        log = voilaLog()
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


def write_tab_output(input_parsed):
    """
    Create tab-delimited output file summarizing all the LSVs detected and quantified with MAJIQ.

    :param input_parsed: parsed input data
    :return:
    """

    if input_parsed.type_summary == constants.COND_TABLE:
        cond_table_tab_output(input_parsed)
    else:
        tab_output(input_parsed)


def cond_table_tab_output(input_parsed):
    majiq_output = input_parsed.majiq_output
    lsvs = majiq_output['lsvs']
    sample_names = majiq_output['sample_names']
    log = voilaLog()

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


def tab_output(input_parsed):
    output_dir = input_parsed.output_dir
    output_html = input_parsed.output_html
    log = voilaLog()
    pairwise_dir = input_parsed.pairwise_dir
    majiq_output = input_parsed.majiq_output
    threshold = input_parsed.threshold
    type_summary = input_parsed.type_summary

    ofile_str = "%s%s.%s" % (output_dir, output_html.rsplit('.html', 1)[0], constants.EXTENSION)
    tlb_categx = {'A5SS': 'prime5', 'A3SS': 'prime3', 'Num. Junctions': 'njuncs', 'Num. Exons': 'nexons', 'ES': 'ES'}

    log.info("Creating Tab-delimited output file in %s..." % ofile_str)

    if pairwise_dir:
        # In deltapsi, add columns with pairwise comparisons between group members
        log.info("Load pairwise comparison files from %s..." % pairwise_dir)
        lmajiq_pairs, group1_name, group2_name = load_dpairs(pairwise_dir, majiq_output)

    with open(ofile_str, 'w+') as ofile:
        headers = ['#Gene Name',
                   'Gene ID',
                   'LSV ID',
                   'E(PSI) per LSV junction',
                   'Var(E(PSI)) per LSV junction',
                   'LSV Type',
                   'A5SS',
                   'A3SS',
                   'ES',
                   'Num. Junctions',
                   'Num. Exons',
                   'De Novo Junctions?',
                   'chr',
                   'strand',
                   'Junctions coords',
                   'Exons coords',
                   'Exons Alternative Start',
                   'Exons Alternative End',
                   'IR coords']
        if 'voila_links' in majiq_output.keys():
            headers.append('Voila link')

        if 'delta' in type_summary:
            headers[3] = 'E(dPSI) per LSV junction'
            headers[4] = 'P(|E(dPSI)|>=%.2f) per LSV junction' % threshold
            psi_headers = ['%s E(PSI)' % majiq_output['meta_exps']['group1'],
                           '%s E(PSI)' % majiq_output['meta_exps']['group2']]
            headers = headers[:5] + psi_headers + headers[5:]

            if pairwise_dir:
                for idx1 in xrange(len(lmajiq_pairs)):
                    for idx2 in xrange(len(lmajiq_pairs[0])):
                        headers.append("%s_%d_%s_%d" % (group1_name, idx1 + 1, group2_name, idx2 + 1))

                exp_names_map = ['#Group names and file names mapping']
                for iexp in xrange(len(lmajiq_pairs)):
                    exp_names_map.append(
                        "#%s_%d=%s" % (group1_name, iexp + 1, lmajiq_pairs[0][0]['meta_exps'][0][iexp]['experiment']))
                for iexp in xrange(len(lmajiq_pairs[0])):
                    exp_names_map.append(
                        "#%s_%d=%s" % (group2_name, iexp + 1, lmajiq_pairs[0][0]['meta_exps'][1][iexp]['experiment']))
                ofile.write('\n'.join(exp_names_map))
                ofile.write('\n')
                ofile.write('#\n')

        # ofile.write("#Tab-delimited file\n#\n")
        ofile.write(constants.DELIMITER.join(headers))
        ofile.write('\n')

        for gene in majiq_output['genes_dict']:
            for llsv_dict in majiq_output['genes_dict'][gene]:
                llsv = llsv_dict
                if type(llsv_dict) == dict:
                    llsv = llsv_dict['lsv']
                lline = []
                lline.extend([llsv.lsv_graphic.name, gene, llsv.lsv_graphic.id])
                lexpected = []
                lconfidence = []
                lexpecs_psi1 = []
                lexpecs_psi2 = []
                for i, bins in enumerate(llsv.bins):
                    if 'delta' in type_summary:
                        lexpected.append(str(-llsv.excl_incl[i][0] + llsv.excl_incl[i][1]))
                        lconfidence.append(str(vlsv.matrix_area(np.array(bins), threshold, collapsed_mat=True).sum()))
                        lexpecs_psi1.append('%.3f' % vlsv.get_expected_psi(np.array(llsv.psi1[i])))
                        lexpecs_psi2.append('%.3f' % vlsv.get_expected_psi(np.array(llsv.psi2[i])))
                    else:
                        lexpected.append(repr(llsv.means[i]))
                        lconfidence.append(repr(llsv.variances[i]))

                lline.append(';'.join(lexpected))
                lline.append(';'.join(lconfidence))
                if 'delta' in type_summary:
                    lline.append(';'.join(lexpecs_psi1))
                    lline.append(';'.join(lexpecs_psi2))

                lline.append(llsv.lsv_graphic.type)
                lline.append(repr(llsv.categories[tlb_categx['A5SS']]))
                lline.append(repr(llsv.categories[tlb_categx['A3SS']]))
                lline.append(repr(llsv.categories[tlb_categx['ES']]))
                lline.append(repr(llsv.categories[tlb_categx['Num. Junctions']]))
                lline.append(repr(llsv.categories[tlb_categx['Num. Exons']]))
                lline.append(str(int(np.any([junc.type_junction == 1 for junc in llsv.lsv_graphic.junctions]))))

                lline.append(llsv.lsv_graphic.chrom)
                lline.append(llsv.lsv_graphic.strand)

                lline.append(';'.join(
                    ['-'.join(str(c) for c in junc.coords) for junc in llsv.lsv_graphic.junctions]))
                lline.append(
                    ';'.join(['-'.join(str(c) for c in exon.coords) for exon in llsv.lsv_graphic.exons]))

                try:
                    lline.append(';'.join(
                        ['|'.join([str(c) for c in exon.alt_starts]) for exon in llsv.lsv_graphic.exons]))
                    lline.append(';'.join(
                        ['|'.join([str(c) for c in exon.alt_ends]) for exon in llsv.lsv_graphic.exons]))
                except TypeError:
                    pass

                lline.append(
                    ';'.join([repr(exon.coords) for exon in llsv.lsv_graphic.exons if exon.intron_retention]))

                if pairwise_dir:
                    llpairwise = []
                    for idx1 in range(len(lmajiq_pairs)):
                        for idx2 in range(len(lmajiq_pairs[0])):
                            lpairwise = []
                            if gene in lmajiq_pairs[idx1][idx2]['genes_dict']:
                                for llsv_tmp in lmajiq_pairs[idx1][idx2]['genes_dict'][gene]:
                                    if llsv_tmp[0].id == llsv.id:
                                        lsv_pair = llsv_tmp[0]
                                        break
                                else:
                                    log.warning("LSV %s present in deltagroup but missing in %s." %
                                                (llsv.id, "%s_%d_%s_%d" % (group1_name, idx1 + 1,
                                                                           group2_name, idx2 + 1)))
                                    lpairwise.append('N/A')
                                    continue
                                for iway in range(len(llsv.bins)):
                                    lpairwise.append(str(sum(lsv_pair.excl_incl[iway])))
                            else:
                                lpairwise.append('N/A')
                            llpairwise.append(';'.join(lpairwise))
                    lline.extend(llpairwise)

                if 'voila_links' in majiq_output.keys():
                    summary_path = majiq_output['voila_links'][llsv.lsv_graphic.name]
                    if not os.path.isabs(summary_path):
                        summary_path = "%s/%s/%s" % (os.getcwd(), output_dir, summary_path)
                    lline.append(constants.URL_COMPOSITE % (summary_path, llsv.lsv_graphic.name))
                ofile.write(constants.DELIMITER.join(lline))
                ofile.write('\n')

    log.info("Delimited output file successfully created in: %s" % ofile_str)


def load_dpsi_tab(tab_files_list, sample_names, thres_change=None, filter_genes=None, filter_lsvs=None,
                  pairwise_dir=None, outdir=None):
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

        lsvs_dict[lsv_idx]['expecs_marks'] = ~np.array(
            idx_most_freq == lsvs_dict[lsv_idx]['njunc'])  # Mark adjusted most changing junction
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

    return lsvs_dict


def create_gff3_txt_files(input_parsed, out_gff3=False):
    """
    Create GFF3 files for each LSV.
    :param input_parsed: parsed input data
    :param out_gff3: output as a GFF3 file
    :return: None
    """

    log = voilaLog()
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

    odir = output_dir + "/static/doc/lsvs"
    utils_voila.create_if_not_exists(odir)
    for gkey, gvalue in majiq_output['genes_dict'].iteritems():
        for lsv_dict in gvalue:
            lsv = lsv_dict
            if type(lsv_dict) == dict:
                lsv = lsv_dict['lsv']
            lsv_file_basename = "%s/%s" % (odir, lsv.lsv_graphic.id)

            try:
                lsv_gff3_str = lsv.get_gff3()
                utils_voila.gff2gtf(lsv_gff3_str.split('\n'), "%s.gtf" % lsv_file_basename)
                if out_gff3:
                    gff_file = "%s.gff3" % (lsv_file_basename)
                    with open(gff_file, 'w') as ofile:
                        ofile.write(header + "\n")
                        ofile.write(lsv_gff3_str + "\n")
            except UnboundLocalError, e:
                log.warning("problem generating GTF file for %s" % lsv.id)
                log.error(e.message)
    log.info("GTF files for LSVs saved in %s" % odir)
    if out_gff3:
        log.info("GFF3 files for LSVs saved in %s" % odir)
