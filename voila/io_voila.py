import csv
import os
from collections import defaultdict
from multiprocessing import Manager, Process, Pool
from multiprocessing.queues import JoinableQueue

import h5py
import numpy as np

import vlsv
from voila import constants
from voila.hdf5 import HDF5, BinsDataSet, Psi1DataSet, Psi2DataSet
from voila.utils import utils_voila
from voila.vlsv import VoilaLsv

__author__ = 'abarrera'
import cPickle as pkl

logger = utils_voila.get_logger(__name__)


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

    def encode_metainfo(self, h, metainfo):
        if type(metainfo) is list:
            for index, mi in enumerate(metainfo):
                self.encode_metainfo(h.create_group(str(index)), mi)
        else:
            for key in metainfo:
                h.attrs[key] = metainfo[key]

    def exclude(self):
        return ['metainfo', 'lsvs']

    def decode_metainfo(self, h):
        if h.keys():
            return [self.decode_metainfo(h[key]) for key in h]
        else:
            return {key: h.attrs[key] for key in h.attrs}

    def to_hdf5(self, h):
        # bins dataset
        bins_length = sum([len(lsv.bins) for lsv in self.lsvs])
        bins_width = len(self.lsvs[0].bins[0])
        BinsDataSet(h, bins_length, bins_width)

        try:
            # psi1 dataset
            psi1_length = sum([len(lsv.psi1) for lsv in self.lsvs])
            Psi1DataSet(h, psi1_length)

            # psi2 dataset
            psi2_length = sum([len(lsv.psi2) for lsv in self.lsvs])
            Psi2DataSet(h, psi2_length)
        except TypeError:
            # this is probably psi data
            pass

        # metainfo
        self.encode_metainfo(h.create_group('metainfo'), self.metainfo)

        # lsvs
        lsv_grp = h.create_group('lsvs')
        for lsv in self.lsvs:
            lsv.to_hdf5(lsv_grp)

        super(VoilaInput, self).to_hdf5(h)

    def from_hdf5(self, h):
        # metainfo
        self.metainfo = self.decode_metainfo(h['metainfo'])

        # lsvs
        self.lsvs = [VoilaLsv((), None).from_hdf5(h['lsvs'][lsv_id]) for lsv_id in h['lsvs']]

        return super(VoilaInput, self).from_hdf5(h)


def voila_input_from_hdf5(hdf5_filename, logger):
    """
    Create VoilaInput object from HDF5 file.  This will process each of the VoilaLsvs in their own thread using the
    Producer Consumer design pattern.
    :param hdf5_filename: HDF5 filename string
    :param logger: instance of logger
    :return: VoilaInput object
    """

    def worker():
        with h5py.File(hdf5_filename, 'r', swmr=True) as h:
            while True:
                id = queue.get()
                manager_lsvs.append(VoilaLsv((), None).from_hdf5(h['lsvs'][id]))
                queue.task_done()

    def producer():
        with h5py.File(hdf5_filename, 'r', swmr=True) as h:
            for id in h['lsvs']:
                queue.put(id)

    def metainfo():
        with h5py.File(hdf5_filename, 'r', swmr=True) as h:
            manager_dict['metainfo'] = voila_input.decode_metainfo(h['metainfo'])

    logger.info('Loading {0}.'.format(hdf5_filename))

    voila_input = VoilaInput()

    queue = JoinableQueue()

    manager = Manager()
    manager_dict = manager.dict()
    manager_lsvs = manager.list()

    metainfo_proc = Process(target=metainfo)
    metainfo_proc.daemon = True
    metainfo_proc.start()

    producer_proc = Process(target=producer)
    producer_proc.daemon = True
    producer_proc.start()

    pool = Pool(None, worker)

    metainfo_proc.join()
    producer_proc.join()
    queue.join()

    pool.close()
    queue.close()

    voila_input.lsvs = list(manager_lsvs)
    voila_input.metainfo = manager_dict['metainfo']

    return voila_input


def dump_voila_input(voila_input, target, logger=None, protocol=-1):
    import os
    base_path, name = os.path.split(target)
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    try:
        pkl.dump(voila_input, target, protocol)
    except pkl.PickleError, e:
        if logger:
            logger.error("Dumping file %s in %s:\n\t%s." % (voila_input, target, e.message), exc_info=1)


def load_voila_input(voila_input_file, logger=None):
    try:
        with open(voila_input_file, 'rb') as vif:
            voila_input = pkl.load(vif)
            return voila_input
    except pkl.PickleError, e:
        if logger:
            logger.error("Loading the file %s:\n\t%s." % (voila_input_file, e.message), exc_info=1)


def load_dpairs(pairwise_dir, majiq_output, logger):
    """
    Load pairwise files from MAJIQ analysis.

    :param str pairwise_dir: directory containing pairwise comparisons produced by MAJIQ.
    :param majiq_output: parsed data from old_majiq.
    :param logger: logger instance.
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
                lmajiq_pairs[idx1][idx2] = utils_voila.get_lsv_delta_exp_data(pairwise_file,
                                                                              show_all=True,
                                                                              gene_name_list=lsv_names,
                                                                              logger=logger)
            except IOError:
                pass
    return lmajiq_pairs, group1_name, group2_name


def write_tab_output(input_parsed):
    """
    Create tab-delimited output file summarizing all the LSVs detected and quantified with MAJIQ.

    :param output_dir: output directory for the file.
    :param output_html: name for the output html file used to create a *.txt version.
    :param majiq_output: parsed data from old_majiq.
    :param type_summary: type of analysis performed.
    :param logger: logger instance.
    :param pairwise_dir: whether pairwise comparisons are included or not.
    :param threshold: minimum change considered as significant (in deltapsi analysis).
    :return: nothing.
    """
    if input_parsed.type_summary == constants.COND_TABLE:
        cond_table_tab_output(input_parsed)
    else:
        tab_output(input_parsed)


def cond_table_tab_output(input_parsed):
    majiq_output = input_parsed.majiq_output
    lsvs = majiq_output['lsvs']
    sample_names = majiq_output['sample_names']

    input_parsed.logger.info('Creating cond-table TSV...')

    tsv_file = os.path.join(input_parsed.output_dir, input_parsed.output_html.split('.html')[0] + '.tsv')
    input_parsed.logger.info(tsv_file)

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
    logger = input_parsed.logger
    pairwise_dir = input_parsed.pairwise_dir
    majiq_output = input_parsed.majiq_output
    threshold = input_parsed.threshold
    type_summary = input_parsed.type_summary

    ofile_str = "%s%s.%s" % (output_dir, output_html.rsplit('.html', 1)[0], constants.EXTENSION)
    tlb_categx = {'A5SS': 'prime5', 'A3SS': 'prime3', 'Num. Junctions': 'njuncs', 'Num. Exons': 'nexons', 'ES': 'ES'}

    logger.info("Creating Tab-delimited output file in %s..." % ofile_str)

    if pairwise_dir:
        # In deltapsi, add columns with pairwise comparisons between group members
        logger.info("Load pairwise comparison files from %s..." % pairwise_dir)
        lmajiq_pairs, group1_name, group2_name = load_dpairs(pairwise_dir, majiq_output, logger=logger)

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
            psi_headers = ['%s E(PSI)' % majiq_output['meta_exps'][0][0]['group'],
                           '%s E(PSI)' % majiq_output['meta_exps'][1][0]['group']]
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
                lline.extend([llsv.lsv_graphic.get_name(), gene, llsv.get_id()])
                lexpected = []
                lconfidence = []
                lexpecs_psi1 = []
                lexpecs_psi2 = []
                for i, bins in enumerate(llsv.get_bins()):
                    if 'delta' in type_summary:
                        lexpected.append(str(-llsv.get_excl_incl()[i][0] + llsv.get_excl_incl()[i][1]))
                        # lconfidence.append(str(utils_voila.get_prob_delta_psi_greater_v(bins, float(lexpected[-1]), threshold)))
                        lconfidence.append(str(vlsv.matrix_area(np.array(bins), threshold, collapsed_mat=True).sum()))
                        lexpecs_psi1.append('%.3f' % vlsv.get_expected_psi(np.array(llsv.psi1[i])))
                        lexpecs_psi2.append('%.3f' % vlsv.get_expected_psi(np.array(llsv.psi2[i])))
                    else:
                        lexpected.append(repr(llsv.get_means()[i]))
                        lconfidence.append(repr(llsv.get_variances()[i]))

                lline.append(';'.join(lexpected))
                lline.append(';'.join(lconfidence))
                if 'delta' in type_summary:
                    lline.append(';'.join(lexpecs_psi1))
                    lline.append(';'.join(lexpecs_psi2))

                lline.append(llsv.get_type())
                lline.append(repr(llsv.get_categories()[tlb_categx['A5SS']]))
                lline.append(repr(llsv.get_categories()[tlb_categx['A3SS']]))
                lline.append(repr(llsv.get_categories()[tlb_categx['ES']]))
                lline.append(repr(llsv.get_categories()[tlb_categx['Num. Junctions']]))
                lline.append(repr(llsv.get_categories()[tlb_categx['Num. Exons']]))
                lline.append(str(int(np.any([junc.get_type() == 1 for junc in llsv.lsv_graphic.get_junctions()]))))

                lline.append(llsv.lsv_graphic.get_chrom())
                lline.append(llsv.lsv_graphic.get_strand())

                lline.append(';'.join(
                    ['-'.join(str(c) for c in junc.get_coords()) for junc in llsv.lsv_graphic.get_junctions()]))
                lline.append(
                    ';'.join(['-'.join(str(c) for c in exon.get_coords()) for exon in llsv.lsv_graphic.get_exons()]))

                try:
                    lline.append(';'.join(
                        ['|'.join([str(c) for c in exon.get_alt_starts()]) for exon in llsv.lsv_graphic.get_exons()]))
                    lline.append(';'.join(
                        ['|'.join([str(c) for c in exon.get_alt_ends()]) for exon in llsv.lsv_graphic.get_exons()]))
                except TypeError:
                    pass

                lline.append(
                    ';'.join([repr(exon.coords) for exon in llsv.lsv_graphic.get_exons() if exon.intron_retention]))

                if pairwise_dir:
                    llpairwise = []
                    for idx1 in range(len(lmajiq_pairs)):
                        for idx2 in range(len(lmajiq_pairs[0])):
                            lpairwise = []
                            if gene in lmajiq_pairs[idx1][idx2]['genes_dict']:
                                for llsv_tmp in lmajiq_pairs[idx1][idx2]['genes_dict'][gene]:
                                    if llsv_tmp[0].get_id() == llsv.get_id():
                                        lsv_pair = llsv_tmp[0]
                                        break
                                else:
                                    logger.warning("LSV %s present in deltagroup but missing in %s." %
                                                   (llsv.get_id(), "%s_%d_%s_%d" % (group1_name, idx1 + 1,
                                                                                    group2_name, idx2 + 1)))
                                    lpairwise.append('N/A')
                                    continue
                                for iway in range(len(llsv.get_bins())):
                                    lpairwise.append(str(sum(lsv_pair.get_excl_incl()[iway])))
                            else:
                                lpairwise.append('N/A')
                            llpairwise.append(';'.join(lpairwise))
                    lline.extend(llpairwise)
                if 'voila_links' in majiq_output.keys():
                    summary_path = majiq_output['voila_links'][llsv.get_gene_name()]
                    if not os.path.isabs(summary_path):
                        summary_path = "%s/%s/%s" % (os.getcwd(), output_dir, summary_path)
                    lline.append(constants.URL_COMPOSITE % (summary_path, llsv.get_gene_name()))
                ofile.write(constants.DELIMITER.join(lline))
                ofile.write('\n')

    logger.info("Delimited output file successfully created in: %s" % ofile_str)


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
                    if fields[0] not in filter_genes and fields[1] not in filter_genes: continue

                if filter_lsvs:
                    if fields[2].upper() not in filter_lsvs: continue

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
        # Update the number of changing
        # lsvs_dict[lsv_idx]['nchangs'] = np.sum([1 for ff in lsvs_dict[lsv_idx]['expecs'] if np.any(np.array([abs(fff) for fff in ff]) > thres_change)])

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
    :param output_dir: output directory for the file.
    :param majiq_output: parsed data from old_majiq.
    :param logger: logger instance.
    :param out_gff3:
    :return: nothing.
    """
    logger = input_parsed.logger
    majiq_output = input_parsed.majiq_output
    output_dir = input_parsed.output_dir

    if input_parsed.type_summary == constants.COND_TABLE:
        logger.info('Skipping generating gff3 format.')
        return

    logger.info("Saving LSVs files in gff3 format ...")
    if 'genes_dict' not in majiq_output or len(majiq_output['genes_dict']) < 1:
        logger.warning("No gene information provided. Genes files are needed to calculate the gff3 files.")
        return

    header = "##gff-version 3"

    odir = output_dir + "/static/doc/lsvs"
    utils_voila.create_if_not_exists(odir)
    for gkey, gvalue in majiq_output['genes_dict'].iteritems():
        for lsv_dict in gvalue:
            lsv = lsv_dict
            if type(lsv_dict) == dict:
                lsv = lsv_dict['lsv']
            lsv_file_basename = "%s/%s" % (odir, lsv.get_id())

            try:
                lsv_gff3_str = lsv.get_gff3(logger=logger)
                utils_voila.gff2gtf(lsv_gff3_str.split('\n'), "%s.gtf" % lsv_file_basename)
                if out_gff3:
                    gff_file = "%s.gff3" % (lsv_file_basename)
                    with open(gff_file, 'w') as ofile:
                        ofile.write(header + "\n")
                        ofile.write(lsv_gff3_str + "\n")
            except UnboundLocalError, e:
                logger.warning("problem generating GTF file for %s" % lsv.get_id())
                logger.error(e.message)
    logger.info("GTF files for LSVs saved in %s" % odir)
    if out_gff3:
        logger.info("GFF3 files for LSVs saved in %s" % odir)
