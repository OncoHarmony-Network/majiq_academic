from __future__ import division
import fnmatch
import logging
import os
import sys
from collections import defaultdict
import json
import shutil
import errno
from voila.io_voila import load_voila_input

from voila.lsv import Lsv
from voila.splice_graphics import ExonGraphic, JunctionGraphic, GeneGraphic, LsvGraphic
from voila.vlsv import extract_bins_info, VoilaLsv


try:
    import cPickle as pkl
except ImportError:
    try:
        import pickle as pkl
    except ImportError:
        print "[Error] :: Neither pickle nor cPickle are installed. Please, check python dependencies."
        import sys
        sys.exit(1)

try:
    import numpy as np
except ImportError:
    print "[Error] :: Numpy not installed. Please, check python dependencies."
    sys.exit(1)


class PickleEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, np.ndarray):
            return list(obj)
        if isinstance(obj, tuple):
            return list(obj)
        if isinstance(obj, np.int64):
            return int(obj)
        if isinstance(obj, Lsv):
            return obj.to_JSON(PickleEncoder)

        return json.JSONEncoder.default(self, obj)


class LsvGraphicEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, np.ndarray):
            return list(obj)
        if isinstance(obj, tuple):
            return list(obj)
        if isinstance(obj, np.int64):
            return int(obj)
        if isinstance(obj, ExonGraphic):
            return obj.to_JSON(PickleEncoder)
        if isinstance(obj, JunctionGraphic):
            return obj.to_JSON(PickleEncoder)
        if isinstance(obj, GeneGraphic):
            return obj.to_JSON(PickleEncoder)
        if isinstance(obj, LsvGraphic):
            return obj.to_JSON(PickleEncoder)

        return json.JSONEncoder.default(self, obj)



def expected_dpsi(bins):
    return sum(np.array(bins) * np.arange(-1+1./len(bins), 1.,2./len(bins)))


def get_prob_delta_psi_greater_v(bins, expected, V=.2):
    """Calculate probability of delta psi outside the acceptable area"""
    bins = np.array(bins)
    step = 2.0 / bins.size
    left = 0
    right = bins.size*2 - 1
    for i, w in enumerate(np.arange(-1+step/2, 1, step)):
        if not left and w > (expected - abs(expected*V)):
            left = i-1
        if right == bins.size*2 and w > (expected + abs(expected*V)):
            right = i
    return (np.sum(bins[:left]) + np.sum(bins[right:]))


def get_lsv_single_exp_data(majiq_bins_file, confidence, gene_name_list=None, lsv_types=None, logger=None):
    """
    Create a dictionary to summarize the information from majiq output file.
    """
    majiq_data = None
    try:
        majiq_data = pkl.load(open(majiq_bins_file, 'rb'))
    except pkl.PickleError:
        logger.error("Pickle could not load the file. Please, check that the file %s is in Pickle format." % majiq_bins_file, exc_info=1)

    except IOError:
        logger.error("%s doesn't exists." % majiq_bins_file, exc_info=1)
        sys.exit(1)

    lsv_list = []

    genes_dict = defaultdict(list)

    nofilter_genes = not gene_name_list and not lsv_types
    if gene_name_list is None:
        gene_name_list = []

    for i, lsv_meta in enumerate(majiq_data[1]):
        # Determine whether include the LSV or not, based on selected filters
        if nofilter_genes or str(lsv_meta[1]).split(':')[0] in gene_name_list or lsv_meta[4].name.upper() in gene_name_list or lsv_meta[2] in lsv_types:
            # print lsv_meta[0], lsv_meta[1], lsv_meta[2]

            bins_array_list = majiq_data[0][i]

            # In 1-way LSVs, create the additional bins set for commodity
            if len(bins_array_list) == 1:
                bins_array_list.append(bins_array_list[-1][::-1])

            try:
                # metadata.append([lsv_meta[0], lsv_meta[1], lsv_meta[2]]) #collapse_lsv(lsv_meta[2])])
                lsv_list.append(Lsv(bins_array_list, lsv_meta, confidence))
                # lsv_list[-1].sort_bins(lsv_meta[4].strand)

            except ValueError, e:
                logger.warning("%s produced an error:\n%s. Skipped." % (bins_array_list, e))
                continue

            genes_dict[str(lsv_meta[1]).split(':')[0]].append([lsv_list[-1], lsv_meta])

    return {'event_list':   lsv_list,
            'metadata':     majiq_data[1],
            'genes_dict':   genes_dict,
            'meta_exps':    majiq_data[2]}



def get_lsv_delta_exp_data(voila_input_file, confidence=.95, threshold=.2, show_all=False, gene_name_list=None, logger=None):
    """
    Load lsv delta psi pickle file. It contains a list with 2 elements:
        [0] List with LSV bins matrices
        [1] List with info per LSV

    :param voila_input_file:
    :param metadata_post:
    :param confidence:
    :param threshold:
    :param logger:
    @return: dictionary
    """

    voila_input = load_voila_input(voila_input_file, logger=logger)
    meta_info = voila_input.metainfo

    genes_dict = defaultdict(list)
    lsv_list = voila_input.lsvs

    for i, vlsv in enumerate(lsv_list):
        include_lsv = show_all
        gene_name_id = str(vlsv.lsv_graphic.id).split(':')[0]
        gene_name = vlsv.lsv_graphic.name.upper()

        if not gene_name_list or gene_name_id in gene_name_list or gene_name in gene_name_list:
            collapsed_bins, excl_inc_perc_list, include_lsv = extract_bins_info(vlsv.bins, threshold, include_lsv)
            if not include_lsv: continue
            vlsv.set_bins_info(collapsed_bins, confidence)

            means = []
            excl_incl = []
            for b in vlsv.get_bins():
                means.append(expected_dpsi(b))
                if means[-1] <= 0:
                    excl_incl.append([-means[-1], 0])
                else:
                    excl_incl.append([0, means[-1]])
            vlsv.set_means(means)
            vlsv.set_excl_incl(excl_incl)
            # lsv_o.sort_bins(lsv_info[i][4].strand)

            # lsv_list.append(lsv)
            if len(vlsv.psi1)<2:
                vlsv.psi1.append(vlsv.psi1[-1][::-1])
            if len(vlsv.psi2)<2:
                vlsv.psi2.append(vlsv.psi2[-1][::-1])

            genes_dict[gene_name_id].append({'lsv': vlsv,
                                             'psi1': VoilaLsv(vlsv.psi1, lsv_graphic=None).set_bins_info(vlsv.psi1),
                                             'psi2': VoilaLsv(vlsv.psi2, lsv_graphic=None).set_bins_info(vlsv.psi2)})

    # logger.info("Number of genes added: %d" % len(genes_dict.keys()))

    return {'genes_dict': genes_dict,
            'meta_exps':  meta_info}


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:  # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        if exc.errno == errno.EEXIST:  # Static folder exists
            shutil.rmtree(dst)
            copyanything(src, dst)
        else:
            raise


def collapse_lsv(lsv_type):
    tab = lsv_type.split('|')
    if len(tab) < 3:
        return lsv_type + '|1'
    min_sss = 20
    min_sst = [20]*20
    res = tab[0]

    ss_list = set()
    ex_list = set()
    for tt in tab[1:]:
        pp= tt.split('e')
        ss2 = int(pp[0])
        min_sss = min(min_sss, ss2)
        try:
            ss3 = int(pp[1].split('.')[1])
        except IndexError, e:
            ss3 = 1
        ext = int(pp[1].split('.')[0])

        min_sst[ext] = min(min_sst[ext], ss3)
    for tt in tab[1:]:
        tab2 = tt.split('e')
        new_ss = int(tab2[0])-min_sss +1
        tab3 = tab2[1].split('.')
        if len(tab3) == 1:
            tab3.append(1)
        new_sst = int(tab3[1])-min_sst[int(tab3[0])] + 1
        ss_list.add(new_ss)
        ex_list.add( '%s.%s'%(tab3[0],new_sst))
        #ss_list += '%s,'%new_ss
        #ex_list += '%s.%s,'%(tab3[0],tab3[1])
    #res += '|%se%s.%s'%(new_ss,tab3[0],tab3[1])
    ss = ','.join([str(x) for x in sorted(ss_list)])
    exs = ','.join([str(x) for x in sorted(ex_list)])
    res += '|%se%s|%s'%(ss,exs, len(tab[1:]))

    return res


def list_files_or_dir(file_or_dir_list, prefix='*', suffix='*', containing='*'):

    if type(file_or_dir_list) != list: return [file_or_dir_list]
    files = []
    for file_or_dir in file_or_dir_list:
        if os.path.isdir(file_or_dir):
            for root, dirnames, filenames in os.walk(file_or_dir):
                for filename in fnmatch.filter(filenames, '%s*%s*%s' % (prefix, containing, suffix)):
                    files.append(os.path.join(root, filename))
            # for file in os.listdir(file_or_dir):
            #     if not suffix or file.endswith(suffix):
            #         files.append(file_or_dir+'/'+file)
        else:
            files.append(file_or_dir)
    return files


def get_logger(logger_name, silent=False, debug=False):
    """
    Returns a logger instance. verbose = False will silence the logger, debug will give
    more information intended for debugging purposes.
    """
    logging_format = "%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s"
    logging.basicConfig(filename=logger_name, format=logging_format)
    logger = logging.getLogger(logger_name)
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    if debug:
        ch.setLevel(logging.DEBUG)
    elif not silent:
        ch.setLevel(logging.INFO)
    else:
        ch.setLevel(logging.WARNING)

    formatter = logging.Formatter("%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


def create_if_not_exists(my_dir, logger=False):
    """Create a directory path if it does not exist"""
    try:
        if logger:
            logger.info("\nCreating directory %s..." % my_dir)
        os.makedirs(my_dir)
    except OSError:
        if logger:
            logger.info("\nDirectory %s already exists..." % my_dir)


def collapse_matrix(matrix):
    """Collapse the diagonals probabilities in 1-D and return them"""
    collapse = []
    #FOR TEST matrix = array([[0, 1, 2, 3, 4, 500], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5], [100, 1, 2, 3, 4, 5], ])

    matrix_corner = matrix.shape[0]+1
    for i in xrange(-matrix_corner, matrix_corner):
        collapse.append(np.diagonal(matrix, offset=i).sum())

    return np.array(collapse)


def gff2gtf(gff_f, out_f=None):
    """Parse a GFF file created by MAJIQ and create a GTF"""

    mrna = None
    with open(out_f, 'w') as wfile:
        with open(gff_f) as gff:
            for gff_l in gff:
                gff_fields = gff_l.strip().split()
                if len(gff_fields)<3: continue
                if gff_fields[2] == 'mRNA':
                    if mrna:
                        to_gtf(wfile, seq_name, source, gene, mrna, start_trans, end_trans, strand, exon_l, frame_l)
                    exon_l = []
                    frame_l = []
                    ids = gff_fields[8].split(';Parent=')
                    gene = ids[1].split(';')[0]
                    mrna = ids[0][5:]
                    seq_name = gff_fields[0]
                    source = gff_fields[1]
                    start_trans = gff_fields[3]
                    end_trans = gff_fields[4]
                    strand = gff_fields[6]
                    len_frame = 0

                if gff_fields[2] == 'exon':
                    exon_l.append([gff_fields[3], gff_fields[4]])
                    frame_l.append((3 - (len_frame % 3)) % 3)
                    len_frame = int(gff_fields[4]) - int(gff_fields[3])
        to_gtf(wfile, seq_name, source, gene, mrna, start_trans, end_trans, strand, exon_l, frame_l)


def to_gtf(wfile, seq_name, source, gene, mrna, start_trans, end_trans, strand, exon_l, frame_l):
    sscore = "0"
    # Iterate over each exon
    exonorcds_list = []
    for i, exon in enumerate(exon_l):
        exonorcds_list.append("\t".join([seq_name, source, "%s", exon[0], exon[1], sscore, strand,
                                         str(frame_l[i]), "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mrna)]))

    if strand == '+':
        first_codon = "\t".join([seq_name, source, "start_codon", start_trans, str(int(start_trans) + 2), sscore,
                                 strand, ".", "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mrna)])
        last_codon = "\t".join([seq_name, source, "stop_codon", str(int(end_trans) + 1), str(int(end_trans) + 3),
                                sscore, strand, ".", "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mrna)])

    else:
        last_codon = "\t".join([seq_name, source, "start_codon", str(int(end_trans) - 2), end_trans, sscore,
                                strand, ".", "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mrna)])
        first_codon = "\t".join([seq_name, source, "stop_codon", str(int(start_trans) - 3), str(int(start_trans) - 1),
                                 sscore, strand, ".", "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mrna)])

    wfile.write(first_codon)
    for eCDS in exonorcds_list:
        wfile.write(eCDS % "CDS")
        wfile.write(eCDS % "exon")
    wfile.write(last_codon)
