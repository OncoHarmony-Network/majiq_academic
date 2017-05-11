from __future__ import division

import errno
import fnmatch
import os
import shutil
from collections import defaultdict

import numpy as np

from voila.utils.voila_log import voila_log
from voila.vlsv import VoilaLsv


def expected_dpsi(bins):
    return sum(np.array(bins) * np.arange(-1 + 1. / len(bins), 1., 2. / len(bins)))


def get_prob_delta_psi_greater_v(bins, expected, V=.2):
    """Calculate probability of delta psi outside the acceptable area"""
    bins = np.array(bins)
    step = 2.0 / bins.size
    left = 0
    right = bins.size * 2 - 1
    for i, w in enumerate(np.arange(-1 + step / 2, 1, step)):
        if not left and w > (expected - abs(expected * V)):
            left = i - 1
        if right == bins.size * 2 and w > (expected + abs(expected * V)):
            right = i
    return np.sum(bins[:left] + np.sum(bins[right:]))


def lsvs_to_gene_dict(voila_lsvs, metainfo, threshold=.2, show_all=False):
    log = voila_log()
    genes_dict = defaultdict(list)

    for i, vlsv in enumerate(voila_lsvs):

        if np.any(np.isnan(vlsv.bins)):
            log.warning("LSV %s bins contain NaNs" % vlsv.get_id())
            continue

        if vlsv.is_delta_psi() and not show_all and not vlsv.is_lsv_changing(threshold):
            continue

        gene_id = vlsv.lsv_id.split(':')[0]

        if vlsv.is_delta_psi():
            genes_dict[gene_id].append({
                'lsv': vlsv,
                'psi1': VoilaLsv(vlsv.psi1, lsv_graphic=None),
                'psi2': VoilaLsv(vlsv.psi2, lsv_graphic=None)
            })
        else:
            genes_dict[gene_id].append(vlsv)

    return {'genes_dict': genes_dict, 'meta_exps': metainfo}


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
    min_sst = [20] * 20
    res = tab[0]

    ss_list = set()
    ex_list = set()
    for tt in tab[1:]:
        pp = tt.split('e')
        ss2 = int(pp[0])
        min_sss = min(min_sss, ss2)
        try:
            ss3 = int(pp[1].split('.')[1])
        except IndexError:
            ss3 = 1
        ext = int(pp[1].split('.')[0])

        min_sst[ext] = min(min_sst[ext], ss3)
    for tt in tab[1:]:
        tab2 = tt.split('e')
        new_ss = int(tab2[0]) - min_sss + 1
        tab3 = tab2[1].split('.')
        if len(tab3) == 1:
            tab3.append(1)
        new_sst = int(tab3[1]) - min_sst[int(tab3[0])] + 1
        ss_list.add(new_ss)
        ex_list.add('%s.%s' % (tab3[0], new_sst))
        # ss_list += '%s,'%new_ss
        # ex_list += '%s.%s,'%(tab3[0],tab3[1])
    # res += '|%se%s.%s'%(new_ss,tab3[0],tab3[1])
    ss = ','.join([str(x) for x in sorted(ss_list)])
    exs = ','.join([str(x) for x in sorted(ex_list)])
    res += '|%se%s|%s' % (ss, exs, len(tab[1:]))

    return res


def list_files_or_dir(file_or_dir_list, prefix='*', suffix='*', containing='*'):
    if type(file_or_dir_list) != list: return [file_or_dir_list]
    files = []
    for file_or_dir in file_or_dir_list:
        if os.path.isdir(file_or_dir):
            for root, dirnames, filenames in os.walk(file_or_dir):
                for filename in fnmatch.filter(filenames, '%s*%s*%s' % (prefix, containing, suffix)):
                    files.append(os.path.join(root, filename))
        else:
            files.append(file_or_dir)
    return files


def create_if_not_exists(my_dir):
    """Create a directory path if it does not exist"""
    try:
        os.makedirs(my_dir)
    except OSError as err:
        if err.errno == errno.EEXIST:
            pass


def gff2gtf(gff_f, out_f):
    """Parse a GFF file created by MAJIQ and create a GTF"""

    mrna = None
    with open(out_f, 'w') as wfile:
        if type(gff_f) == list:
            gff = gff_f
        else:
            gff = open(gff_f)

        for gff_l in gff:
            gff_fields = gff_l.strip().split()
            if len(gff_fields) < 3:
                continue
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


def debug(text):
    print text
    return ''
