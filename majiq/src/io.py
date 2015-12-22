import pickle
import random
import gc
import os
import sys
from collections import namedtuple
import gzip
import urllib
import math

import numpy as np
import pysam
import ConfigParser

import majiq.src.config as majiq_config
from voila.io_voila import VoilaInput
from voila.vlsv import VoilaLsv
from majiq.grimoire.gene import Gene, Transcript
import majiq.grimoire.exon as majiq_exons
from majiq.grimoire.junction import Junction


def create_if_not_exists(my_dir, logger=False):
    """Create a directory path if it does not exist"""
    try:
        if logger:
            logger.info("\nCreating directory %s..." % my_dir)
        os.makedirs(my_dir)
    except OSError:
        if logger:
            logger.info("\nDirectory %s already exists..." % my_dir)


def load_bin_file(filename, logger=None):
    if not os.path.exists(filename):
        if logger:
            logger.error('Path %s for loading does not exist' % filename)
        return

    fop = open(filename, 'rb')

    fast_pickler = pickle.Unpickler(fop)
    # fast_pickler.fast = 1
    data = fast_pickler.load()
    fop.close()
    return data


def dump_bin_file(data, filename):
    with open(filename, 'wb') as ofp:
        fast_pickler = pickle.Pickler(ofp, protocol=2)
        # fast_pickler.fast = 1
        fast_pickler.dump(data)


# pickle.dump(data, protocol=2)


def __cross_junctions(read):
    """
      This part will parse the jI from STAR
      # Column 17: jI:B:I,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1- based)
      # jI:B:I,-1 means that the read doesn't cross any junction
    """

    jlist = []
    cross = False
    try:
        list_junc = read.opt('jI')
        if len(list_junc) > 1 and list_junc[0] != -1:
            for idx in range(0, len(list_junc), 2):
                junc_start = int(list_junc[idx]) - 1
                junc_end = int(list_junc[idx + 1]) + 1
                jlist.append((junc_start, junc_end))
            # end for idx ...
            cross = True
            # end else ...
    except KeyError:
        # if len(jlist) != 0: print "STAR ::",jlist
        # print"THIS IS NOT a WELL defined STAR output"
        off = 0
        for op, num in read.cigar:
            if op in [0, 5, 6, 7, 8]:
                off += num
            elif op in [1, 5]:
                off += 0
            elif op == 2:
                off += num
            elif op == 3:
                jlist.append((read.pos + off, read.pos + off + num + 1))
                off += num
                cross = True
                # if len(jlist) !=0 : print "NOSTAR:", jlist, read.cigar

    return cross, jlist


def __is_unique(read):
    unique = True
    try:
        if read.opt('NH') > 1:
            unique = False
    except KeyError:
        if read.flag & 0x100 == 1:
            unique = False
    return unique


def __get_num_reads(read):
    try:
        nreads = int(read.opt('HI'))
    except KeyError:
        nreads = 1
    # return 1
    return nreads


def _match_strand(read, gene_strand):
    res = True
    if majiq_config.strand_specific:
        if (read.flag & 0x10 == 0x10 and gene_strand == '+') or (read.flag & 0x10 == 0x00 and gene_strand == '-'):
            res = True
        else:
            res = False
    return res


def count_mapped_reads(filename, exp_idx):
    stats = pysam.flagstat(filename)
    mapped_reads = int(stats[2].split()[0])
    majiq_config.num_mapped_reads[exp_idx] = mapped_reads


# def is_neg_strand(read):
#     res = False
#     if read.flag & 0x10 == 0x10:
#         # print "FLAG",read.flag
#         res = True
#
#     if majiq_config.strand_specific:
#         res = not res
#
#     return res


def get_junc_from_list(coords, list_elem):
    res = None
    for xx in list_elem:
        newcoord = xx.get_coordinates()
        if newcoord[0] == coords[0] and newcoord[1] == coords[1]:
            res = xx
            break
    return res


def rnaseq_intron_retention(filenames, gene_list, chnk, permissive=True, nondenovo=False, logging=None):
    samfile = [pysam.Samfile(xx, "rb") for xx in filenames]
    num_bins = 10
    for gne in gene_list:
        intron_list = gne.get_all_introns()
        strand = gne.get_strand()
        chrom = gne.get_chromosome()
        for exon1, exon2 in intron_list:
            ex1_end = exon1.get_coordinates()[1]
            ex2_start = exon2.get_coordinates()[0]
            intron_start = ex1_end + 1
            intron_end = ex2_start - 1


            intron_len = intron_end - intron_start
            if intron_len <= 0:
                continue

            if intron_len <= 1000:
                nchunks = 1
            else:
                nchunks = num_bins

            # we want to take just the middle part not the reads that are crossing the junctions
            # since 8 is the overlapping number of nucleotites we accept, the inner part is the
            # real intron size - (readlen-8)/*start part*/ - (readlen-8)/*end part*/

            chunk_len = intron_len / nchunks

            bmap = np.ones(shape=intron_len, dtype=bool)
            index_list = []
            for ii in range(nchunks):
                start = ii * chunk_len
                end = min(intron_len, (ii + 1) * chunk_len)
                index_list.append((start, end))

            intron_parts = np.zeros(shape=nchunks, dtype=np.float)
            junc1 = None
            junc2 = None

            for name, ind_list in majiq_config.tissue_repl.items():
                n_exp = 0
                repl_thresh = int(math.ceil(len(ind_list) * 0.5))

                for idx, exp_index in enumerate(ind_list):
                    try:
                        read_iter = samfile[exp_index].fetch(chrom, intron_start + 8, intron_end - 8)
                    except ValueError:
                        # logging.info('There are no reads in %s:%d-%d' % (chrom, ex1_end, ex1_end+1))
                        continue

                    for read in read_iter:
                        is_cross, junc_list = __cross_junctions(read)
                        if not _match_strand(read, gene_strand=strand):
                            continue

                        unique = __is_unique(read)
                        r_start = read.pos
                        nreads = __get_num_reads(read)

                        if not unique:
                            intron_idx = r_start - (ex1_end + 1)
                            if not (0 <= intron_idx <= intron_len):
                                continue
                            bmap[intron_idx] = False
                            continue

                        if is_cross:
                            jvals = [xx for xx, yy in junc_list if not (yy < intron_start or xx > intron_end)]
                            if len(jvals) > 0:
                                continue

                        nc = read.seq.count('C') + read.seq.count('c')
                        ng = read.seq.count('g') + read.seq.count('G')
                        gc_content = float(nc + ng) / float(len(read.seq))
                        readlen = len(read.seq)
                        offset = readlen - 8

                        if intron_start - r_start > readlen:
                            r_start = intron_start - (readlen - 16) - 1

                        if r_start < ex1_end - 8:
                            if junc1 is None:
                                junc1 = Junction(ex1_end, intron_start, exon1, None, gne, readN=0)
                            junc1.update_junction_read(exp_index, nreads, r_start, gc_content, unique)

                        elif (ex2_start - offset - 1) < r_start < ex2_start:
                            if junc2 is None:
                                junc2 = Junction(intron_end, ex2_start, exon2, None, gne, readN=0)
                            junc2.update_junction_read(exp_index, nreads, r_start, gc_content, unique)

                        else:
                            # section 3
                            intron_idx = r_start - (ex1_end + 1)
                            rel_start = intron_idx / chunk_len
                            indx = -1 if rel_start > nchunks else rel_start
                            if not bmap[intron_idx]:
                                bmap[intron_idx] = True
                            intron_parts[indx] += nreads

                    if junc1 is None or junc2 is None:
                        continue

                    cov1 = junc1.get_coverage(exp_index).sum()
                    cov2 = junc2.get_coverage(exp_index).sum()

                    # intron_parts /= chunk_len

                    intron_body_covered = True
                    comp_chunk = nchunks
                    intron_covered = 0

                    if intron_len > readlen:
                        for ii in range(nchunks):
                            # for ii in intron_parts:
                            # num_positions = np.count_nonzero(bmap[index_list[ii][0]:index_list[ii][1]])
                            num_positions = np.count_nonzero(bmap[index_list[ii][0]:index_list[ii][1]])
                            nii = intron_parts[ii]
                            if nii == 0:
                                val = 0
                            elif num_positions == 0:
                                continue
                            else:
                                val = float(nii) / num_positions
                            if val < majiq_config.MIN_INTRON:
                                intron_body_covered = False
                                break

                    if cov1 >= majiq_config.min_denovo and cov2 >= majiq_config.min_denovo and intron_body_covered:
                        n_exp += 1

                if n_exp >= repl_thresh:
                    exnum = majiq_exons.new_exon_definition(intron_start, intron_end,
                                                            None, junc1, junc2, gne, nondenovo=nondenovo,
                                                            isintron=True)
                    if exnum == -1:
                        for exp_index in ind_list:
                            if not junc2 is None:
                                junc2.reset_coverage(exp_index)
                            if not junc1 is None:
                                junc1.reset_coverage(exp_index)
                    else:

                        junc1.add_donor(exon1)
                        for ex in exon1.exonRead_list:
                            st, end = ex.get_coordinates()
                            if end == junc1.get_coordinates()[0]:
                                ex.add_5prime_junc(junc1)
                                break

                        junc2.add_acceptor(exon2)
                        for ex in exon2.exonRead_list:
                            st, end = ex.get_coordinates()
                            if st == junc2.get_coordinates()[1]:
                                ex.add_3prime_junc(junc2)
                                break

                        if exnum == 1:
                            logging.info("NEW INTRON RETENTION EVENT %s, %d-%d" % (gne.get_name(),
                                                                                   intron_start,
                                                                                   intron_end))
                else:
                    for exp_index in ind_list:
                        if not junc2 is None:
                            junc2.reset_coverage(exp_index)
                        if not junc1 is None:
                            junc1.reset_coverage(exp_index)
        gne.prepare_exons()

    for ss in samfile:
        ss.close()
    gc.collect()
    return


def read_sam_or_bam(filenames, gene_list, chnk, nondenovo=False, logging=None):
    counter = [0] * 6
    samfile = [pysam.Samfile(xx, "rb") for xx in filenames]
    temp_ex = []
    non_unique_num = 0
    skip_gene = 0
    non_skip = 0

    for gne in gene_list:
        junctions = []
        strt, end = gne.get_coordinates()
        j_list = gne.get_all_junctions()
        ex_list = gne.get_exon_list()
        strand = gne.get_strand()
        chrom = gne.get_chromosome()

        for exp_index in range(len(filenames)):

            #readlen = config.readLen[exp_index]
            try:
                read_iter = samfile[exp_index].fetch(chrom, strt, end)
            except ValueError:
                logging.info('There are no reads in %s:%d-%d' % (chrom, strt, end))
                continue
            for read in read_iter:

                if not _match_strand(read, gene_strand=strand):
                    continue
                unique = __is_unique(read)
                if not unique:
                    non_unique_num += 1
                    continue
                nreads = __get_num_reads(read)
                gne.add_read_count(nreads, exp_index)
                is_cross, junc_list = __cross_junctions(read)
                r_start = read.pos
                if r_start < strt or r_start > end:
                    continue

                for ex_idx in range(len(ex_list)):
                    ex_start, ex_end = ex_list[ex_idx].get_coordinates()
                    if ex_start <= r_start <= ex_end:
                        ex_list[ex_idx].update_coverage(exp_index, nreads)
                        temp_ex.append(ex_list[ex_idx])
                        break

                if not is_cross:
                    continue
                nc = read.seq.count('C') + read.seq.count('c')
                ng = read.seq.count('g') + read.seq.count('G')
                gc_content = float(nc + ng) / float(len(read.seq))
                readlen = len(read.seq)
                for (junc_start, junc_end) in junc_list:
                    if junc_start - r_start > readlen:
                        r_start = junc_start - (readlen - 16) - 1
                    elif junc_start - r_start >= readlen - 8 or junc_start - r_start <= 8:
                        continue

                    if junc_end - junc_start < 10:
                        counter[0] += 1
                        continue

                    found = False

                    for jj in j_list:
                        (j_st, j_ed) = jj.get_coordinates()
                        if j_st > junc_start or (j_st == junc_start and j_ed > junc_end):
                            break
                        elif j_st < junc_start or (j_st == junc_start and j_ed < junc_end):
                            continue
                        elif junc_start == j_st and junc_end == j_ed:
                            ''' update junction and add to list'''
                            found = True
                            counter[3] += 1
                            jj.update_junction_read(exp_index, nreads, r_start, gc_content, unique)
                            if not (junc_start, '5prime', jj) in junctions:
                                junctions.append((junc_start, '5prime', jj))
                                junctions.append((junc_end, '3prime', jj))
                            break
                            # end elif junc_start == ...
                    # end for jj in j_list

                    if not found:
                        if nondenovo:
                            continue

                        ''' update junction and add to list'''
                        junc = None
                        for (coord, t, jnc) in junctions:
                            if jnc.start == junc_start and jnc.end == junc_end:
                                jnc.update_junction_read(exp_index, nreads, r_start, gc_content, unique)
                                if not (junc_start, '5prime', jnc) in junctions:
                                    junctions.append((junc_start, '5prime', jnc))
                                    junctions.append((junc_end, '3prime', jnc))
                                junc = jnc
                                break
                                # end if (j.start) == ...
                        # end for (coord,t,j) ...
                        if junc is None:
                            '''mark a new junction '''
                            bb = gne.check_antisense_junctions(junc_start, junc_end)
                            if not bb:
                                counter[4] += 1
                                junc = Junction(junc_start, junc_end, None, None, gne, readN=nreads)
                                junc.update_junction_read(exp_index, nreads, r_start, gc_content, unique)
                                junctions.append((junc_start, '5prime', junc))
                                junctions.append((junc_end, '3prime', junc))
                                # end if not found ...
                                # end for junc ...
                                #            print "JJJunctions", junctions
        if len(junctions) > 0:
            majiq_exons.detect_exons(gne, junctions, None)
        gne.prepare_exons()

    for ss in samfile:
        ss.close()
    gc.collect()

    logging.debug("INVALID JUNC", counter[0])
    logging.debug("READ WRONG GENE", counter[1])
    logging.debug("READ IN GENE", counter[2])
    logging.debug("READ FOUND JUNC", counter[3])
    logging.debug("READ NEW JUNC", counter[4])
    logging.debug("READ ALL JUNC", counter[5])
    logging.debug("Non Unique", non_unique_num)

    logging.debug("Skipped genes without exons", skip_gene)
    logging.debug(" Non skipped", non_skip)
    return


def read_bed_pcr(filename, list_genes):
    input_f = open(filename, 'r')
    readlines = input_f.readlines()
    alt_exon = []
    pre_chrom = ''
    gene_list = {}
    lnum = 0
    while lnum < len(readlines):

        event = {}
        rl = readlines[lnum]
        if rl.startswith('#'):
            continue
        tab = rl.strip().split()
        event['name'] = tab[3]
        event['chrom'] = tab[0]
        event['strand'] = tab[5]

        score = tab[4]
        lnum += 1

        strand = event['strand']

        if event['chrom'] != pre_chrom:
            try:
                gene_list = list_genes[event['chrom']]
            except KeyError:
                continue

            pre_chrom = event['chrom']
            idx = {'+': 0, '-': 0}
        name = event['name']

        exon_start = int(tab[1])
        exon_end = int(tab[2])

        while idx[strand] < len(gene_list[strand]):
            gn = gene_list[strand][idx[strand]]
            (g_start, g_end) = gn.get_coordinates()
            if exon_end < g_start:
                break
            elif exon_start > g_end:
                idx[strand] += 1
                continue
            ex = gn.exist_exon(exon_start, exon_end)
            if ex is None:
                break
            ex.set_pcr_score(name, score, alt_exon)

            break


gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)


def __parse_gff_attributes(attribute_string):
    """Parse the GFF3 attribute column and return a dict
    :param attribute_string:
    """  #
    if attribute_string == ".":
        return {}
    ret = {}
    for attribute in attribute_string.split(";"):
        key, value = attribute.split("=")
        key = urllib.unquote(key)
        if key in ret:
            key = 'extra_%s' % key
            if not key in ret:
                ret[key] = []
            ret[key].append(urllib.unquote(value))
        else:
            ret[key] = urllib.unquote(value)
    return ret


def __parse_gff3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    # Parse with transparent decompression
    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            # If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            # Normalize data
            normalized_info = {
                "seqid": None if parts[0] == "." else urllib.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib.unquote(parts[0]),
                "type": None if parts[2] == "." else urllib.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.unquote(parts[7]),
                "attributes": __parse_gff_attributes(parts[8])
            }
            #Alternatively, you can emit the dictionary here, if you need mutabwility:
            #    yield normalized_info
            yield GFFRecord(**normalized_info)


def _prepare_and_dump_old(genes, logging=None):
    n_genes = 0
    for chrom in genes.keys():
        temp_ex = []
        for strand, gg in genes[chrom].items():
            n_genes += len(gg)
            genes[chrom][strand] = sorted(gg)
            for gene in genes[chrom][strand]:
                gene.collapse_exons()
                temp_ex.extend(gene.get_exon_list())
        if not logging is None:
            logging.info("Calculating gc_content chromosome %s........." % chrom)
        majiq_exons.set_exons_gc_content(chrom, temp_ex)
        gc.collect()
        temp_dir = "%s/tmp/%s" % (majiq_config.outDir, chrom)
        create_if_not_exists(temp_dir)
        # ipdb.set_trace()
        # objgraph.show_most_common_types(limit=20)
        if not logging is None:
            logging.info("Creating temporal annotation %s" % chrom)
        fname = '%s/annot_genes.pkl' % temp_dir
        dump_bin_file(genes[chrom], fname)

    tmp_chrom = "%s/tmp/chromlist.pkl" % majiq_config.outDir
    dump_bin_file(genes.keys(), tmp_chrom)
    if not logging is None:
        logging.debug("Number of Genes", n_genes)


def __annot_dump(nthrd, temp_ex, lsv_list, logging=None):
    for chrom, ex_list in temp_ex.items():
        majiq_exons.set_exons_gc_content(chrom, ex_list)
    gc.collect()
    temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, nthrd)
    create_if_not_exists(temp_dir)
    if not logging is None:
        logging.info("Creating temporal annotation chunk %s (%d genes)" % (nthrd, len(lsv_list)))
    fname = '%s/annot_genes.pkl' % temp_dir
    dump_bin_file(lsv_list, fname)


def __get_overlaped(gn, temp_ex, dumped_genes):
    lsv_list = []
    num_gns = 0
    over_genes = gn.get_overlapped_genes()
    if not over_genes is None:
        for extra_gn_id in over_genes:
            if extra_gn_id in dumped_genes:
                continue
            extra_gn = majiq_config.gene_tlb[extra_gn_id]
            extra_gn.collapse_exons()
            temp_ex.extend(extra_gn.get_exon_list())
            lsv_list.append(extra_gn)
            dumped_genes.append(extra_gn_id)
            a, b = __get_overlaped(extra_gn, temp_ex, dumped_genes)
            lsv_list.extend(a)
            num_gns += (b + 1)

    return lsv_list, num_gns


def _prepare_and_dump(logging=None):
    list_genes = sorted(majiq_config.gene_tlb.values())
    if not logging is None:
        logging.debug("Number of Genes", len(list_genes))

    chunk_size = len(list_genes) / majiq_config.num_final_chunks
    temp_ex = {}
    nthrd = 0
    csize = chunk_size
    lsv_list = []

    dumped_genes = []

    for gidx, gn in enumerate(list_genes):
        if gn.get_id() in dumped_genes:
            continue
        gn.collapse_exons()
        csize -= 1

        chrom = gn.get_chromosome()
        if not chrom in temp_ex:
            temp_ex[chrom] = []
        temp_ex[chrom].extend(gn.get_exon_list())
        lsv_list.append(gn)
        dumped_genes.append(gn.get_id())
        a, b = __get_overlaped(gn, temp_ex[chrom], dumped_genes)
        csize -= b
        lsv_list.extend(a)

        if csize <= 0:
            __annot_dump(nthrd, temp_ex, lsv_list, logging)

            lsv_list = []
            csize = chunk_size - 1
            nthrd += 1
            temp_ex = {}
            if not chrom in temp_ex:
                temp_ex[chrom] = []
            temp_ex[chrom].extend(gn.get_exon_list())

    if len(lsv_list) > 0:
        __annot_dump(nthrd, temp_ex, lsv_list, logging)


def read_gff(filename, pcr_filename, nthreads, logging=None):
    """
    :param filename: GFF input filename
    :param pcr_filename: BED file name with the PCR validations
    :param logging: logger object
    :return: :raise RuntimeError:
    """
    all_genes = {}
    gene_id_dict = {}
    trcpt_id_dict = {}
    last_end = {}
    for record in __parse_gff3(filename):
        chrom = record.seqid
        strand = record.strand
        start = record.start
        end = record.end

        if record.type == 'gene':
            gene_name = record.attributes['Name']
            gene_id = record.attributes['ID']

            if not chrom in all_genes:
                all_genes[chrom] = {'+': [], '-': []}

            gn = Gene(gene_id, gene_name, chrom, strand, start, end)
            gn.exist_antisense_gene(all_genes[chrom])

            if gene_id in majiq_config.gene_tlb and gn != majiq_config.gene_tlb[gene_id]:
                raise RuntimeError('Two Different Genes with the same name %s' % gene_name)
            majiq_config.gene_tlb[gene_id] = gn
            all_genes[chrom][strand].append(gn)
            gene_id_dict[record.attributes['ID']] = gn

        elif record.type == 'mRNA' or record.type == 'transcript':
            transcript_name = record.attributes['ID']
            parent = record.attributes['Parent']
            try:
                gn = gene_id_dict[parent]
                trcpt = Transcript(transcript_name, gn, start, end)
                gn.add_transcript(trcpt)
                trcpt_id_dict[record.attributes['ID']] = trcpt
                last_end[record.attributes['ID']] = (None, None)
            except KeyError:
                if not logging is None:
                    logging.info("Error, incorrect gff. mRNA %s doesn't have valid gene %s"
                                 % (transcript_name, parent))
                raise

        elif record.type == 'exon':
            parent_tx_id = record.attributes['Parent']
            try:
                parent_tx = trcpt_id_dict[parent_tx_id]
                gn = parent_tx.get_gene()
                txex = gn.new_annotated_exon(start, end, parent_tx)
                parent_tx.add_exon(txex)

            except KeyError:
                if not logging is None:
                    logging.info("Error, incorrect gff. exon %s doesn't have valid mRNA %s" % (record.attributes['ID'],
                                                                                               parent_tx_id))
                    # end elif
    # end for
    for tid, trcpt in trcpt_id_dict.items():

        exon_list = trcpt.prepare_exon_list()
        gn = trcpt.get_gene()
        pre_end = None
        pre_txex = None
        for ex in exon_list:
            start, end = ex.get_coordinates()
            junc = gn.new_annotated_junctions(pre_end, start, trcpt)
            ex.add_3prime_junc(junc)
            if not pre_txex is None:
                pre_txex.add_5prime_junc(junc)
            pre_end = end
            pre_txex = ex

        junc = gn.new_annotated_junctions(pre_end, None, trcpt)
        pre_txex.add_5prime_junc(junc)

    # end for
    import sys

    print sys.getsizeof(all_genes)

    _prepare_and_dump(logging)
    if pcr_filename is not None:
        read_bed_pcr(pcr_filename, all_genes)

    chr_list = all_genes.keys()
    del all_genes
    return chr_list


# Quantifier i/o
def load_data_lsv(path, group_name, logger=None):
    """Load data from the preprocess step. Could change to a DDBB someday"""
    data = pickle.load(open(path))
    lsv_cov_list = []
    #lsv_gc = []
    lsv_info = []
    const_info = []
    num_pos = data[1][0].junction_list.shape[1]

    meta_info = data[0]
    meta_info['group'] = group_name
    for lsv in data[1]:
        try:
            lsv_info.append([lsv.coords, lsv.id, lsv.type, 0, lsv.visual])
        except AttributeError, e:
            lsv_info.append([lsv.coords, lsv.id, lsv.type, 0])

        cov = lsv.junction_list.toarray()
        lsv_cov_list.append(cov)
        #gc = lsv.gc_factor.toarray()
        #lsv_gc.append(gc)

    clist = random.sample(data[2], min(5000, len(data[2])))
    const_list = np.zeros(shape=(len(clist), num_pos), dtype=np.dtype('int'))
    #const_gc = np.zeros(shape=(len(clist), num_pos), dtype=np.dtype('float'))
    for cidx, const in enumerate(clist):
        const_info.append(const.id)
        const_list[cidx, :] = const.coverage.toarray()
        #const_gc[cidx, :] = const.gc_factor.toarray()

    return meta_info, [lsv_cov_list, lsv_info], [const_list, const_info]


def dump_lsvs_voila(pickle_path, posterior_matrix, lsvs_info, meta_info, psi_list1=None, psi_list2=None):
    """Create VoilaLSVs objects readable by voila."""
    vlsvs = []
    psi1, psi2 = None, None
    for ii, bins in enumerate(posterior_matrix):
        lsv_graphic = lsvs_info[ii][-1]
        if psi_list1:
            psi1, psi2 = psi_list1[ii], psi_list2[ii]
        vlsvs.append(VoilaLsv(bins, lsv_graphic=lsv_graphic, psi1=psi1, psi2=psi2))

    pickle.dump(VoilaInput(vlsvs, meta_info), open(pickle_path, 'w'))


def read_multi_dpsi_conf(filename):
    config = ConfigParser.ConfigParser()
    config.read(filename)
    # TODO: check if filename exists
    list_of_deltas = ConfigSectionMap(config, "deltas")
    list_of_groups = ConfigSectionMap(config, "groups")
    info = ConfigSectionMap(config, "info")

    groups = {}
    files_dict = {}
    for kk, vv in list_of_groups.items():
        fls_lst = vv.split(',')
        groups[kk.lower()] = []
        for fls in fls_lst:
            fl_name = os.path.split(fls)[1]
            groups[kk].append(fl_name)
            abs_fls = "%s/%s" % (info['path'], fl_name)
            files_dict[fl_name] = abs_fls

    deltas = []
    print list_of_deltas
    for kk, vv in list_of_deltas.items():
        deltas.append((kk.lower(), vv.lower()))
    return groups, files_dict, deltas



def ConfigSectionMap(Config, section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1