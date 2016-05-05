import gzip
import urllib
from collections import namedtuple
import h5py
import math
import numpy as np
import majiq.src.config as majiq_config
from exon import detect_exons, new_exon_definition
from gene import Gene, Transcript
from junction import Junction


# READING BAM FILES
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


def rnaseq_intron_retention(gne, samfile_list, chnk, permissive=True, nondenovo=False, logging=None):

    # filenames, gene_list, chnk, permissive=True, nondenovo=False, logging=None)
    num_bins = 10
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
                    read_iter = samfile_list[exp_index].fetch(chrom, intron_start + 8, intron_end - 8)
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
                exnum = new_exon_definition(intron_start, intron_end,
                                            None, junc1, junc2, gne, nondenovo=nondenovo,
                                            isintron=True)
                if exnum == -1:
                    pass
                    # for exp_index in ind_list:
                    #     if junc2 is not None:
                    #         junc2.reset_coverage(exp_index)
                    #     if junc1 is not None:
                    #         junc1.reset_coverage(exp_index)
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
                        logging.debug("NEW INTRON RETENTION EVENT %s, %d-%d" % (gne.get_name(),
                                                                               intron_start,
                                                                               intron_end))
            else:
                pass
                # for exp_index in ind_list:
                #     if not junc2 is None:
                #         junc2.reset_coverage(exp_index)
                #     if not junc1 is None:
                #         junc1.reset_coverage(exp_index)
    gne.prepare_exons()


def read_sam_or_bam(gne, samfile_list, chnk, counter, nondenovo=False, info_msg='0', logging=None):

    junctions = []
    strt, end = gne.get_coordinates()
    j_list = gne.get_all_junctions()
    ex_list = gne.get_exon_list()
    strand = gne.get_strand()
    chrom = gne.get_chromosome()
    for exp_index, samfl in enumerate(samfile_list):

        try:
            read_iter = samfl.fetch(chrom, strt, end)
        except ValueError:
            logging.error('\t[%s]There are no reads in %s:%d-%d' % (info_msg, chrom, strt, end))
            continue
        for read in read_iter:
            r_start = read.pos
            unique = __is_unique(read)
            if not _match_strand(read, gene_strand=strand) or r_start < strt or not unique:
                continue

            nreads = __get_num_reads(read)
            gne.add_read_count(nreads, exp_index)
            is_cross, junc_list = __cross_junctions(read)

            for ex_idx in range(len(ex_list)):
                ex_start, ex_end = ex_list[ex_idx].get_coordinates()
                if ex_start <= r_start <= ex_end:
                    ex_list[ex_idx].update_coverage(exp_index, nreads)
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

                if not found and not nondenovo:
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
        detect_exons(gne, junctions, None)
    gne.prepare_exons()

    logging.debug("INVALID JUNC", counter[0])
    logging.debug("READ WRONG GENE", counter[1])
    logging.debug("READ IN GENE", counter[2])
    logging.debug("READ FOUND JUNC", counter[3])
    logging.debug("READ NEW JUNC", counter[4])
    logging.debug("READ ALL JUNC", counter[5])


# ANNOTATION DB FUNCTIONS

gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)


def __parse_gff_attributes(attribute_string):
    """
    Parse the GFF3 attribute column and return a dict
    :param attribute_string:
    """
    if attribute_string == ".":
        return {}
    ret = {}
    for attribute in attribute_string.split(";"):
        key, value = attribute.split("=")
        key = urllib.unquote(key)
        if key in ret:
            key = 'extra_%s' % key
            if key not in ret:
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
            # Alternatively, you can emit the dictionary here, if you need mutabwility:
            #    yield normalized_info
            yield GFFRecord(**normalized_info)


def read_gff(filename, list_of_genes, logging=None):
    """
    :param filename: GFF input filename
    :param list_of_genes: List of genes that will be updated with all the gene_id detected on the gff file
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

            if chrom not in all_genes:
                all_genes[chrom] = {'+': [], '-': []}

            gn = Gene(gene_id, gene_name, chrom, strand, start, end)
            gn.exist_antisense_gene(all_genes[chrom])

            if gene_id in majiq_config.gene_tlb and gn != majiq_config.gene_tlb[gene_id]:
                raise RuntimeError('Two Different Genes with the same name %s' % gene_name)
            majiq_config.gene_tlb[gene_id] = gn
            list_of_genes.append(gene_id)
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
                if logging is not None:
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

                    logging.WARNING("Error, incorrect gff. exon at line %s doesn't have valid mRNA %s" % (0,
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
            if pre_txex is not None:
                pre_txex.add_5prime_junc(junc)
            pre_end = end
            pre_txex = ex

        junc = gn.new_annotated_junctions(pre_end, None, trcpt)
        pre_txex.add_5prime_junc(junc)

    _prepare_and_dump(filename="%s/tmp/db.hdf5" % majiq_config.outDir, logging=logging)
    del all_genes


#######
# HDF5 API
#######

def _prepare_and_dump(filename, logging=None):

    if logging is not None:
        logging.debug("Number of Genes in DB", len(majiq_config.gene_tlb))
    db_f = h5py.File(filename, 'w', compression='gzip', compression_opts=9)
    for gidx, gn in enumerate(majiq_config.gene_tlb.values()):
        gn.collapse_exons()

        gn.to_hdf5(db_f)
