import gzip
import urllib.parse as urllib
from collections import namedtuple
from majiq.src.constants import *
import h5py
import numpy as np
import pysam
import os
import pickle
import traceback
import sys

from majiq.grimoire.exon import detect_exons, new_exon_definition, set_exons_gc_content
from majiq.grimoire.gene import Gene, Transcript
from majiq.grimoire.junction import Junction
from majiq.src.config import Config
import majiq.src.utils as majiq_utils
from majiq.src.normalize import gc_factor_calculation
from voila.io_voila import VoilaInput
from voila.vlsv import VoilaLsv
from voila.splice_graphics import LsvGraphic

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
    majiq_config = Config()
    res = True
    if majiq_config.strand_specific:
        if (read.flag & 0x10 == 0x10 and gene_strand == '+') or (read.flag & 0x10 == 0x00 and gene_strand == '-'):
            res = True
        else:
            res = False
    return res


def _check_read(read):
    return __is_unique(read) and _match_strand(read, gg_strand)


def gc_content_per_file(args_vals, output_gc_vals, outdir):
    global gg_strand
    try:
        gc_pairs = {'GC': [], 'COV': []}
        db_f = h5py.File(get_build_temp_db_filename(outdir), 'r')
        for exp_idx, ff in args_vals:
            samfile = pysam.AlignmentFile(ff, "rb")
            for gene_name in db_f.keys():
                chromsome = db_f[gene_name].attrs['chromosome']
                gg_strand = db_f[gene_name].attrs['strand']

                for ex in db_f[gene_name]['exons']:
                    gc_val = db_f[gene_name]['exons/%s' % ex].attrs['gc_content']
                    st = db_f[gene_name]['exons/%s' % ex].attrs['start']
                    end = db_f[gene_name]['exons/%s' % ex].attrs['end']

                    if gc_val == 0 or end - st < 30:
                        continue
                    nreads = samfile.count(reference=chromsome, start=st, end=end,
                                           until_eof=False, read_callback=_check_read)
                    # if nreads > 0:
                    #     ex.set_in_data()
                    gc_pairs['GC'].append(gc_val)
                    gc_pairs['COV'].append(nreads)
            samfile.close()

            factor, meanbins = gc_factor_calculation(gc_pairs, nbins=10)
            output_gc_vals[exp_idx] = (factor, meanbins)


    except:
        traceback.print_exc()
        sys.stdout.flush()
        raise


def rnaseq_intron_retention(gne, samfl, chnk, permissive=True, nondenovo=False, logging=None):

    # filenames, gene_list, chnk, permissive=True, nondenovo=False, logging=None)
    num_bins = NUM_INTRON_BINS
    intron_list = gne.get_all_introns()
    strand = gne.get_strand()
    chrom = gne.get_chromosome()
    majiq_config = Config()
    for exon1, exon2 in intron_list:
        ex1_end = exon1.get_coordinates()[1]
        ex2_start = exon2.get_coordinates()[0]
        intron_start = ex1_end + 1
        intron_end = ex2_start - 1

        intron_len = intron_end - intron_start
        if intron_len <= 0:
            continue

        try:
            read_iter = samfl.fetch(chrom, intron_start + 8, intron_end - 8)

        except ValueError:
            # logging.info('There are no reads in %s:%d-%d' % (chrom, ex1_end, ex1_end+1))
            continue

        nchunks = 1 if intron_len <= MIN_INTRON_LEN else num_bins

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
            offset = readlen - MIN_BP_OVERLAP

            if intron_start - r_start > readlen:
                r_start = intron_start - (readlen - MIN_BP_OVERLAP*2) - 1

            if r_start < ex1_end - MIN_BP_OVERLAP:
                if junc1 is None:
                    junc1 = Junction(ex1_end, intron_start, exon1, None, gne.get_id(), retrieve=True)
                junc1.update_junction_read(nreads, r_start, gc_content, unique)

            elif (ex2_start - offset - 1) < r_start < ex2_start:
                if junc2 is None:
                    junc2 = Junction(intron_end, ex2_start, exon2, None, gne.get_id(), retrieve=True)
                junc2.update_junction_read(nreads, r_start, gc_content, unique)

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

        cov1 = junc1.get_coverage().sum()
        cov2 = junc2.get_coverage().sum()

        # intron_parts /= chunk_len

        intron_body_covered = True

        if intron_len > majiq_config.readLen:
            for ii in range(nchunks):
                num_positions = np.count_nonzero(bmap[index_list[ii][0]:index_list[ii][1]])
                nii = intron_parts[ii]
                if nii == 0:
                    val = 0
                elif num_positions == 0:
                    continue
                else:
                    val = float(nii) / num_positions
                if val < majiq_config.min_intronic_cov:
                    intron_body_covered = False
                    break

        if cov1 >= majiq_config.min_denovo and cov2 >= majiq_config.min_denovo and intron_body_covered:
            exnum = new_exon_definition(intron_start, intron_end,
                                        junc1, junc2, gne, nondenovo=nondenovo,
                                        isintron=True)
            if exnum == -1:
                continue
            logging.debug("NEW INTRON RETENTION EVENT %s, %d-%d" % (gne.get_name(), intron_start, intron_end))
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
    gne.prepare_exons()


def open_rnaseq(samfile):
    return pysam.Samfile(samfile, "rb")


def close_rnaseq(samfl):
    samfl.close()


def read_sam_or_bam(gne, samfl, counter,  h5py_file, nondenovo=False, info_msg='0', logging=None):

    junctions = []
    strt, end = gne.get_coordinates()
    j_list = gne.get_all_junctions(filter=False)
    ex_list = gne.get_exon_list()
    strand = gne.get_strand()
    chrom = gne.get_chromosome()
    majiq_config = Config()
    try:
        read_iter = samfl.fetch(chrom, strt, end, multiple_iterators=True)
       # kk = samfl.pileup(reference=chrom, start=strt, end=end)
        for read in read_iter:
            r_start = read.pos
            unique = __is_unique(read)
            if not _match_strand(read, gene_strand=strand) or r_start < strt or not unique:
                continue

            nreads = __get_num_reads(read)
            gne.add_read_count(nreads)
            is_cross, junc_list = __cross_junctions(read)

            if majiq_config.gcnorm:
                for ex_idx in range(len(ex_list)):
                    ex_start, ex_end = ex_list[ex_idx].get_coordinates()
                    if ex_start <= r_start <= ex_end:
                        ex_list[ex_idx].update_coverage(nreads)
                        break

            if not is_cross:
                continue

            nc = read.seq.count('C') + read.seq.count('c')
            ng = read.seq.count('g') + read.seq.count('G')
            gc_content = float(nc + ng) / float(len(read.seq))
            readlen = len(read.seq)
            for (junc_start, junc_end) in junc_list:
                if junc_start - r_start > readlen:
                    r_start_offset = junc_list[0][0] - r_start
                    r_start = junc_start - r_start_offset
                    if junc_start - r_start >= readlen - MIN_BP_OVERLAP or junc_start - r_start <= MIN_BP_OVERLAP:
                        continue
                elif junc_start - r_start >= readlen - MIN_BP_OVERLAP or junc_start - r_start <= MIN_BP_OVERLAP:
                    continue

                if junc_end - junc_start < MIN_JUNC_LENGTH:
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
                        jj.update_junction_read(nreads, r_start, gc_content, unique)
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
                            jnc.update_junction_read(nreads, r_start, gc_content, unique)
                            if not (junc_start, '5prime', jnc) in junctions:
                                junctions.append((junc_start, '5prime', jnc))
                                junctions.append((junc_end, '3prime', jnc))
                            junc = jnc
                            break
                            # end if (j.start) == ...
                    # end for (coord,t,j) ...

                    if junc is None:
                        '''mark a new junction '''
                        bb = gne.check_antisense_junctions_hdf5(junc_start, junc_end, h5py_file)
                        if not bb:
                            counter[4] += 1
                            junc = Junction(junc_start, junc_end, None, None, gne.get_id(), retrieve=True)
                            junc.update_junction_read(nreads, r_start, gc_content, unique)
                            junctions.append((junc_start, '5prime', junc))
                            junctions.append((junc_end, '3prime', junc))
                            # end if not found ...
                            # end for junc ...
                            #            print "JJJunctions", junctions

        if len(junctions) > 0:
            detect_exons(gne, junctions, None)
    except ValueError:
        logging.error('\t[%s]There are no reads in %s:%d-%d' % (info_msg, chrom, strt, end))
    finally:
        gne.prepare_exons()

    # logging.debug("INVALID JUNC", counter[0])
    # logging.debug("READ WRONG GENE", counter[1])
    # logging.debug("READ IN GENE", counter[2])
    # logging.debug("READ FOUND JUNC", counter[3])
    # logging.debug("READ NEW JUNC", counter[4])
    # logging.debug("READ ALL JUNC", counter[5])


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

accepted_transcripts = ['mRNA', 'transcript']
transcript_id_keys = ['ID']
gene_name_keys = ['Name', 'gene_name']
gene_id_keys = ['ID', 'gene_id']


def read_gff(filename, list_of_genes, logging=None):
    """
    :param filename: GFF input filename
    :param list_of_genes: List of genes that will be updated with all the gene_id detected on the gff file
    :param logging: logger object
    :return: :raise RuntimeError:
    """

    majiq_utils.monitor('PRE_GFF')
    majiq_config = Config()
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

            for gname_k in gene_name_keys:
                try:
                    gene_name = record.attributes[gname_k]
                    break
                except KeyError:
                    continue
            else:
                if logging is not None:
                    logging.info("Error, Gene doesn't contain one of the Name attribute "
                                 "information values: %s" % gene_name_keys)

            for gid_k in gene_id_keys:
                try:
                    gene_id = record.attributes[gid_k]
                    break
                except KeyError:
                    continue
            else:
                if logging is not None:
                    logging.info("Error, Gene doesn't contain one of the ID attribute "
                                 "information values: %s" % gene_id_keys)

            if chrom not in all_genes:
                all_genes[chrom] = {'+': [], '-': []}

            gn = Gene(gene_id, gene_name, chrom, strand, start, end)
            gn.exist_antisense_gene(all_genes[chrom])

            if gene_id in majiq_config.gene_tlb and gn != majiq_config.gene_tlb[gene_id]:
                raise RuntimeError('Two Different Genes with the same name %s' % gene_name)
            majiq_config.gene_tlb[gene_id] = gn
            list_of_genes.append(gene_id)
            all_genes[chrom][strand].append(gn)
            gene_id_dict[gene_id] = gn

        elif record.type in accepted_transcripts:
            for tid_k in transcript_id_keys:
                try:
                    transcript_name = record.attributes[tid_k]
                    break
                except KeyError:
                    continue
            else:
                if logging is not None:
                    logging.info("Error, Transcript doesn't contain one of the ID attribute "
                                 "information values: %s" % transcript_id_keys)

            try:
                parent = record.attributes['Parent']
            except KeyError:
                if logging is not None:
                    logging.info("Error, incorrect gff. mRNA %s doesn't have valid parent attribute"
                                 % transcript_name)

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
                    logging.WARNING("Error, incorrect gff. exon at line %s "
                                    "doesn't have valid mRNA %s" % (0, parent_tx_id))
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

    for chrom in all_genes.keys():
        exon_list = []
        for strand in all_genes[chrom].keys():
            for gn in all_genes[chrom][strand]:
                gn.collapse_exons()
                exon_list.extend(gn.get_exon_list())
        if majiq_config.gcnorm:
            try:
                set_exons_gc_content(chrom, exon_list)
            except RuntimeWarning:
                continue

    # majiq_utils.monitor('GC_CONTENT')
    # if majiq_config.gcnorm:
    #     get_exon_gc_content(gc_pairs, sam_list, all_genes)

    _prepare_and_dump(filename="%s/tmp/db.hdf5" % majiq_config.outDir, logging=logging)
    majiq_utils.monitor('POST_GFF')
    del all_genes


#######
# HDF5 API
#######

def _prepare_and_dump(filename, logging=None):
    majiq_config = Config()
    if logging is not None:
        logging.debug("Number of Genes in DB", len(majiq_config.gene_tlb))
    db_f = h5py.File(filename, 'w')
    for gidx, gn in enumerate(majiq_config.gene_tlb.values()):
        gn.to_hdf5(db_f)


def get_const_junctions(filename, logging=None):
    if not os.path.exists(filename):
        if logging is not None:
            logging.error('File % doesn\'t exists' % filename)
            raise UserWarning
    else:
        db_f = h5py.File(filename, 'r')
        cc = db_f[JUNCTIONS_DATASET_NAME][()]
        db_f.close()
        return np.array(cc)


def extract_lsv_summary(files):

    lsvid2idx = {}
    lsv_types = {}
    idx_junc = {}
    total_idx = 0
    simpl_juncs = []

    lsv_dict_graph = {}
    for fidx, ff in enumerate(files):
        simpl_juncs.append([[0, 0.0] for xx in idx_junc.keys()])
        data = h5py.File(ff, 'r')
        for lsvid in data['LSVs']:
            lsv = data['LSVs/%s' % lsvid]
            lsvgraph = LsvGraphic.easy_from_hdf5(data['LSVs/%s/visual' % lsvid])
            cov = data[JUNCTIONS_DATASET_NAME][lsv.attrs['coverage']]

            lsv_types[lsvid] = lsvgraph.lsv_type
            ljunc = lsvgraph.junction_ids()

            #JV
            lsv_dict_graph[lsvid] = lsvgraph

            cov = [(cov != 0).sum(axis=1), cov.sum(axis=1)]
            lsvid2idx[lsvid] = []
            for jidx, jj in enumerate(ljunc):
                try:
                    indx = idx_junc[jj]
                    simpl_juncs[fidx][indx] = [cov[0][jidx], cov[1][jidx]]
                except KeyError:
                    idx_junc[jj] = total_idx
                    indx = total_idx
                    total_idx += 1
                    simpl_juncs[fidx].append([cov[0][jidx], cov[1][jidx]])
                    [simpl_juncs[dx].append([0, 0.0]) for dx in range(fidx)]
                lsvid2idx[lsvid].append(indx)

    simpl_juncs = np.array(simpl_juncs)

    metas = read_meta_info(files)

    return lsvid2idx, lsv_types, simpl_juncs, metas, lsv_dict_graph


def load_data_lsv(path, group_name, logger=None):
    """Load data from the preprocess step. Could change to a DDBB someday"""
    data = h5py.File(path, 'r')
    lsv_cov_list = []
    lsv_info = []

    meta_info = dict()
    meta_info['group'] = group_name
    meta_info['sample_id'] = data.attrs['sample_id']
    meta_info['fitfunc'] = data.attrs['fitfunc']

    try:
        for lsvid in data['LSVs'].keys():
            lsv = data['LSVs/%s' % lsvid]
            lsv_info.append([lsv.attrs['id'], lsv.attrs['type'], lsv['visual']])
            lsv_cov_list.append(data[JUNCTIONS_DATASET_NAME][lsv.attrs['coverage']])
            sh = data[JUNCTIONS_DATASET_NAME][lsv.attrs['coverage']].shape
            if sh[0] < 2:
                print("WRONG LSV %s" % lsvid)

    except KeyError:
        logger.info("No LSVs in file")
        raise

    return meta_info, [lsv_cov_list, lsv_info]


def load_lsvgraphic_from_majiq(h5df_grp, lsv_id):
    try:
        return h5df_grp['/LSVs/%s/visual' % lsv_id]
    except KeyError:
        return None

def read_meta_info(list_of_files):
    meta = {'experiments': []}
    for fl in list_of_files:
        with h5py.File(fl, 'r') as fp :
            meta['experiments'].append(fp.attrs['sample_id'])
            try:
                if meta['genome'] != fp.attrs['genome']:
                    raise RuntimeError('Combining experiments from different genome assemblies. Exiting')
            except KeyError:
                meta['genome'] = fp.attrs['genome']
                continue

    return meta


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


def open_hdf5_file(filename, **kwargs):
    return h5py.File(filename, 'r')


def close_hdf5_file(fp):
    return fp.close()



# bootstrap files

from majiq.grimoire.lsv import quant_lsv


def get_extract_lsv_list(list_of_lsv_id, file_list):
    result = []
    fitfunc = []

    for lsv_id in list_of_lsv_id:
        lsv_cov = []
        lsv_type = None
        for fidx, fname in enumerate(file_list):

            with open_hdf5_file(fname) as data:
                if len(fitfunc) < (fidx+1):
                    fitfunc.append(data.attrs['fitfunc'])

                try:
                    if lsv_type is None:

                        lsv_type = data['LSVs/%s' % lsv_id].attrs['type']

                    assert data['LSVs/%s' % lsv_id].attrs['type'] == lsv_type, "ERROR lsv_type doesn't match for %s" % lsv_id
                    lsv_cov.append(data[JUNCTIONS_DATASET_NAME][data['LSVs/%s' % lsv_id].attrs['coverage']])
                except KeyError:
                    lsv_cov.append(None)

#        lsv_cov = np.array(lsv_cov)
        qq = quant_lsv(lsv_id, lsv_type, lsv_cov)
        result.append(qq)
    return result, fitfunc


def add_lsv_to_bootstrapfile(lsv_id, lsv_type, samples, num_exp, lock_per_file, outdir, name):

    for ii in range(num_exp):
        vals = {'samples': samples[ii], 'id': lsv_id, 'type': lsv_type}
        file_name = '%s/%s.%d.boots.hdf5' % (outdir, name, ii)
        lock_per_file[ii].acquire()
        with h5py.File(file_name, 'r+') as f:
            lsv_idx = f.attrs['lsv_idx']
            lsv_idx = boots_write(f, vals, lsv_idx)
            f.attrs['lsv_idx'] = lsv_idx
        lock_per_file[ii].release()


def create_bootstrap_file(file_list, outdir, name, m=100):
    import datetime
    for ii, ff in enumerate(file_list):
        f = h5py.File('%s/%s.%d.boots.hdf5' % (outdir, name, ii), 'w')
        f.create_dataset('junctions', (5000, m), maxshape=(None, m))
        # fill meta info
        f.attrs['sample_id'] = ff
        f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
        f.attrs['VERSION'] = VERSION
        f.attrs['lsv_idx'] = 0
        f.close()


def load_bootstrap_samples(lsv_id, file_list, weight=True):
    lsv_samples = []
    lsv_type = None
    for ii, f in enumerate(file_list):
        samples = f['junctions'][f['LSVs/%s' % lsv_id].attrs['coverage']]
        lsv_type = f['LSVs/%s' % lsv_id].attrs['type']
        if weight:
            #wght = f['weights'][f["LSVs/%s" % lsv_id].attrs['weight_idx']]
            wght = f["LSVs/%s" % lsv_id].attrs['weight']
            samples *= wght
        lsv_samples.append(samples)
    return lsv_samples, lsv_type


def store_weights_bootstrap(lsv_list, wgts, file_list, outdir, name):
    for ii, ff in enumerate(file_list):
        file_name = '%s/%s.%d.boots.hdf5' % (outdir, name, ii)
        with h5py.File(file_name, 'r+') as f:

            for idx, lsv in lsv_list.items():
                f["LSVs/%s" % lsv.id].attrs['weight'] = wgts[idx, ii]
            f.close()


def open_bootstrap_samples(num_exp, directory, name):
    result = []
    for ii in range(num_exp):
        file_name = '%s/%s.%d.boots.hdf5' % (directory, name, ii)
        result.append(h5py.File(file_name, 'r+'))
    return result


def close_bootstrap_file(file_list, outdir, name, m=100):
    for ii, ff in enumerate(file_list):
        file_name = '%s/%s.%d.boots.hdf5' % (outdir, name, ii)
        f = h5py.File(file_name, 'r+')
        lsv_idx = f.attrs['lsv_idx']
        f['junctions'].resize((lsv_idx, m))
        f.close()


def boots_write(hg_grp, vals, lsv_idx, dpsi=False):

    njunc = vals['samples'].shape[0]
    if lsv_idx + njunc > 2:
        shp = hg_grp['junctions'].shape
        shp_new = shp[0] + 5000
        hg_grp['junctions'].resize((shp_new, shp[1]))

    hg_grp['junctions'][lsv_idx:lsv_idx+njunc] = vals['samples']

    h_lsv = hg_grp.create_group("LSVs/%s" % vals['id'])
    h_lsv.attrs['id'] = vals['id']
    h_lsv.attrs['type'] = vals['type']
    h_lsv.attrs['coverage'] = hg_grp['junctions'].regionref[lsv_idx:lsv_idx + njunc]

    return lsv_idx + njunc