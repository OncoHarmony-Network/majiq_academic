import os
import logging
import sys
from multiprocessing import Pool
import random
from contextlib import contextmanager as ctx

import numpy as np

from majiq.src.normalize import gc_factor_calculation
import majiq.src.config as majiq_config
from majiq.grimoire.junction import MajiqJunction
import majiq.src.io as majiq_io
from voila.splice_graphics import ExonGraphic
from voila.splice_graphics import GeneGraphic
from voila.splice_graphics import JunctionGraphic
from voila import constants as voila_const
import resource


def monitor(msg):
    print "MONITOR", msg, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000, 'MB'
    sys.stdout.flush()


def create_if_not_exists(my_dir, logger=False):
    """Create a directory path if it does not exist"""
    try:
        if logger:
            logger.info("\nCreating directory %s..." % my_dir)
        os.makedirs(my_dir)
    except OSError:
        if logger:
            logger.info("\nDirectory %s already exists..." % my_dir)


def get_logger(logger_name, silent=False, debug=False): 
    """
    Returns a logger instance. verbose = False will silence the logger, debug will give 
    more information intended for debugging purposes.
    """
    logging_format = "%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s"
    logging.basicConfig(filename="%s" % logger_name, format=logging_format)
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

"""
splicegraph generation functions
"""

def generate_visualization_output(allgenes, temp_dir):
    # gene_list = {}
    for name, ind_list in majiq_config.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):
            # gene_list[mglobals.exp_list[exp_idx]] = []
            gene_list = []
            for gg in allgenes:
                junc_list = []
                junc_l = []
                alt_empty_ends = []
                alt_empty_starts = []
                for jj in gg.get_all_junctions():

                    cc = jj.get_coordinates()
                    if jj.get_donor() is None:
                        alt_empty_ends.append(cc[1])
                        continue
                    if jj.get_acceptor() is None:
                        alt_empty_starts.append(cc[0])
                        continue
                    num_reads = jj.get_read_num(exp_idx)
                    if jj.is_annotated() and num_reads == 0:
                        if (jj.get_read_num(-1) - num_reads) > 0:
                            jtype = voila_const.JUNCTION_TYPE_DB_OTHER_RNASEQ
                        else:
                            jtype = voila_const.JUNCTION_TYPE_DB
                    elif jj.is_annotated() and num_reads > 0:
                        jtype = voila_const.JUNCTION_TYPE_DB_RNASEQ
                    elif not jj.is_annotated() and num_reads > majiq_config.MINREADS:
                        jtype = voila_const.JUNCTION_TYPE_RNASEQ
                    else:
                        jtype = voila_const.JUNCTION_TYPE_RNASEQ
                        # continue

                    ir_type = None
                    if jj.get_donor().is_intron():
                        ir_type = voila_const.IR_TYPE_END
                    elif jj.get_acceptor().is_intron():
                        ir_type = voila_const.IR_TYPE_START

                    junc_l.append(jj.get_coordinates())
                    junc_list.append(JunctionGraphic(cc, jtype, num_reads, transcripts=jj.get_transcript_list(),
                                                     ir=ir_type))

                junc_l = np.asarray(junc_l)
                exon_list = []
                #print gg.get_id(),
                #print gg.get_exon_list()
                for ex in gg.get_exon_list():
                    cc = ex.get_coordinates()
                    a3 = []
                    alt_start = []
                    for ss3 in set(ex.ss_3p_list):
                        if ss3 in alt_empty_starts:
                            alt_start.append(ss3)
                            # continue
                        for jidx, jjl in enumerate(junc_l):
                            if ss3 == jjl[1]:
                                a3.append(jidx)

                    a5 = []
                    alt_ends = []
                    for ss5 in set(ex.ss_5p_list):
                        if ss5 in alt_empty_starts:
                            alt_ends.append(ss5)
                            # continue
                        for jidx, jjl in enumerate(junc_l):
                            if ss5 == jjl[0]:
                                a5.append(jidx)

                    ex_reads = ex.get_total_read_num(exp_idx)
                    if ex.is_miss_start():
                        visual_type = voila_const.EXON_TYPE_MISSING_START
                    elif ex.is_miss_end():
                        visual_type = voila_const.EXON_TYPE_MISSING_END
                    elif ex.is_annotated() and ex_reads == 0.0:
                        visual_type = voila_const.EXON_TYPE_DB
                    elif ex.is_annotated() and ex_reads > 0.0:
                        visual_type = voila_const.EXON_TYPE_DB_RNASEQ
                    elif not ex.is_annotated() and ex_reads > 0.0:
                        visual_type = voila_const.EXON_TYPE_RNASEQ
                    else:
                        visual_type = voila_const.EXON_TYPE_RNASEQ
                    # continue
                    extra_coords = []
                    if ex.is_annotated():
                        if ex.start < ex.db_coord[0]:
                            extra_coords.append([ex.start, ex.db_coord[0] - 1])
                        if ex.end > ex.db_coord[1]:
                            extra_coords.append([ex.db_coord[1] + 1, ex.end])
                    eg = ExonGraphic(a3, a5, cc, type_exon=visual_type, coords_extra=extra_coords,
                                     intron_retention=ex.get_ir(), alt_starts=alt_start, alt_ends=alt_ends)
                    exon_list.append(eg)
                ggraph = GeneGraphic(id=gg.get_id(), name=gg.get_name(), strand=gg.get_strand(), exons=exon_list,
                                     junctions=junc_list, chrom=gg.get_chromosome())
                # gene_list[mglobals.exp_list[exp_idx]].append(ggraph)
                gene_list.append(ggraph)

            filename = '%s/%s.splicegraph.pkl' % (temp_dir, majiq_config.exp_list[exp_idx])
            majiq_io.dump_bin_file(gene_list, filename)


"""
Majiq file generation
"""
def prepare_lsv_table(lsv_list, non_as, temp_dir):

    #out_temp = dict()
    for name, ind_list in majiq_config.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):
            majiq_table_as = [0] * len(lsv_list[name])
            majiq_table_nonas = [0] * len(non_as[name])

            for iix, lsv in enumerate(lsv_list[name]):
                majiq_table_as[iix] = lsv.to_majiqLSV(exp_idx)
            for jix, jn in enumerate(non_as[name]):
                majiq_table_nonas[jix] = MajiqJunction(jn, exp_idx)

            out_temp = (majiq_table_as, majiq_table_nonas)
            fname = "%s/%s.majiq.pkl" % (temp_dir, majiq_config.exp_list[exp_idx])
            majiq_io.dump_bin_file(out_temp, fname)


def merge_and_create_majiq_file(exp_idx, pref_file):

    """
    :param exp_idx: Index of experiment in config file
    :param pref_file: Prefix for the majiq name
    """

    experiment = majiq_config.exp_list[exp_idx]
    all_visual = []
    as_table = []
    nonas_table = []
    for chnk in range(majiq_config.num_final_chunks):
        temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
        temp_filename = '%s/%s.splicegraph.pkl' % (temp_dir, experiment)
        if os.path.exists(temp_filename):
            visual_gene_list = majiq_io.load_bin_file(temp_filename)
            all_visual.append(visual_gene_list)
    fname = '%s/%s%s.splicegraph' % (majiq_config.outDir, pref_file, experiment)
    visual = np.concatenate(all_visual)
    majiq_io.dump_bin_file(visual, fname)
    del all_visual
    del visual

    for chnk in range(majiq_config.num_final_chunks):
        temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
        filename = "%s/%s.majiq.pkl" % (temp_dir, experiment)
        if os.path.exists(filename):
            temp_table = majiq_io.load_bin_file(filename)
            as_table.append(temp_table[0])
            nonas_table.append(temp_table[1])

    if len(as_table) == 0:
        return

    info = dict()
    info['experiment'] = experiment
    # info['GC_bins'] = majiq_config.gc_bins[exp_idx]
    # info['GC_bins_val'] = majiq_config.gc_bins_val[exp_idx]
    info['genome'] = majiq_config.genome
    info['num_reads'] = majiq_config.num_mapped_reads[exp_idx]

    at = np.concatenate(as_table)
    for lsv in at:
        lsv.set_gc_factor(exp_idx)


    nat = np.concatenate(nonas_table)
    clist = random.sample(nat, min(5000, len(nat)))
    for jnc in clist:
        jnc.set_gc_factor(exp_idx)

    fname = '%s/%s%s.majiq' % (majiq_config.outDir, pref_file, experiment)
    majiq_io.dump_bin_file((info, at, clist), fname)


def print_junc_matrices(mat, tlb=None, fp=None):
    if not fp is None:
        out = open('./junc_matrix.tab', 'a+')
    else:
        out = sys.stdout
    out.write("\n=== BEGIN %s === \n\n" % id)
    N, M = mat.shape
    header = [0]*N
    if not tlb is None:
        out.write("Nan\t")
        for ex, (p1, p2) in tlb.items():
            for n in p1:
                out.write("%d\t" % (ex+1))
            for n in p2:
#                header[nid] = "%d:%d\t"%(ex,nid)
                header[n] = "%d" % (ex+1)
    out.write("\n")
    for ii in np.arange(N):
        if not tlb is None:
            out.write("%s\t" % header[ii])
        for jj in np.arange(M):
            val = mat[ii, jj]
            out.write("%s\t" % val)
        out.write("\n")
    out.write("\n=== END %s === \n\n" % id)
    if fp is None:
        out.close()


def get_validated_pcr_lsv(candidates, out_dir):

    pcr_list = []
    print "get_validated_pcr_lsv", len(candidates[0])
    for lsv in candidates[0]:
        if not lsv.has_pcr_score():
            continue
        alt_coord = lsv.exon.get_pcr_candidate()
        score = lsv.get_pcr_score()
        for jidx, jj in enumerate(lsv.junctions):
            if lsv.is_Ssource():
                excoord = jj.get_acceptor().get_coordinates()
            else:
                excoord = jj.get_donor().get_coordinates()
            if excoord[1] > alt_coord[0] and excoord[0] < alt_coord[1]:
                name = "%s#%s" % (lsv.id, jidx)
                pcr_lsv = [lsv.exon.get_pcr_name(), name, score]
                pcr_list.append(pcr_lsv)
                print "PCR", ' '.join(pcr_lsv)
    fname = '%s/pcr.pkl' % out_dir
    majiq_io.dump_bin_file(pcr_list, fname)


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


def gff2gtf(gff_f, out_f=None):
    """Parse a GFF file created by MAJIQ and create a GTF"""

    mrna = None
    with file_or_stdout(out_f) as wfile:
        with open(gff_f) as gff:
            for gff_l in gff:
                gff_fields = gff_l.strip().split()
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


@ctx
def file_or_stdout(file_name):
    if file_name is None:
        yield sys.stdout
    else:
        with open(file_name, 'w') as out_file:
            yield out_file


def gather_files(out_dir, prefix='', gff_out=None, pcr_out=None, nthreads=1, logger=None):

    #GATHER
    logger.info("Gather outputs")
    if prefix != '':
        prefix = '%s.' % prefix

    if majiq_config.gcnorm:
        gc_factor_calculation(10)

    if nthreads > 1:
        nthr = min(nthreads, 4)
        pool = Pool(processes=nthr, maxtasksperchild=1)

    for name, ind_list in majiq_config.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):
            if nthreads > 1:
                pool.apply_async(merge_and_create_majiq_file, [exp_idx, prefix])
            else:
                merge_and_create_majiq_file(exp_idx, prefix)

    if nthreads > 1:
        pool.close()
        pool.join()

    if not gff_out is None:
        logger.info("Gather lsv and generate gff")
        fp = open('%s/%s' % (out_dir, gff_out), 'w+')
        for chnk in range(majiq_config.num_final_chunks):
            temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
            yfile = '%s/temp_gff.pkl' % temp_dir
            if not os.path.exists(yfile):
                continue
            gff_list = majiq_io.load_bin_file(yfile)
            for gff in gff_list:
                fp.write("%s\n" % gff)
        fp.close()

    if not pcr_out is None:
        logger.info("Gather pcr results")
        fp = open('%s/pcr_match.tab' % majiq_config.outDir, 'w+')
        for chnk in range(majiq_config.num_final_chunks):
            temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
            yfile = '%s/pcr.pkl' % temp_dir
            if not os.path.exists(yfile):
                continue
            pcr_l = majiq_io.load_bin_file(yfile)
            for pcr in pcr_l:
                fp.write("%s\n" % pcr)
        fp.close()


#ANALYSIS FUNCTIONS

def analyze_denovo_junctions(genes, output):

    denovo_list = [[] for xx in range(majiq_config.num_experiments)]
    annot_list = [[] for xx in range(majiq_config.num_experiments)]

    for gg in genes:
        jlist = gg.get_all_junctions()
        for jj in jlist:
            for tissue, list_idx in majiq_config.tissue_repl.items():
                for exp_n in list_idx:
                    if jj.is_annotated():
                        annot_list[exp_n].append(jj)
                    else:
                        denovo_list[exp_n].append(jj)

    majiq_io.dump_bin_file([majiq_config.tissue_repl, annot_list, denovo_list], output)


def histogram_for_exon_analysis(genes, output):

    denovo_list = []
    annotated_list = []
    #for strand, gglist in genes.items():
    for gg in genes:
        ex_list = gg.get_exon_list()
        for ex in ex_list:
            lngth = ex.get_length()
            if ex.is_annotated():
                annotated_list.append(lngth)
            else:
                denovo_list.append(lngth)

    majiq_io.dump_bin_file([annotated_list, denovo_list], output)




