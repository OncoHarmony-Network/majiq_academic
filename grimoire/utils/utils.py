import matplotlib
matplotlib.use('Agg')
import os
import logging
import numpy as np
import pickle
import sys
from itertools import izip
from scipy.stats.mstats import mquantiles
import scipy.sparse
from matplotlib import pyplot
import grimoire.mglobals as mglobals
from grimoire.junction import MajiqJunc
from voila.splice_graphics.exonGraphic import ExonGraphic
from voila.splice_graphics.geneGraphic import GeneGraphic 
from voila.splice_graphics.junctionGraphic import JunctionGraphic 


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


def __gc_factor_ind(val, exp_idx):
    res = 0
    for ii, jj in enumerate(mglobals.gc_bins[exp_idx]):
        if val < jj:
            res = ii
    return res


def prepare_lsv_table(lsv_list, non_as, temp_dir):

    out_temp = list()
    for name, ind_list in mglobals.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):

            majiq_table_as = np.zeros(shape=(len(lsv_list[exp_idx])), dtype=np.dtype('object'))
            majiq_table_nonas = np.zeros(shape=(len(non_as[exp_idx])), dtype=np.dtype('object'))

            for iix, lsv in enumerate(lsv_list[exp_idx]):
                majiq_table_as[iix] = lsv.to_majiqLSV(exp_idx)
            for jix, jn in enumerate(non_as[exp_idx]):
                majiq_table_nonas[jix] = MajiqJunc(jn, exp_idx)

            out_temp.append((majiq_table_as,majiq_table_nonas))
    file_pi = open("%s/majiq.pkl" % temp_dir, 'w+')
    pickle.dump(out_temp, file_pi)
    file_pi.close()


def merge_and_create_majiq_file(chr_list, pref_file):

    """

    :param chr_list:
    :param pref_file:
    """
    if pref_file != '':
        pref_file = '%s.' % pref_file

    all_visual = [list() for xx in xrange(mglobals.num_experiments)]
    as_table = [list() for xx in xrange(mglobals.num_experiments)]
    nonas_table = [list() for xx in xrange(mglobals.num_experiments)]

    for chrom in chr_list:
        temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
        temp_filename = '%s/splicegraph.pkl' % temp_dir
        if not os.path.exists(temp_filename):
            continue
        temp_file = open(temp_filename, 'rb')
        visual_gene_list = pickle.load(temp_file)

        filename = "%s/majiq.pkl" % temp_dir
        if not os.path.exists(filename):
            continue
        file_pi2 = open(filename, 'rb')
        temp_table = pickle.load(file_pi2)

        for name, ind_list in mglobals.tissue_repl.items():
            for idx, exp_idx in enumerate(ind_list):
                all_visual[exp_idx].extend(visual_gene_list[mglobals.exp_list[exp_idx]])
                as_table[exp_idx].append(temp_table[exp_idx][0])
                nonas_table[exp_idx].append(temp_table[exp_idx][1])

    for name, ind_list in mglobals.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):

            file_pi = open('%s/%s%s.splicegraph' % (mglobals.outDir, pref_file, mglobals.exp_list[exp_idx]), 'w+')
            pickle.dump(all_visual[exp_idx], file_pi)
            file_pi.close()

            info = dict()
            info['experiment'] = mglobals.exp_list[exp_idx]
            info['GC_bins'] = mglobals.gc_bins[exp_idx]
            info['GC_bins_val'] = mglobals.gc_bins_val[exp_idx]
            info['genome'] = mglobals.genome
            info['num_reads'] = mglobals.num_mapped_reads[exp_idx]

            if len(as_table) == 0:
                continue
            at = np.concatenate(as_table[exp_idx])
            for lsv in at:
                lsv.set_gc_factor(exp_idx)
            nat = np.concatenate(nonas_table[exp_idx])
            for jnc in nat:
                jnc.set_gc_factor(exp_idx)

            file_pi = open('%s/%s%s.majiq' % (mglobals.outDir, pref_file, mglobals.exp_list[exp_idx]), 'w+')
            pickle.dump((info, at, nat), file_pi)
            file_pi.close()


def set_exons_gc_content(chrom, exon_list):

    fastadir_path = "%s/%s/" % (mglobals.genome_path, mglobals.genome)

    #print "Loading chromosome... %s"%chrom
    chrom_path = fastadir_path + chrom + ".fa"
    if not os.path.exists(chrom_path):
        return
    chrom_file = open(chrom_path)
    loaded_chrom = []
    for chrom_line in chrom_file:
        if not chrom_line.startswith(">"):
            loaded_chrom.append(chrom_line.strip("\n"))
    loaded_chrom = ''.join(loaded_chrom)

    for exon in exon_list:
        strt, end = exon.get_coordinates()
        sequence = loaded_chrom[strt:end]
        #reverse the sequence if the strand is reverse
        sequence = sequence.lower()
        if exon.get_strand() == "-":
            new_seq = []
            for char in sequence[::-1]:
                if char == 'g':
                    new_seq.append('c')
                elif char == 'c':
                    new_seq.append('g')
                elif char == 'a':
                    new_seq.append('t')
                elif char == 't':
                    new_seq.append('a')
                else:
                    new_seq.append(char)
            sequence = ''.join(new_seq)
        if len(sequence) == 0:
            print "KKKKseq", exon.get_coordinates(), sequence
        exon.set_gc_content(sequence)


def generate_visualization_output(allgenes, temp_dir):
    gene_list = {}
    for name, ind_list in mglobals.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):
            gene_list[mglobals.exp_list[exp_idx]] = []
            for genes_l in allgenes.values():
                for gg in genes_l:
                    junc_list = []
                    junc_l = []
                    alt_empty_ends = []
                    alt_empty_starts = []
                    for jj in gg.get_all_junctions():

                        cc = jj.get_coordinates()
                        if jj.get_donor() is None:
                            alt_empty_ends.append(cc[0])
                            continue
                        if jj.get_acceptor() is None:
                            alt_empty_starts.append(cc[1])
                            continue

                        if jj.is_annotated() and jj.get_read_num(exp_idx) == 0:
                            jtype = 2
                        elif jj.is_annotated() and jj.get_read_num(exp_idx) > 0:
                            jtype = 0
                        elif not jj.is_annotated() and jj.get_read_num(exp_idx) > mglobals.MINREADS:
                            jtype = 1
                        else:
                            jtype = 1
                            continue
                        junc_l.append(jj.get_coordinates())
                        junc_list.append(JunctionGraphic(cc, jtype, jj.get_read_num(exp_idx)))
                    junc_l = np.asarray(junc_l)
                    exon_list = []
                    for ex in gg.get_exon_list():
                        cc = ex.get_coordinates()
                        a3 = []
                        alt_start = []
                        for ss3 in set(ex.ss_3p_list):
                            for jidx, jjl in enumerate(junc_l):
                                if ss3 == jjl[1]:
                                    a3.append(jidx)
                                    break
                                elif ss3 in alt_empty_starts:
                                    alt_start.append(ss3)
                                    break

                        a5 = []
                        alt_ends = []
                        for ss5 in set(ex.ss_5p_list):
                            for jidx, jjl in enumerate(junc_l):
                                if ss5 == jjl[0]:
                                    a5.append(jidx)
                                    break
                                elif ss5 in alt_empty_starts:
                                    alt_ends.append(ss5)
                                    break

                        if ex.annotated and ex.coverage[exp_idx].sum() == 0.0:
                            visual_type = 2
                        elif ex.annotated and ex.coverage[exp_idx].sum() > 0.0:
                            visual_type = 0
                        elif not ex.annotated and ex.coverage[exp_idx].sum() > 0.0:
                            visual_type = 1
                        else:
                            visual_type = 1
    #                        continue
                        extra_coords = []
                        if ex.annotated:
                            if ex.start < ex.db_coord[0]:
                                extra_coords.append([ex.start, ex.db_coord[0]-1])
                            if ex.end > ex.db_coord[1]:
                                extra_coords.append([ex.db_coord[1]+1, ex.end])
                        eg = ExonGraphic(a3, a5, cc, type_exon=visual_type, coords_extra=extra_coords,
                                         intron_retention=ex.ir, alt_starts=alt_start, alt_ends=alt_ends)
                        exon_list.append(eg)
                    gene_list[mglobals.exp_list[exp_idx]].append(GeneGraphic(gg.get_id(), gg.get_strand(), exon_list,
                                                                             junc_list, gg.get_chromosome()))

    file_pi = open('%s/splicegraph.pkl' % temp_dir, 'w+')
    pickle.dump(gene_list, file_pi)
    file_pi.close()


def prepare_junctions_gc(junc, exp_idx):

    gc = scipy.sparse.lil_matrix(shape=(mglobals.readLen - 16+1))
    gci = np.zeros(shape=(mglobals.readLen - 16+1))
    for jj in range(mglobals.readLen - 16+1):
        if not junc is None and junc.get_gc_content()[exp_idx, jj] != 0:
            #gci[jj] = __gc_factor_ind(junc.get_gc_content()[exp_idx,jj],exp_idx)
            pass
            gc[jj] = mglobals.gc_factor[exp_idx](junc.get_gc_content()[exp_idx, jj])

    if not junc is None:
        junc.add_gc_content_positions(gc)
    return


def print_junc_matrices(mat, tlb=None, fp=None):
    if fp is None:
        out = open('./junc_matrix.tab', 'a+')
    else:
        out = sys.stdout
    out.write("\n=== BEGIN %s === \n\n" % id)
    N, M = mat.shape
    header = [0]*N
    if not tlb is None:
        out.write("Nan\t")
        for ex, (p1, p2) in tlb.items():
            for nid, n in enumerate(p1):
                out.write("%d\t" % (ex+1))
            for nid, n in enumerate(p2):
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
            if lsv.is_Ssource:
                excoord = jj.get_acceptor().get_coordinates()
            else:
                excoord = jj.get_donor().get_coordinates()
            if excoord[1] > alt_coord[0] and excoord[0] < alt_coord[1]:
                name = "%s#%s" % (lsv.id, jidx)
                pcr_lsv = [lsv.exon.get_pcr_name(), name, score]
                pcr_list.append(pcr_lsv)
                print "PCR", ' '.join(pcr_lsv)
    op = open('%s/pcr.pkl' % out_dir, 'w+')
    pickle.dump(pcr_list, op)
    op.close()


def prepare_gc_content(gene_list, temp_dir):
    gc_pairs = {'GC': [[] for xx in xrange(mglobals.num_experiments)],
                'COV': [[] for xx in xrange(mglobals.num_experiments)]}
    for strand, glist in gene_list.items():
        for gn in glist:
            for ex in gn.get_exon_list():
                gc_val = ex.get_gc_content()
                st, end = ex.get_coordinates()
                if gc_val == 0 or end-st < 30:
                    continue
                for exp_n in xrange(mglobals.num_experiments):
                    cov = ex.get_coverage(exp_n)
                    if cov < 1:
                        continue
                    gc_pairs['GC'][exp_n].append(gc_val)
                    gc_pairs['COV'][exp_n].append(cov)

    file_pi = open('%s/gccontent.temppkl' % temp_dir, 'w+')
    pickle.dump(gc_pairs, file_pi)
    file_pi.close()


def gc_factor_calculation(chr_list, nb):

    local_bins = np.zeros(shape=(mglobals.num_experiments, nb+1), dtype=np.dtype('float'))
    local_meanbins = np.zeros(shape=(mglobals.num_experiments, nb),   dtype=np.dtype('float'))
    local_factor = np.zeros(shape=(mglobals.num_experiments, nb),   dtype=np.dtype('float'))

    gc_pairs = {'GC': [[] for xx in xrange(mglobals.num_experiments)],
                'COV': [[] for xx in xrange(mglobals.num_experiments)]}

    # read local files
    for chrom in chr_list:
        temp_dir = "%s/tmp/%s" % (mglobals.outDir, chrom)
        yfile = '%s/gccontent.temppkl' % temp_dir
        if not os.path.exists(yfile):
            continue
        gc_c = pickle.load(open(yfile, 'rb'))
        for exp_n in xrange(mglobals.num_experiments):
            gc_pairs['GC'][exp_n].extend(gc_c['GC'][exp_n])
            gc_pairs['COV'][exp_n].extend(gc_c['COV'][exp_n])

    #print mglobals.tissue_repl
    for tissue, list_idx in mglobals.tissue_repl.items():
        for exp_n in list_idx:
            count = gc_pairs['COV'][exp_n]
            gc = gc_pairs['GC'][exp_n]

            if len(gc) == 0:
                continue

            count, gc = izip(*sorted(izip(count, gc), key=lambda x: x[1]))

            num_regions = len(count)
            nperbin = num_regions / nb

            quant_median = [0.0]*8
            mean_bins = [0]*nb
            bins = [0]*nb

            for ii in range(nb):
                lb = ii * nperbin
                if ii == nb-1:
                    ub = num_regions
                else:
                    ub = (ii+1) * nperbin

                a = np.asarray(count[lb:ub])
                t = np.asarray(gc[lb:ub])

                try:
                    local_bins[exp_n, ii] = t.min()
                except ValueError:
                    local_bins[exp_n, ii] = 0
                if ii == nb -1:
                    local_bins[exp_n, ii+1] = np.max(t)

                #mean_bins[ii] = np.median(t)
                mean_bins[ii] = np.mean(t)
                bins[ii] = mquantiles(a, prob=np.arange(0.1, 0.9, 0.1))
                print "quantiles", bins[ii]

            for qnt in range(8):
                qnt_bns = np.ndarray(len(bins))
                for idx, bb in enumerate(bins):
                    qnt_bns[idx] = bb[qnt]
                print "BINS", qnt_bns
                #quant_median[qnt]=np.median(qnt_bns)
                quant_median[qnt] = np.mean(qnt_bns)

            #print quant_median
            gc_factor = np.zeros(nb, dtype=np.dtype('float'))
            for ii in range(nb):
                offst = np.zeros(len(quant_median), dtype=np.dtype('float'))
                for idx, xx in enumerate(quant_median):
                    offst[idx] = float(bins[ii][idx]) / float(xx)
                gc_factor[ii] = 1/np.mean(offst)

            #print 'MMMMM', gc_factor
            local_meanbins[exp_n] = mean_bins
            local_factor[exp_n] = gc_factor

    mglobals.set_gc_factors(local_bins, local_factor, local_meanbins)


def plot_gc_content():

    idx = 0
    for tissue, list_idx in mglobals.tissue_repl.items():
        pyplot.figure(idx)
        for exp_n in list_idx:
#            f = interpolate.interp1d(mglobals.gc_means[exp_n], mglobals.gc_bins_vaL[exp_n])
#            print mglobals.gc_means[exp_n]
            mn = mglobals.gc_means[exp_n].min()
            mx = mglobals.gc_means[exp_n].max()
            xx = np.arange(mn, mx, 0.001)
            yy = mglobals.gc_factor[exp_n](xx)
            print "XX", xx
            print "Yy", yy
            pyplot.plot(xx, yy, label=mglobals.exp_list[exp_n])
            pyplot.axis((0.3, 0.7, 0.5, 1.5))
            pyplot.title("Gc factor")
            pyplot.grid()
            pyplot.legend(loc='upper left')
#        pyplot.show()
        pyplot.savefig('%s/gcontent_%s.png' % (mglobals.outDir, tissue))
        idx += 1


def recreate_gene_tlb(gene_list):

    for strand, glist in gene_list.items():
        for gn in glist:
            mglobals.gene_tlb[gn.get_id()] = gn


def clear_gene_tlb():
    mglobals.gene_tlb.clear()


def to_gtf(wfile, seq_name, source, gene, mRNA, start_trans, end_trans, strand, exon_l, frame_l):
    SSCORE = "0"
    # Iterate over each exon
    exonOrCDS_list = []
    for i, exon in enumerate(exon_l):
        exonOrCDS_list.append("\t".join([seq_name, source, "%s", exon[0], exon[1], SSCORE, strand,
                                         str(frame_l[i]), "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mRNA)]))

    if strand == '+':
        first_codon = "\t".join([seq_name, source, "start_codon", start_trans, str(int(start_trans) + 2), SSCORE,
                                 strand, ".", "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mRNA)])
        last_codon = "\t".join([seq_name, source, "stop_codon", str(int(end_trans) + 1), str(int(end_trans) + 3),
                                SSCORE, strand, ".", "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mRNA)])

    else:
        last_codon = "\t".join([seq_name, source, "start_codon", str(int(end_trans) - 2), end_trans, SSCORE,
                                strand, ".", "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mRNA)])
        first_codon = "\t".join([seq_name, source, "stop_codon", str(int(start_trans) - 3), str(int(start_trans) - 1),
                                 SSCORE, strand, ".", "gene_id \"%s\"; transcript_id \"%s\";\n" % (gene, mRNA)])

    wfile.write(first_codon)
    for eCDS in exonOrCDS_list:
        wfile.write(eCDS % "CDS")
        wfile.write(eCDS % "exon")
    wfile.write(last_codon)


def gff2gtf(gff_f, out_f=None):
    """Parse a GFF file created by MAJIQ and create a GTF"""

    mRNA = None
    with file_or_stdout(out_f) as wfile:
        with open(gff_f) as gff:
            for gff_l in gff:
                gff_fields = gff_l.strip().split()
                if gff_fields[2] == 'mRNA':
                    if mRNA:
                        to_gtf(wfile, seq_name, source, gene, mRNA, start_trans, end_trans, strand, exon_l, frame_l)
                    exon_l = []
                    frame_l = []
                    ids = gff_fields[8].split(';Parent=')
                    gene = ids[1].split(';')[0]
                    mRNA = ids[0][5:]
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
        to_gtf(wfile, seq_name, source, gene, mRNA, start_trans, end_trans, strand, exon_l, frame_l)

from contextlib import contextmanager as ctx
@ctx
def file_or_stdout(file_name):
    if file_name is None:
        yield sys.stdout
    else:
        with open(file_name, 'w') as out_file:
            yield out_file
