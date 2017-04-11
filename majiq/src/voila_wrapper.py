import numpy as np

from majiq.src import config as majiq_config
from voila import constants as voila_const
from voila.splice_graphics import JunctionGraphic, ExonGraphic, GeneGraphic, SpliceGraph
from majiq.src.constants import *


def init_splicegraph(filename):
    with SpliceGraph(filename, 'w') as sg:
        sg.erase_splice_graph_file()
        sg.add_experiment_names(majiq_config.exp_list)


def gene_to_splicegraph(gne, lock):

    lock.acquire()
    with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'r+') as sg:

        junc_list = []
        junc_l = []
        alt_empty_ends = []
        alt_empty_starts = []
        for jj in gne.get_all_junctions(filter=False):

            cc = jj.get_coordinates()

            if jj.get_donor() is None:
                alt_empty_ends.append(cc[1])
                continue
            if jj.get_acceptor() is None:
                alt_empty_starts.append(cc[0])
                continue
            num_reads = jj.get_coverage().sum(axis=1)
            jj.get_acceptor().coverage += num_reads.sum()
            jj.get_donor().coverage += num_reads.sum()
            jtype_list = []
            for exp_idx in xrange(majiq_config.num_experiments):
                if jj.annotated and num_reads[exp_idx] == 0:
                    if (num_reads.sum() - num_reads[exp_idx]) > 0:
                        jtype = voila_const.JUNCTION_TYPE_DB_OTHER_RNASEQ
                    else:
                        jtype = voila_const.JUNCTION_TYPE_DB
                elif jj.annotated and num_reads[exp_idx] > 0:
                    jtype = voila_const.JUNCTION_TYPE_DB_RNASEQ
                elif not jj.annotated and num_reads[exp_idx] > majiq_config.MINREADS:
                    jtype = voila_const.JUNCTION_TYPE_RNASEQ
                else:
                    jtype = voila_const.JUNCTION_TYPE_RNASEQ
                    # continue
                jtype_list.append(jtype)

            ir_type = voila_const.NONE_IR_TYPE
            if jj.get_donor().is_intron():
                ir_type = voila_const.IR_TYPE_END
            elif jj.get_acceptor().is_intron():
                ir_type = voila_const.IR_TYPE_START

            junc_l.append(jj.get_coordinates())
            junc_list.append(JunctionGraphic(start=cc[0], end=cc[1], junction_type_list=jtype_list,
                                             reads_list=num_reads, transcripts=jj.get_transcript_list(),
                                             intron_retention=ir_type))

        junc_l = np.asarray(junc_l)
        exon_list = []
        for ex in gne.get_exon_list():
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

            ex_reads = ex.get_coverage()
            visual_types = []
            for exp_idx in xrange(majiq_config.num_experiments):

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
                visual_types.append(visual_type)
            # continue
            extra_coords = []
            if ex.is_annotated():
                if ex.start < ex.db_coord[0]:
                    extra_coords.append([ex.start, ex.db_coord[0] - 1])
                if ex.end > ex.db_coord[1]:
                    extra_coords.append([ex.db_coord[1] + 1, ex.end])
            eg = ExonGraphic(a3, a5, start=cc[0], end=cc[1], exon_type_list=visual_types, coords_extra=extra_coords,
                             intron_retention=ex.get_ir(), alt_starts=alt_start, alt_ends=alt_ends)
            exon_list.append(eg)
        ggraph = GeneGraphic(gene_id=gne.get_id(), name=gne.get_name(), strand=gne.get_strand(), exons=exon_list,
                             junctions=junc_list, chromosome=gne.get_chromosome())

        sg.add_gene(ggraph)

    lock.release()
