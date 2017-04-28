import numpy as np

from majiq.src import config_old as majiq_config
from voila import constants as voila_const
from voila.splice_graphics import JunctionGraphic, ExonGraphic, GeneGraphic


def generate_visualization_output(allgenes, temp_dir, out_queue):
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

                    ir_type = voila_const.NONE_IR_TYPE
                    if jj.get_donor().is_intron():
                        ir_type = voila_const.IR_TYPE_END
                    elif jj.get_acceptor().is_intron():
                        ir_type = voila_const.IR_TYPE_START

                    junc_l.append(jj.get_coordinates())
                    junc_list.append(JunctionGraphic(cc, jtype, num_reads, transcripts=jj.get_transcript_list(),
                                                     ir=ir_type))

                junc_l = np.asarray(junc_l)
                exon_list = []
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
                out_queue.put([2, ggraph, exp_idx], block=True)
                gene_list.append(ggraph)


            # filename = '%s/%s.splicegraph.pkl' % (temp_dir, majiq_config.exp_list[exp_idx])
            # majiq_io_base.dump_bin_file(gene_list, filename)