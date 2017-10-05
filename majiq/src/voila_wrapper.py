from majiq.src.config import Config
from majiq.src.constants import *
from voila import constants as voila_const
from voila.api import SpliceGraph


def init_splicegraph(filename):
    majiq_config = Config()
    # erase splice graph file
    with SpliceGraph(filename, 'w') as sg:
        sg.add_experiment_names(majiq_config.exp_list)


def gene_to_splicegraph(dict_of_genes, dict_junctions, exon_dict, list_introns, majiq_config, lock):
    for gne_id, gne in dict_of_genes.items():
        junc_list = []
        junc_l = {}
        alt_empty_starts = []
        alt_empty_ends = []
        jidx = 0

        with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'a') as sg:
            for jid in sorted(dict_junctions[gne_id].keys()):
                jj = dict_junctions[gne_id][jid]
                if jj.start == FIRST_LAST_JUNC:
                    alt_empty_starts.append(jj.end)
                    continue
                if jj.end == FIRST_LAST_JUNC:
                    alt_empty_ends.append(jj.start)
                    continue

                jtype_list = []

                num_reads = [jj.nreads] * majiq_config.num_experiments
                for exp_idx in range(majiq_config.num_experiments):
                    if jj.annot and num_reads[exp_idx] == 0:
                        jtype = voila_const.JUNCTION_TYPE_DB_OTHER_RNASEQ if sum(
                            num_reads) > 0 else voila_const.JUNCTION_TYPE_DB
                    elif jj.annot and num_reads[exp_idx] > 0:
                        jtype = voila_const.JUNCTION_TYPE_DB_RNASEQ
                    else:
                        jtype = voila_const.JUNCTION_TYPE_RNASEQ
                    jtype_list.append(jtype)

                # TODO: add transcripts
                junc_l[(jj.start, jj.end)] = jidx

                # todo: review this code edit.
                # junc_list.append(JunctionGraphic(start=jj.start, end=jj.end, junction_type_list=jtype_list,
                #                                  reads_list=num_reads, transcripts=[],
                #                                  intron_retention=jj.intronic))
                junc_list.append(
                    sg.junction('{0}:{1}-{2}'.format(gne_id, jj.start, jj.end), start=jj.start, end=jj.end,
                                junction_type_list=jtype_list,
                                reads_list=num_reads, transcripts=[],
                                intron_retention=jj.intronic))

                jidx += 1

            exon_list = []
            for ex in sorted(exon_dict[gne_id], key=lambda x: (x.start, x.end)):

                covered = False

                a3 = []
                alt_start = []
                for jj in set(ex.ib):
                    covered = covered or (jj.nreads > 0)
                    if jj.end in alt_empty_starts:
                        alt_start.append(jj.end)
                    if jj.start != FIRST_LAST_JUNC:
                        a3.append(junc_l[(jj.start, jj.end)])

                a5 = []
                alt_ends = []
                for jj in set(ex.ob):
                    covered = covered or (jj.nreads > 0)
                    if jj.start in alt_empty_ends:
                        alt_ends.append(jj.start)
                    if jj.end != FIRST_LAST_JUNC:
                        a5.append(junc_l[(jj.start, jj.end)])

                visual_types = []
                for exp_idx in range(majiq_config.num_experiments):
                    if ex.start == EMPTY_COORD:
                        visual_type = voila_const.EXON_TYPE_MISSING_START
                    elif ex.end == EMPTY_COORD:
                        visual_type = voila_const.EXON_TYPE_MISSING_END
                    elif ex.annot and not covered:
                        visual_type = voila_const.EXON_TYPE_DB
                    elif ex.annot and covered:
                        visual_type = voila_const.EXON_TYPE_DB_RNASEQ
                    elif not ex.annot and covered:
                        visual_type = voila_const.EXON_TYPE_RNASEQ
                    else:
                        visual_type = voila_const.EXON_TYPE_RNASEQ
                    visual_types.append(visual_type)
                # continue
                extra_coords = []
                if ex.annot:
                    if ex.start < ex.db_coords[0]:
                        extra_coords.append([ex.start, ex.db_coords[0] - 1])
                    if ex.end > ex.db_coords[1]:
                        extra_coords.append([ex.db_coords[1] + 1, ex.end])

                # todo: review this code edit.
                # eg = ExonGraphic(a3, a5, start=ex.start, end=ex.end, exon_type_list=visual_types, coords_extra=extra_coords,
                #                  intron_retention=False, alt_starts=alt_start, alt_ends=alt_ends)
                # exon_list.append(eg)

                exon_list.append(
                    sg.exon('{0}:{1}-{2}'.format(gne_id, ex.start, ex.end), a3=a3, a5=a5, start=ex.start, end=ex.end,
                            exon_type_list=visual_types, coords_extra=extra_coords, intron_retention=False,
                            alt_starts=alt_start, alt_ends=alt_ends)
                )

            # for info in list_introns.keys():
            #     jtype_list = [voila_const.EXON_TYPE_RNASEQ] * majiq_config.num_experiments
            #     junc_list.append(JunctionGraphic(start=info[3]-1, end=info[3], junction_type_list=jtype_list,
            #                                      reads_list=0, transcripts=[],
            #                                      intron_retention=True))
            #     a3 = [len(junc_list)]
            #     junc_list.append(JunctionGraphic(start=info[4], end=info[4]+1, junction_type_list=jtype_list,
            #                                      reads_list=0, transcripts=[],
            #                                      intron_retention=True))
            #
            #     a5 = [len(junc_list)]
            #     eg = ExonGraphic(a3, a5, start=info[3], end=info[4],
            #                      exon_type_list=voila_const.EXON_TYPE_RNASEQ,
            #                      coords_extra=(),
            #                      intron_retention=True)
            #     exon_list.append(eg)

            # todo: review this code edit.
            # ggraph = GeneGraphic(gene_id=gne_id, name=gne['name'], strand=gne['strand'], exons=exon_list,
            #                      junctions=junc_list, chromosome=gne['chromosome'])
            #
            # # lock.acquire()
            # with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'r+') as sg:
            #     sg.add_gene(ggraph)
            # # lock.release()
            g = sg.gene(gne_id, name=gne['name'], strand=gne['strand'], exons=exon_list,
                        junctions=junc_list, chromosome=gne['chromosome'])

