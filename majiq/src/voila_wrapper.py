from majiq.src.config import Config
from majiq.src.constants import *
from voila.api import SpliceGraph


def update_splicegraph_junctions(dict_junctions, junc_mtrx, outDir, exp, lock):
    jsum = junc_mtrx.sum(axis=1)

    with SpliceGraph(get_builder_splicegraph_filename(outDir)) as sg:
        for gid, jlist in dict_junctions.items():
            for xx in jlist.values():
                jg = sg.junction(gid, xx.start, xx.end)

                if xx.index > 0:
                    cov = jsum[xx.index]
                    jg.update_reads(cov, exp)


def init_splicegraph(filename):
    majiq_config = Config()
    with SpliceGraph(filename, delete=True) as sg:
        sg.add_experiment_names(majiq_config.exp_list)


def gene_to_splicegraph(dict_of_genes, dict_junctions, exon_dict, list_introns, majiq_config, lock):
    for gne_id, gne in dict_of_genes.items():
        junc_l = {}
        alt_empty_starts = []
        alt_empty_ends = []
        jidx = 0

        with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir)) as sg:
            sg.gene(gne_id).add(
                name=gne['name'],
                strand=gne['strand'],
                chromosome=gne['chromosome']
            )

            for jid in sorted(dict_junctions[gne_id].keys()):
                jj = dict_junctions[gne_id][jid]
                # if jj.intronic: continue
                if jj.start == FIRST_LAST_JUNC:
                    alt_empty_starts.append(jj.end)
                    continue
                if jj.end == FIRST_LAST_JUNC:
                    alt_empty_ends.append(jj.start)
                    continue

                # TODO: add transcripts
                junc_l[(jj.start, jj.end)] = jidx

                sg.junction(gne_id, jj.start, jj.end).add(
                    annotated=jj.annot,
                    intron_retention=jj.intronic
                )

                jidx += 1

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

                extra_coords = []
                if ex.annot:
                    if ex.start < ex.db_coords[0]:
                        extra_coords.append([ex.start, ex.db_coords[0] - 1])
                    if ex.end > ex.db_coords[1]:
                        extra_coords.append([ex.db_coords[1] + 1, ex.end])

                # ex_start = ex.start if ex.start > -1 else ex.end - 10
                # ex_end = ex.end if ex.end > -1 else ex.start + 10
                sg.exon(gne_id, ex.start, ex.end).add(
                    coords_extra=extra_coords,
                    intron_retention=ex.intron,
                    annotated=ex.annot,
                    alt_starts=alt_start,
                    alt_ends=alt_ends
                )

            # for info in list_introns[gne_id]:
            #     # intr_coord = int(info.start) - 1
            #     # junc_list.append(sg.junction('%s:%s-%s' % (gne_id, intr_coord, info.start), start=intr_coord,
            #     #                              end=info.start, transcripts=[],
            #     #                              intron_retention=voila_const.IR_TYPE_START))
            #     #
            #     # a3 = [len(junc_list)-1]
            #     # intr_coord = int(info.end)+1
            #     # junc_list.append(sg.junction('%s:%s-%s' % (gne_id, info.end, intr_coord), start=info.end, end=intr_coord,
            #     #                  transcripts=[], intron_retention=voila_const.IR_TYPE_END))
            #     #
            #     # a5 = [len(junc_list)-1]
            #
            #     # intr_coord = int(info.start) - 1
            #     # a3 = [junc_l[(intr_coord, info.start)]]
            #     # intr_coord = int(info.end) + 1
            #     # a5 = [junc_l[(info.end, intr_coord)]]
            #
            #     eg = sg.exon(gne_id, info.start, info.end)
            #     if not eg:
            #         eg.add(
            #             annotated=info.annot,
            #             intron_retention=True
            #         )
            #     else:
            #         print(dict(eg))

                    # exon_list.append(eg)

                    # todo: review this code edit.
                    # ggraph = GeneGraphic(gene_id=gne_id, name=gne['name'], strand=gne['strand'], exons=exon_list,
                    #                      junctions=junc_list, chromosome=gne['chromosome'])
                    #
                    # # lock.acquire()
                    # with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'r+') as sg:
                    #     sg.add_gene(ggraph)
                    # # lock.release()

                    # sg.gene(gne_id).add(
                    #     name=gne['name'],
                    #     strand=gne['strand'],
                    #     exons=exon_list,
                    #     junctions=junc_list,
                    #     chromosome=gne['chromosome']
                    # )

                    # gene.junctions = junc_list
                    # gene.exons = exon_list
