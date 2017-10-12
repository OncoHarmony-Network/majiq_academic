from majiq.src.config import Config
from majiq.src.constants import *
from voila import constants as voila_const
from voila.api import SpliceGraph


def update_splicegraph_junctions(dict_junctions, junc_mtrx, outDir, exp, lock):


# ggraph = GeneGraphic(gene_id=gne_id, name=gne['name'], strand=gne['strand'], exons=exon_list,
#                      junctions=junc_list, chromosome=gne['chromosome'])
#
    jsum = junc_mtrx.sum(axis=1)

    lock.acquire()
    with SpliceGraph(get_builder_splicegraph_filename(outDir), 'r+') as sg:
        for gid, jlist in dict_junctions.items():
            for xx in jlist.values():
                jg = sg.junction("%s:%s-%s" % (gid, xx.start, xx.end))
                cov = 0

                if xx.index > 0:
                    cov = jsum[xx.index]
                jg.update_reads(exp, cov)
    lock.release()


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

                # TODO: add transcripts
                junc_l[(jj.start, jj.end)] = jidx
                read_list = [1] * majiq_config.num_experiments

                junc_list.append(
                    sg.junction('{0}:{1}-{2}'.format(gne_id, jj.start, jj.end), start=jj.start, end=jj.end, reads_list=read_list,
                                transcripts=[], annotated=jj.annot, intron_retention=jj.intronic))
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

                extra_coords = []
                if ex.annot:
                    if ex.start < ex.db_coords[0]:
                        extra_coords.append([ex.start, ex.db_coords[0] - 1])
                    if ex.end > ex.db_coords[1]:
                        extra_coords.append([ex.db_coords[1] + 1, ex.end])

                # ex_start = ex.start if ex.start > -1 else ex.end-10
                # ex_end = ex.end if ex.end > -1 else ex.start + 10
                exon_list.append(
                    sg.exon('{0}:{1}-{2}'.format(gne_id, ex.start, ex.end),
                            a3=a3, a5=a5, start=ex.start, end=ex.end, coords_extra=extra_coords, intron_retention=False,
                            annotated=ex.annot, alt_starts=alt_start, alt_ends=alt_ends)
                )

            for info in list_introns[gne_id]:
                intr_coord = int(info.start)-1
                junc_list.append(sg.junction('%s:%s-%s' % (gne_id, intr_coord, info.start), start=intr_coord,
                                             end=info.start, transcripts=[], intron_retention=True))

                a3 = [len(junc_list)-1]
                intr_coord = int(info.end)+1
                junc_list.append(sg.junction('%s:%s-%s' % (gne_id, info.end, intr_coord), start=info.end, end=intr_coord,
                                 transcripts=[], intron_retention=True))

                a5 = [len(junc_list)-1]
                eg = sg.exon('{0}:{1}-{2}'.format(gne_id, info.start, info.end),
                             a3=a3, a5=a5, start=info.start, end=info.end, coords_extra=(), intron_retention=True)

                exon_list.append(eg)

            # todo: review this code edit.
            # ggraph = GeneGraphic(gene_id=gne_id, name=gne['name'], strand=gne['strand'], exons=exon_list,
            #                      junctions=junc_list, chromosome=gne['chromosome'])
            #
            # # lock.acquire()
            # with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'r+') as sg:
            #     sg.add_gene(ggraph)
            # # lock.release()
            sg.gene(gne_id, name=gne['name'], strand=gne['strand'], exons=exon_list,
                        junctions=junc_list, chromosome=gne['chromosome'])

