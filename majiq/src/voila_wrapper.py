from majiq.src.config import Config, get_builder_splicegraph_filename
from majiq.src.constants import FIRST_LAST_JUNC
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
    # erase splice graph file
    with SpliceGraph(filename, delete=True) as sg:
        sg.add_experiment_names(majiq_config.exp_list)


# def gene_to_splicegraph(dict_of_genes, dict_junctions, exon_dict, list_introns, majiq_config, lock):
def gene_to_splicegraph(gne_id, gne, dict_junctions, exon_dict, list_introns, majiq_config, lock):
    # for gne_id, gne in dict_of_genes.items():
    junc_list = []
    junc_l = {}
    alt_empty_starts = []
    alt_empty_ends = []
    jidx = 0

    with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir), 'a') as sg:

        sg.gene(gne_id).add(
            name=gne['name'],
            strand=gne['strand'],
            chromosome=gne['chromosome']
        )

        for jid in sorted(dict_junctions.keys()):
            jj = dict_junctions[jid]

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

        for ex in sorted(exon_dict, key=lambda x: (x.start, x.end)):
            if ex.intron:
                continue

            covered = False
            alt_start = []
            for jj in set(ex.ib):
                covered = covered or (jj.nreads > 0)
                if jj.end in alt_empty_starts:
                    alt_start.append(jj.end)

            alt_ends = []
            for jj in set(ex.ob):
                covered = covered or (jj.nreads > 0)
                if jj.start in alt_empty_ends:
                    alt_ends.append(jj.start)

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
                intron_retention=False,
                annotated=ex.annot,
                alt_starts=alt_start,
                alt_ends=alt_ends
            )

        for info in list_introns:
            sg.exon(gne_id, info.start, info.end).add(
                annotated=info.annot,
                intron_retention=True
            )
