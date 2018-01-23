import majiq.src.io as majiq_io
from majiq.grimoire.exon import detect_exons, expand_introns

from majiq.src.config import Config
from majiq.src.constants import *
from voila.api import SpliceGraph


def generate_splicegraph(majiq_config, elem_dict, dict_of_genes):
    init_splicegraph(get_builder_splicegraph_filename(majiq_config.outDir))

    for gne_id, gene_obj in dict_of_genes.items():
        list_exons = []
        dict_junctions = {}
        list_introns = []
        majiq_io.from_matrix_to_objects(gne_id, elem_dict[gne_id], dict_junctions, list_exons, list_introns,
                                        default_index=0)

        detect_exons(dict_junctions, list_exons)
        if majiq_config.ir:
            expand_introns(gne_id, list_introns, list_exons, dict_junctions, default_index=0)
        gene_to_splicegraph(gne_id, gene_obj, dict_junctions, list_exons, list_introns, majiq_config)


def update_splicegraph_junction(sg, gene_id, start, end, nreads, exp):
    # jg = sg.junction("%s:%s-%s" % (gene_id, start, end))
    # jg.update_reads(exp, nreads)
    if -2 not in (start, end):
        sg.junction(gene_id, start, end).update_reads(exp, nreads)


def init_splicegraph(filename):
    majiq_config = Config()
    # erase splice graph file
    with SpliceGraph(filename, delete=True) as sg:
        sg.experiment_names = majiq_config.exp_list
        sg.genome = majiq_config.genome


# def gene_to_splicegraph(dict_of_genes, dict_junctions, exon_dict, list_introns, majiq_config, lock):
def gene_to_splicegraph(gne_id, gne, dict_junctions, exon_dict, list_introns, majiq_config):
    # for gne_id, gne in dict_of_genes.items():
    junc_list = []
    junc_l = {}
    alt_empty_starts = []
    alt_empty_ends = []
    jidx = 0

    with SpliceGraph(get_builder_splicegraph_filename(majiq_config.outDir)) as sg:
        sg.gene(gne_id).add(name=gne['name'], strand=gne['strand'], chromosome=gne['chromosome'])

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
            read_list = [0] * majiq_config.num_experiments

            # junc_list.append(
            #     sg.junction('{0}:{1}-{2}'.format(gne_id, jj.start, jj.end), start=jj.start, end=jj.end,
            #                 reads_list=read_list, transcripts=[], annotated=jj.annot, intron_retention=jj.intronic))
            sg.junction(gne_id, jj.start, jj.end).add(annotated=jj.annot, intron_retention=jj.intronic)

            jidx += 1

        exon_list = []
        for ex in sorted(exon_dict, key=lambda x: (x.start, x.end)):
            if ex.intron:
                continue

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

            ex_start = ex.start if ex.start > -1 else ex.end - 10
            ex_end = ex.end if ex.end > -1 else ex.start + 10
            # exon_list.append(
            #     sg.exon('{0}:{1}-{2}'.format(gne_id, ex.start, ex.end),
            #             a3=a3, a5=a5, start=ex_start, end=ex_end, coords_extra=extra_coords, intron_retention=False,
            #             annotated=ex.annot, alt_starts=alt_start, alt_ends=alt_ends)
            # )
            sg.exon(gne_id, ex.start, ex.end).add(coords_extra=extra_coords, intron_retention=False, annotated=ex.annot,
                                                  alt_starts=alt_start, alt_ends=alt_ends)

        for info in list_introns:
            if info.skip:
                continue
            intr_coord = int(info.start) - 1
            a3 = [junc_l[(intr_coord, info.start)]]
            intr_coord = int(info.end) + 1
            a5 = [junc_l[(info.end, intr_coord)]]
            # eg = sg.exon('{0}:{1}-{2}'.format(gne_id, info.start, info.end), annotated=info.annot,
            #              a3=a3, a5=a5, start=info.start, end=info.end, coords_extra=(), intron_retention=True)
            sg.exon(gne_id, info.start, info.end).add(annotated=info.annot, intron_retention=True)

            # exon_list.append(eg)
        # sg.gene(gne_id, name=gne['name'], strand=gne['strand'], exons=exon_list,
        #         junctions=junc_list, chromosome=gne['chromosome'])
