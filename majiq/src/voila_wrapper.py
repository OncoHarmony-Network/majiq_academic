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
    if FIRST_LAST_JUNC not in (start, end):
        sg.junction(gene_id, start, end).update_reads(exp, nreads)


def init_splicegraph(filename):
    majiq_config = Config()

    # erase splice graph file
    with SpliceGraph(filename, delete=True) as sg:
        sg.experiment_names = majiq_config.exp_list
        sg.genome = majiq_config.genome


def gene_to_splicegraph(gne_id, gne, dict_junctions, exon_dict, list_introns, majiq_config):
    alt_empty_starts = []
    alt_empty_ends = []

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
            sg.junction(gne_id, jj.start, jj.end).add(annotated=jj.annot, intron_retention=jj.intronic)

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

            exon = sg.exon(gne_id, ex.start, ex.end)
            if not exon:
                exon.add(coords_extra=extra_coords, intron_retention=False, annotated=ex.annot, alt_starts=alt_start,
                         alt_ends=alt_ends)
            else:
                print('duplicate exon:', gne_id, ex.start, ex.end)

        for info in list_introns:
            if info.skip:
                continue

            exon = sg.exon(gne_id, info.start, info.end)
            if not exon:
                sg.exon(gne_id, info.start, info.end).add(annotated=info.annot, intron_retention=True)
            else:
                print('duplicate exon (introns):', gne_id, info.start, info.end)
