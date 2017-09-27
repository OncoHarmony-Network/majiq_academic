from majiq.src.config import Config
from majiq.src.constants import *
from voila import constants as voila_const
from voila.api import SpliceGraphs
from voila.splice_graphics import JunctionGraphic, ExonGraphic, GeneGraphic


def init_splicegraph(filename):
    majiq_config = Config()
    # erase splice graph file
    with SpliceGraphs(filename, 'w') as sg:
        sg.add_experiment_names(majiq_config.exp_list)


def gene_to_splicegraph(dict_of_genes, dict_junctions, exon_dict, majiq_config, lock):
    junc_list = []
    junc_l = []
    for gne_id, gne in dict_of_genes.items():
        alt_empty_starts = []
        alt_empty_ends = []
        for jj in dict_junctions[gne_id].values():
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
                    jtype = voila_const.JUNCTION_TYPE_DB_OTHER_RNASEQ if sum(num_reads) > 0 else voila_const.JUNCTION_TYPE_DB
                elif jj.annot and num_reads[exp_idx] > 0:
                    jtype = voila_const.JUNCTION_TYPE_DB_RNASEQ
                else:
                    jtype = voila_const.JUNCTION_TYPE_RNASEQ
                jtype_list.append(jtype)

            #TODO: add transcripts
            junc_l.append((jj.start, jj.end))
            junc_list.append(JunctionGraphic(start=jj.start, end=jj.end, junction_type_list=jtype_list,
                                             reads_list=num_reads, transcripts=[],
                                             intron_retention=jj.intronic))

        exon_list = []
        for ex in exon_dict[gne_id]:
            covered = False
            a3 = []
            alt_start = []
            for jji in set(ex.ib):
                covered = covered or (jji.nreads > 0)
                if jji.end in alt_empty_starts:
                    alt_start.append(jji.end)
                for jidx, jjl in enumerate(junc_l):
                    if jji.end == jjl[1]:
                        a3.append(jidx)

            a5 = []
            alt_ends = []

            for jjo in set(ex.ob):
                covered = covered or (jjo.nreads > 0)
                if jjo.start in alt_empty_starts:
                    alt_ends.append(jjo.start)
                for jidx, jjl in enumerate(junc_l):
                    if jjo.start == jjl[0]:
                        a5.append(jidx)

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
            eg = ExonGraphic(a3, a5, start=ex.start, end=ex.end, exon_type_list=visual_types, coords_extra=extra_coords,
                             intron_retention=False, alt_starts=alt_start, alt_ends=alt_ends)
            exon_list.append(eg)

        ggraph = GeneGraphic(gene_id=gne_id, name=gne['name'], strand=gne['strand'], exons=exon_list,
                             junctions=junc_list, chromosome=gne['chromosome'])

    # lock.acquire()
    with SpliceGraphs(get_builder_splicegraph_filename(majiq_config.outDir), 'r+') as sg:
        sg.add_gene(ggraph)
    # lock.release()
