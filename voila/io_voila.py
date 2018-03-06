import csv
import os
from os.path import join

import numpy as np

from voila import vlsv
from voila.api import SpliceGraph
from voila.api.view_matrix import ViewDeltaPsi, ViewPsi
from voila.api.view_splice_graph import ViewGene
from voila.constants import JUNCTION_TYPE_RNASEQ
from voila.utils import utils_voila
from voila.utils.voila_log import voila_log
from voila.view.html import Html

__author__ = 'abarrera'


def load_dpairs(pairwise_dir, majiq_output):
    """
    Load pairwise files from MAJIQ analysis.

    :param str pairwise_dir: directory containing pairwise comparisons produced by MAJIQ.
    :param majiq_output: parsed data from old_majiq.
    :return: list of deltapsi lsvs
    :return: name of condition 1
    :return: name of condition 2
    """
    meta_exps = majiq_output['meta_exps']
    lmajiq_pairs = [[None for i in range(len(meta_exps[1]))] for j in range(len(meta_exps[0]))]

    lsv_names = majiq_output['genes_dict'].keys()

    group1_name = meta_exps[0][0]['group']
    group2_name = meta_exps[1][0]['group']

    for idx1 in range(len(meta_exps[0])):
        for idx2 in range(len(meta_exps[1])):
            pairwise_file = "%s/%s_%d_%s_%d.deltapsi.pickle" % (
                pairwise_dir, group1_name, idx1 + 1, group2_name, idx2 + 1)
            try:
                lmajiq_pairs[idx1][idx2] = utils_voila.get_lsv_delta_exp_data(
                    pairwise_file,
                    show_all=True,
                    gene_name_list=lsv_names
                )
            except IOError:
                pass
    return lmajiq_pairs, group1_name, group2_name


def filter_exons(exons):
    for exon in exons:
        if exon.start == -1:
            yield 'nan', exon.end
        elif exon.end == -1:
            yield exon.start, 'nan'
        else:
            yield exon.start, exon.end


def delta_psi_tab_output(args, voila_links):
    def semicolon_join(value_list):
        return ';'.join(str(x) for x in value_list)

    log = voila_log()
    log.info("Creating Tab-delimited output file")
    output_html = Html.get_output_html(args, args.voila_file)
    tsv_file = join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

    with ViewDeltaPsi(args) as m, SpliceGraph(args.splice_graph) as sg:
        metadata = m.metadata

        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction']
        group1 = metadata['group_names'][0]
        group2 = metadata['group_names'][1]
        fieldnames = fieldnames[:3] + ['E(dPSI) per LSV junction',
                                       'P(|dPSI|>=%.2f) per LSV junction' % args.threshold,
                                       '%s E(PSI)' % group1, '%s E(PSI)' % group2]

        fieldnames += ['LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                       'strand', 'Junctions coords', 'Exons coords', 'IR coords']

        if voila_links:
            fieldnames.append('Voila link')

        with open(tsv_file, 'w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for gene_id in m.view_gene_ids():
                gene = sg.gene(gene_id).get

                for lsv_id in m.view_gene_lsvs(gene_id):
                    log.debug('Write TSV row for {0}'.format(lsv_id))
                    lsv = m.delta_psi(lsv_id)

                    lsv_junctions = tuple(ViewGene(gene).lsv_junctions(lsv))
                    lsv_exons = tuple(ViewGene(gene).lsv_exons(lsv))
                    group_means = tuple(lsv.group_means)
                    excl_incl = tuple(lsv.excl_incl)

                    row = {
                        '#Gene Name': gene.name,
                        'Gene ID': gene_id,
                        'LSV ID': lsv_id,
                        'LSV Type': lsv.lsv_type,
                        'A5SS': lsv.prime5,
                        'A3SS': lsv.prime3,
                        'ES': lsv.exon_skipping,
                        'Num. Junctions': lsv.junction_count,
                        'Num. Exons': lsv.exon_count,
                        'chr': gene.chromosome,
                        'strand': gene.strand,
                        'De Novo Junctions': semicolon_join(
                            int(not junc.annotated) for junc in lsv_junctions
                        ),
                        'Junctions coords': semicolon_join(
                            '{0}-{1}'.format(junc.start, junc.end) for junc in lsv_junctions
                        ),
                        'Exons coords': semicolon_join(
                            '{0}-{1}'.format(start, end) for start, end in filter_exons(lsv_exons)
                        ),
                        'IR coords': semicolon_join(
                            '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                        ),
                        'E(dPSI) per LSV junction': semicolon_join(
                            excl_incl[i][1] - excl_incl[i][0] for i in
                            range(np.size(lsv.bins, 0))
                        ),
                        'P(|dPSI|>=%.2f) per LSV junction' % args.threshold: semicolon_join(
                            vlsv.matrix_area(np.array(bin), args.threshold, collapsed_mat=True).sum() for bin in
                            lsv.bins
                        ),
                        '%s E(PSI)' % group1: semicolon_join(
                            '%.3f' % i for i in group_means[0]
                        ),
                        '%s E(PSI)' % group2: semicolon_join(
                            '%.3f' % i for i in group_means[1]
                        )
                    }

                    if voila_links:
                        summary_path = voila_links[gene_id]
                        if not os.path.isabs(summary_path):
                            summary_path = join(os.getcwd(), args.output, summary_path)
                        row['Voila link'] = "file://{0}".format(summary_path)

                    writer.writerow(row)

    log.info("Delimited output file successfully created in: %s" % tsv_file)


def psi_tab_output(args, voila_links):
    def semicolon_join(value_list):
        return ';'.join(str(x) for x in value_list)

    log = voila_log()
    log.info("Creating Tab-delimited output file")

    output_html = Html.get_output_html(args, args.voila_file)
    tsv_file = join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

    with ViewPsi(args) as m, SpliceGraph(args.splice_graph) as sg:
        fieldnames = ['#Gene Name', 'Gene ID', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction',
                      'LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr',
                      'strand', 'Junctions coords', 'Exons coords', 'IR coords']
        if voila_links:
            fieldnames.append('Voila link')

        with open(tsv_file, 'w') as tsv:
            writer = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            for gene_id in m.view_gene_ids():
                gene = sg.gene(gene_id).get
                view_gene = ViewGene(gene)

                for lsv_id in m.view_gene_lsvs(gene_id):
                    lsv = m.psi(lsv_id)
                    lsv_junctions = tuple(view_gene.lsv_junctions(lsv))
                    lsv_exons = tuple(view_gene.lsv_exons(lsv))

                    row = {
                        '#Gene Name': gene.name,
                        'Gene ID': gene_id,
                        'LSV ID': lsv_id,
                        'LSV Type': lsv.lsv_type,
                        'A5SS': lsv.prime5,
                        'A3SS': lsv.prime3,
                        'ES': lsv.exon_skipping,
                        'Num. Junctions': lsv.junction_count,
                        'Num. Exons': lsv.exon_count,
                        'chr': gene.chromosome,
                        'strand': gene.strand,
                        'De Novo Junctions': semicolon_join(
                            int(not junc.annotated) for junc in lsv_junctions
                        ),
                        'Junctions coords': semicolon_join(
                            '{0}-{1}'.format(junc.start, junc.end) for junc in lsv_junctions
                        ),
                        'Exons coords': semicolon_join(
                            '{0}-{1}'.format(start, end) for start, end in filter_exons(lsv_exons)
                        ),
                        'IR coords': semicolon_join(
                            '{0}-{1}'.format(e.start, e.end) for e in lsv_exons if e.intron_retention
                        ),
                        'E(PSI) per LSV junction': semicolon_join(lsv.means),
                        'Var(E(PSI)) per LSV junction': semicolon_join(lsv.variances)
                    }

                    if voila_links:
                        summary_path = voila_links[gene_id]
                        if not os.path.isabs(summary_path):
                            summary_path = join(os.getcwd(), args.output, summary_path)
                        row['Voila link'] = "file://{0}".format(summary_path)

                    log.debug('Write TSV row for {0}'.format(lsv_id))

                    writer.writerow(row)

    log.info("Delimited output file successfully created in: %s" % tsv_file)


# def generic_feature_format_txt_files(args, out_gff3=False):
#     """
#     Create GFF3 files for each LSV.
#     :param majiq_output: majiq data
#     :param args: parsed input data
#     :param out_gff3: output as a GFF3 file
#     :return: None
#     """
#
#     log = voila_log()
#     output_dir = args.output
#
#     if out_gff3:
#         log.info("Create GFF files for LSVs")
#     else:
#         log.info("Create GTF files for LSVs")
#
#     header = "##gff-version 3"
#
#     odir = join(output_dir, "static/doc/lsvs")
#     utils_voila.create_if_not_exists(odir)
#
#     with Voila(args.voila_file, 'r') as v:
#         for lsv in v.get_voila_lsvs(args):
#             lsv_file_basename = "%s/%s" % (odir, lsv.lsv_id)
#             try:
#                 lsv_gff3_str = lsv.to_gff3()
#                 utils_voila.gff2gtf(lsv_gff3_str.split('\n'), "%s.gtf" % lsv_file_basename)
#
#                 # not accessible from command line
#                 if out_gff3:
#                     gff_file = "%s.gff3" % (lsv_file_basename)
#                     with open(gff_file, 'w') as ofile:
#                         ofile.write(header + "\n")
#                         ofile.write(lsv_gff3_str + "\n")
#
#             except UnboundLocalError as e:
#                 log.warning("problem generating GTF file for %s" % lsv.id)
#                 log.error(e)


def get_lsv_info_fieldnames():
    return ['Gene Name', 'Gene ID', 'LSV ID']


def get_lsv_info(lsv):
    return {
        'Gene Name': lsv.name,
        'Gene ID': lsv.lsv_id.split(':')[0],
        'LSV ID': lsv.lsv_id,
    }


def get_lsv_extra_info_fieldnames():
    return ['LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'De Novo Junctions', 'chr', 'strand',
            'Junctions coords', 'Exons coords', 'Exons Alternative Start', 'Exons Alternative End', 'IR coords']


def get_lsv_extra_info(lsv):
    return {
        'LSV Type': lsv.lsv_type,
        'A5SS': lsv.categories['prime5'],
        'A3SS': lsv.categories['prime3'],
        'ES': lsv.categories['ES'],
        'Num. Junctions': lsv.categories['njuncs'],
        'Num. Exons': lsv.categories['nexons'],
        'chr': lsv.chromosome,
        'strand': lsv.strand,
        'De Novo Junctions': semicolon_join(
            int(junc.junction_type == JUNCTION_TYPE_RNASEQ) for junc in lsv.junctions
        ),
        'Junctions coords': semicolon_join(
            '{0}-{1}'.format(junc.start, junc.end) for junc in lsv.junctions
        ),
        'Exons coords': semicolon_join(
            '{0}-{1}'.format(e.start, e.end) for e in lsv.exons
        ),
        'Exons Alternative Start': semicolon_join(
            '|'.join(str(a) for a in e.alt_starts) for e in lsv.exons if e.alt_starts
        ),
        'Exons Alternative End': semicolon_join(
            '|'.join(str(a) for a in e.alt_ends) for e in lsv.exons if e.alt_ends
        ),
        'IR coords': semicolon_join(
            '{0}-{1}'.format(e.start, e.end) for e in lsv.exons if e.intron_retention
        )
    }


def semicolon_join(value_list):
    return ';'.join(str(x) for x in value_list)


def het_tab_output(args):
    voila_log().info('Creating HET TSV file...')

    output_html = get_output_html(args, args.voila_file)
    tsv_file = join(args.output, output_html.rsplit('.html', 1)[0] + '.tsv')

    with Voila(args.voila_file, 'r') as v:

        metainfo = v.get_metainfo()

        fieldnames = get_lsv_info_fieldnames() + list(metainfo['stat_names']) + get_lsv_extra_info_fieldnames()

        lsv_fieldnames = ['Junction ID'] + list(itertools.chain.from_iterable(metainfo['experiment_names']))

        with open(tsv_file, 'w') as tsvfile:
            tsv = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
            tsv.writeheader()

            for lsv in v.get_voila_lsvs(args):
                row = get_lsv_info(lsv)
                row.update(get_lsv_extra_info(lsv))

                for stat_name, junction_stat in zip(metainfo['stat_names'], lsv.het.junction_stats):
                    row[stat_name] = semicolon_join(junction_stat)

                tsv.writerow(row)

                rows = {}

        for lsv in v.get_voila_lsvs(args):

            for group, experiment_names in zip(lsv.het.groups, metainfo['experiment_names']):
                for experiment_index, experiment_name in enumerate(experiment_names):
                    for junction_index, junction_id in enumerate(lsv.junction_ids()):

                        psi = group.get_psi(experiment_index=experiment_index, junction_index=junction_index)

                        try:
                            rows[junction_id][experiment_name] = psi
                        except KeyError:
                            rows[junction_id] = {experiment_name: psi}

            with open(join(args.output, lsv.lsv_id.replace(':', '_') + '.tsv'), 'w') as tsvfile:
                tsv = csv.DictWriter(tsvfile, fieldnames=lsv_fieldnames, delimiter='\t')
                tsv.writeheader()
                for junction_id, row_dict in rows.items():
                    row = {'Junction ID': junction_id}
                    row.update({column: value for column, value in row_dict.items()})
                    tsv.writerow(row)
