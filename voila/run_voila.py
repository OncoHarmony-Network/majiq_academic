import argparse
import copy
import errno
import fileinput
import json
import os
import sys
import textwrap
import time
from collections import defaultdict, namedtuple, OrderedDict
from os.path import basename

from jinja2 import Environment, FileSystemLoader, escape

import voila.constants as constants
import voila.io_voila as io_voila
import voila.module_locator as module_locator
import voila.utils.utils_voila as utils_voila
from voila.splice_graphics import SpliceGraph
from voila.utils.voilaLog import voilaLog

EXEC_DIR = module_locator.module_path() + "/"


class VoilaNoLSVsException(Exception):
    def __init__(self):
        message = "There are no LSVs detected.  It could be the threshold is too high or the filters are " \
                  "incorrectly set."
        voilaLog().error(message)
        super(VoilaNoLSVsException, self).__init__(message)


def table_marks_set(size):
    """
    Calculate the number of elements to show in LSV tables.

    :param size: total number of LSVs.
    :return: set of total number of elements to show.
    """
    # TODO: Customizable by script options
    ideal_set = (10, 20, 50, 100)
    index = 0
    for _ in ideal_set:
        if size < ideal_set[index]:
            break
        index += 1

    return ideal_set[0:index]


def render_summary(output_dir, output_html, majiq_output, type_summary, threshold=None, extra_args=None):
    """Render a HTML summary using the Jinja2 template system in the output directory specified.

    :param output_dir: output directory for the summaries.
    :param output_html: name for the output html files.
    :param majiq_output: parsed data from old_majiq.
    :param type_summary: type of analysis performed.
    :param threshold: minimum change considered as significant (in deltapsi analysis).
    :param extra_args: additional arguments needed for certain summaries.
    :return: nothing.
    """
    log = voilaLog()

    log.info("Creating the interactive HTML5 summary in %s ..." % output_dir)

    def to_json(value):
        return json.dumps(value.to_dict()).replace('"', '\'')

    def to_json_especial(value):
        return escape(json.dumps(value, cls=utils_voila.LsvGraphicEncoder).replace('\"', '\''))

    def is_combined(spliceg_id):
        return spliceg_id.startswith(constants.COMBINED_PREFIX)

    env = Environment(extensions=["jinja2.ext.do"], loader=FileSystemLoader(os.path.join(EXEC_DIR, "templates/")))
    env.filters.update({'to_json': to_json, 'to_json_especial': to_json_especial, 'debug': utils_voila.debug,
                        'is_combined': is_combined})
    template_file_name = type_summary.replace("-", "_") + "_summary_template.html"
    sum_template = env.get_template(template_file_name)

    if type_summary == constants.ANALYSIS_PSI:
        count_pages = 0
        gene_keys = sorted(majiq_output['genes_dict'].keys())

        log.info("Number of genes detected in Voila: %d." % len(gene_keys))
        log.info("Number of LSVs detected in Voila: %d." % sum(
            [len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
        log.info("Creating HTML5 with splice graphs summaries ...")
        links_dict = {}

        # Subfolder for summary pages
        summaries_subfolder = "%s/%s" % (output_dir, constants.SUMMARIES_SUBFOLDER)
        utils_voila.create_if_not_exists(summaries_subfolder)

        try:
            comb_spliceg_cond1 = \
                [gg for gg in majiq_output['genes_exp'][0].keys() if gg.startswith(constants.COMBINED_PREFIX)][0]
        except IndexError:
            comb_spliceg_cond1 = majiq_output['genes_exp'][0].keys()[0]

        while count_pages * constants.MAX_GENES < len(gene_keys):
            prev_page = None
            next_page = None

            subset_keys = gene_keys[count_pages * constants.MAX_GENES: constants.MAX_GENES * (count_pages + 1)]
            genes_dict = OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)
            log.info("Processing %d out of %d genes ..." % (
                min((count_pages + 1) * constants.MAX_GENES, len(gene_keys)), len(majiq_output['genes_dict'])))
            if (count_pages + 1) * constants.MAX_GENES < len(majiq_output['genes_dict']):
                next_page = str(count_pages + 1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            full_path = "%s/%s" % (summaries_subfolder, name_page)
            voila_output = open(full_path, 'w')
            voila_output.write(
                sum_template.render(tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                                    genes_dict=genes_dict,
                                    prevPage=prev_page,
                                    nextPage=next_page,
                                    namePage=name_page,
                                    lexps=majiq_output['meta_exps'],
                                    genes_exps_list=majiq_output['genes_exp'],
                                    comb_spliceg_cond1=comb_spliceg_cond1
                                    ))
            voila_output.close()
            for g_key, glsv_list in genes_dict.iteritems():
                links_dict[glsv_list[0].name] = "%s/%s" % (constants.SUMMARIES_SUBFOLDER, name_page)
            count_pages += 1
        majiq_output['voila_links'] = links_dict

        # Generate index
        log.info("Creating HTML5 index summary ...")
        sum_template = env.get_template("index_single_summary_template.html")
        voila_output = open(output_dir + "index.html", 'w')
        voila_output.write(sum_template.render(lsvList=majiq_output['lsv_list'],
                                               tableMarks=table_marks_set(len(majiq_output['lsv_list'])),
                                               lexps=majiq_output['meta_exps'],
                                               links_dict=links_dict,
                                               maxLsvs=constants.MAX_LSVS_PSI_INDEX
                                               ))
        voila_output.close()

    elif type_summary == constants.ANALYSIS_DELTAPSI:
        count_pages = 0

        gene_keys = sorted(majiq_output['genes_dict'].keys())
        log.info("Number of genes detected in Voila: %d." % len(gene_keys))
        log.info("Number of LSVs detected in Voila: %d." % sum(
            [len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
        log.info("Creating HTML5 with splice graphs summaries ...")
        links_dict = {}

        # Subfolder for summary pages
        summaries_subfolder = "%s/%s" % (output_dir, constants.SUMMARIES_SUBFOLDER)
        utils_voila.create_if_not_exists(summaries_subfolder)

        try:
            comb_spliceg_cond1 = \
                [gg for gg in majiq_output['genes_exp'][0].keys() if gg.startswith(constants.COMBINED_PREFIX)][0]
        except IndexError:
            comb_spliceg_cond1 = majiq_output['genes_exp'][0].keys()[0]
        try:
            comb_spliceg_cond2 = \
                [gg for gg in majiq_output['genes_exp'][1].keys() if gg.startswith(constants.COMBINED_PREFIX)][0]
        except IndexError:
            comb_spliceg_cond2 = majiq_output['genes_exp'][1].keys()[0]

        while count_pages * constants.MAX_GENES < len(gene_keys):
            prev_page = None
            next_page = None

            subset_keys = gene_keys[count_pages * constants.MAX_GENES: constants.MAX_GENES * (count_pages + 1)]
            genes_dict = OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)

            log.info("Processing %d out of %d genes ..." % (
                min((count_pages + 1) * constants.MAX_GENES, len(gene_keys)), len(majiq_output['genes_dict'])))
            if (count_pages + 1) * constants.MAX_GENES < len(majiq_output['genes_dict']):
                next_page = str(count_pages + 1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            full_path = "%s/%s" % (summaries_subfolder, name_page)
            with open(full_path, 'w') as voila_output:
                voila_output.write(
                    sum_template.render(
                        tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                        genes_dict=genes_dict,
                        genes_exps_list=majiq_output['genes_exp'],
                        prevPage=prev_page,
                        nextPage=next_page,
                        namePage=name_page,
                        threshold=threshold,
                        lexps=majiq_output['meta_exps'],
                        comb_spliceg_cond1=comb_spliceg_cond1,
                        comb_spliceg_cond2=comb_spliceg_cond2
                    ))

            for glsv_list in genes_dict.values():
                links_dict[glsv_list[0]['lsv'].name] = "%s/%s" % (constants.SUMMARIES_SUBFOLDER, name_page)

            count_pages += 1

        majiq_output['voila_links'] = links_dict

        # Generate index
        log.info("Creating HTML5 index summary ...")
        sum_template = env.get_template("index_delta_summary_template.html")

        with open(output_dir + "index.html", 'w') as voila_output:
            voila_output.write(sum_template.render(
                lsvList=majiq_output['lsv_list'],
                tableMarks=table_marks_set(len(majiq_output['lsv_list'])),
                threshold=threshold,
                lexps=majiq_output['meta_exps'],
                links_dict=links_dict,
                maxLsvs=constants.MAX_LSVS_DELTAPSI_INDEX
            ))

    elif type_summary == constants.LSV_THUMBNAILS:
        voila_output = open(output_dir + output_html, 'w')
        voila_output.write(sum_template.render(lsvList=majiq_output,
                                               collapsed=int(extra_args['collapsed'])))

    elif type_summary == constants.SPLICE_GRAPHS:
        count_pages = 0
        genes = majiq_output['gene_dicts'][majiq_output['gene_dicts'].keys()[0]].values()
        log.info("Number of genes detected in Voila: %d." % len(genes))
        while count_pages * constants.MAX_GENES < len(genes):
            prev_page = None
            next_page = None

            log.info("Processing %d out of %d genes ..." % (
                min((count_pages + 1) * constants.MAX_GENES, len(genes)), len(genes)))
            if (count_pages + 1) * constants.MAX_GENES < len(genes):
                next_page = str(count_pages + 1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            voila_output = open(output_dir + name_page, 'w')
            voila_output.write(sum_template.render(
                genes=genes[count_pages * constants.MAX_GENES:(count_pages + 1) * constants.MAX_GENES],
                gene_dicts=majiq_output['gene_dicts'],
                first_sample=[sam for sam in majiq_output['gene_dicts'].keys()][0],
                prevPage=prev_page,
                nextPage=next_page,
                namePage=name_page
            ))
            voila_output.close()
            count_pages += 1

    elif type_summary == constants.COND_TABLE:
        voila_output = open(output_dir + output_html, 'w')
        voila_output.write(sum_template.render(lsvs=majiq_output['lsvs'],
                                               sample_names=majiq_output['sample_names'],
                                               tableMarks=table_marks_set(len(majiq_output['lsvs'].keys())),
                                               cond_pair=majiq_output['cond_pair'],
                                               thres=majiq_output['thres']
                                               ))

    else:
        log.error("summary type not recognized %s." % type_summary, exc_info=1)

    log.info("Copying static files from Voila sources ...")
    utils_voila.copyanything(EXEC_DIR + "templates/static", output_dir + "static")
    utils_voila.copyanything(EXEC_DIR + "templates/static", "%s%s/static" % (output_dir, constants.SUMMARIES_SUBFOLDER))

    log.info("HTML5 Summary successfully created in %s." % output_dir)


def combine_gg(gg_comb_dict, gg_new):
    """Combine a set of gene graphs (new) with an already existing collection."""
    gg = gg_comb_dict[gg_new.gene_id]
    log = voilaLog()
    if gg is None:
        gg_comb_dict[gg_new.gene_id] = copy.deepcopy(gg_new)
        return

    for i, eg in enumerate(gg.exons):
        eg.type_exon = min(eg.type_exon, gg_new.exons[i].type_exon)

    for j, jg in enumerate(gg.junctions):
        jg.type_junction = min(jg.type_junction, gg_new.junctions[j].type_junction)
        jg.num_reads += gg_new.junctions[j].num_reads


def parse_gene_graphics(splice_graph_file, metainfo, gene_names_list=(), gene_ids_list=()):
    """
    Load and combine splice graph files.

    :param splicegraph_flist: list of splice graph files or directory containing splice graphs.
    :param gene_names_list: list of genes of interest.
    :param condition_names: ids for condition 1 [and condition 2, in deltapsi].
    :return: list of genes graphic per condition.
    """
    if not splice_graph_file:
        return

    log = voilaLog()
    genes_exp1_exp2 = []
    gene_names_list = [str(x).lower() for x in gene_names_list]

    log.info("Parsing splice graph information files ...")

    # Combined SpliceGraph data structures
    # gg_combined = defaultdict(lambda: None)
    # gg_combined_name = "%s%s" % (constants.COMBINED_PREFIX, condition_names[grp_i])

    with SpliceGraph(splice_graph_file, 'r') as sg:
        genes_g = sg.get_genes_list()
        gene_experiments_list = sg.get_experiments_list()

    for experiments in [metainfo['experiments1'], metainfo.get('experiments2', [])]:
        genes_exp = defaultdict()
        for experiment in experiments:

            genes_graphic = {}
            genes_g.sort()

            for gene_obj in genes_g:

                if gene_obj.name.lower() in gene_names_list or gene_obj.gene_id in gene_ids_list:
                    experiment_index = gene_experiments_list.index(experiment)

                    genes_graphic[gene_obj.gene_id] = {
                        'json': json.dumps(gene_obj.get_experiment(experiment_index)).replace("\"", "'"),
                        'coordinates': gene_obj.get_coordinates()
                    }

            ggenes_set = set(genes_graphic.keys())

            if not len(ggenes_set):
                log.warning("No gene matching the splice graph file %s." % splice_graph_file)

            if gene_names_list and len(gene_names_list) != len(ggenes_set):
                log.warning(
                    "Different number of genes in splicegraph (%d) and majiq (%d) files! Hint: Are you sure "
                    "you are using bins and splicegraph files from the same execution?" % (
                        len(ggenes_set), len(gene_names_list)))

            # genes_exp[os.path.basename(splice_graph_f)] = genes_graphic
            genes_exp[experiment] = genes_graphic

            # Add combined SpliceGraph (when more than one sample)
            # if len(genes_exp.keys()) > 1:
            #     for gkey, gg_comb in gg_combined.iteritems():
            #         gg_combined[gkey] = [
            #             json.dumps(gg_comb, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'"),
            #             gg_comb.get_strand(),
            #             gg_comb.get_coordinates(),
            #             gg_comb.get_chrom(),
            #             gg_comb.get_name()
            #         ]

        genes_exp1_exp2.append(OrderedDict(sorted(genes_exp.items(), key=lambda t: t[0])))

    log.info("Splice graph information files correctly loaded.")

    return genes_exp1_exp2


def parse_gene_graphics_obj(splicegraph_flist, gene_name_list, condition_names=('group1', 'group2')):
    """
    Load and combine splice graph files. Returns GeneGraphic objects.

    :param splicegraph_flist: list of splice graph files or directory containing splice graphs.
    :param gene_name_list: list of genes of interest.
    :param condition_names: ids for condition 1 [and condition 2, in deltapsi].
    :return: list of genes graphic per condition.
    """
    genes_exp1_exp2 = []
    log = voilaLog()
    log.info("Parsing splice graph information files ...")
    for grp_i, gene_flist in enumerate(splicegraph_flist):
        genes_exp = {}
        splice_files = utils_voila.list_files_or_dir(gene_flist, suffix=constants.SUFFIX_SPLICEGRAPH)

        # Check that the folders have splicegraphs
        if not len(splice_files):
            log.error("No file with extension .%s found in %s." % (constants.SUFFIX_SPLICEGRAPH, gene_flist))
            sys.exit(1)

        for splice_graph_f in splice_files:
            with SpliceGraph(splice_graph_f, 'r') as sg:
                genes = sg.get_genes_list()
                experiments = sg.get_experiments_list()
            genes.sort()

            for experiment_number, experiment_name in enumerate(experiments):
                for gene in genes:
                    if not gene_name_list or gene.gene_id in gene_name_list or gene.name.upper() in gene_name_list:
                        experiment_dict = gene.get_experiment(experiment_number)
                        drop_down_name = '{0}: {1}'.format(basename(splice_graph_f), experiment_name)
                        try:
                            genes_exp[drop_down_name][gene.gene_id] = experiment_dict
                        except KeyError:
                            genes_exp[drop_down_name] = {gene.gene_id: experiment_dict}

        genes_exp1_exp2.append(OrderedDict(sorted(genes_exp.items(), key=lambda t: t[0])))
    log.info("Splice graph information files correctly loaded.")
    return genes_exp1_exp2


def parse_input(args):
    """This method generates an html summary from a old_majiq output file and the rest of the arguments."""

    log = voilaLog()

    type_summary = args.type_analysis
    output_dir = args.output_dir

    if not output_dir.endswith('/'):
        output_dir += '/'
    utils_voila.create_if_not_exists(output_dir)

    log.info("Processing %s summary." % type_summary)

    threshold = None
    pairwise = None
    voila_file = None

    if type_summary != constants.COND_TABLE:
        voila_file = args.majiq_quantifier
        output_html = os.path.splitext(os.path.split(voila_file)[1])[0] + "_" + type_summary.replace("-", "_") + '.html'

    majiq_output = None
    meta_postprocess = {}

    if type_summary == constants.ANALYSIS_PSI:

        gene_names_list = args.gene_names

        voila_input = io_voila.VoilaInput.from_hdf5_file(voila_file,
                                                         lsv_types=args.lsv_types,
                                                         lsv_ids=args.lsv_ids,
                                                         gene_names=gene_names_list)

        majiq_output = utils_voila.lsvs_to_gene_dict(voila_input)

        gene_ids_list = majiq_output['genes_dict'].keys()

        if not len(voila_input.lsvs):
            raise VoilaNoLSVsException()

        # Get gene info
        majiq_output['genes_exp'] = parse_gene_graphics(args.splice_graph,
                                                        voila_input.metainfo,
                                                        gene_names_list=gene_names_list,
                                                        gene_ids_list=gene_ids_list)

        majiq_output['lsv_list'] = [ll for g in majiq_output['genes_dict'].viewvalues() for ll in g]

    if type_summary == constants.ANALYSIS_DELTAPSI:
        threshold = args.threshold
        voila_input = io_voila.VoilaInput.from_hdf5_file(voila_file,
                                                         lsv_ids=args.lsv_ids,
                                                         gene_names=args.gene_names,
                                                         lsv_types=args.lsv_types)

        majiq_output = utils_voila.lsvs_to_gene_dict(voila_input,
                                                     threshold=args.threshold,
                                                     show_all=args.show_all
                                                     )

        gene_ids_list = majiq_output['genes_dict'].keys()

        if not len(voila_input.lsvs):
            raise VoilaNoLSVsException()

        # Get gene info
        majiq_output['genes_exp'] = parse_gene_graphics(
            args.splice_graph,
            voila_input.metainfo,
            gene_ids_list=gene_ids_list,
            gene_names_list=args.gene_names
        )

        majiq_output['lsv_list'] = [ll['lsv'] for g in majiq_output['genes_dict'].viewvalues() for ll in g]

    if type_summary == constants.LSV_THUMBNAILS:
        try:
            majiq_output = []
            with open(voila_file, 'r') as types_file:
                for line in types_file:
                    majiq_output.append(line.rstrip())
        except IOError, e:
            log.error(e.message, exc_info=1)

        meta_postprocess['collapsed'] = args.collapsed

        render_summary(output_dir, output_html, majiq_output, type_summary, threshold, meta_postprocess)

        return

    if type_summary == constants.SPLICE_GRAPHS:
        gene_names_list = None
        if args.gene_names:
            gene_names_list = []
            for gene_name in fileinput.input(args.gene_names):
                gene_names_list.append(gene_name.rstrip().upper())

        log.info("Loading %s." % voila_file)
        majiq_output = {
            'gene_dicts': parse_gene_graphics_obj(
                [[voila_file]],
                gene_names_list,
                condition_names=['samples']
            )[0]
        }
        render_summary(output_dir, output_html, majiq_output, type_summary)
        return

    if type_summary == constants.COND_TABLE:
        cond_pair = args.cond_pair
        sample_files = args.sample_files
        sample_names = args.sample_names
        thres_change = args.thres_change

        gene_names_list = None
        if args.gene_names:
            gene_names_list = []
            for gene_name in fileinput.input(args.gene_names):
                gene_names_list.append(gene_name.rstrip().upper())

        lsv_name_list = None
        if args.lsv_names:
            lsv_name_list = []
            for lsv_name in fileinput.input(args.lsv_names):
                lsv_name_list.append(lsv_name.rstrip().upper())
        print thres_change
        output_html = "%s_%s_comp_table_%.2f.html" % (cond_pair[0], cond_pair[1], thres_change)
        lsvs_dict = io_voila.load_dpsi_tab(sample_files, sample_names, thres_change=thres_change,
                                           filter_genes=gene_names_list, filter_lsvs=lsv_name_list,
                                           pairwise_dir=args.pair_dir, outdir=args.output_dir)
        log.info("LSVs added to the table: %d" % len(lsvs_dict.keys()))
        majiq_output = {'lsvs': lsvs_dict, 'sample_names': sample_names, 'cond_pair': cond_pair, 'thres': thres_change}
        render_summary(output_dir, output_html, majiq_output, type_summary)

    input_parsed = namedtuple('InputParsed',
                              'output_dir output_html majiq_output type_summary threshold meta_postprocess pairwise_dir')
    return input_parsed(output_dir=output_dir, output_html=output_html, majiq_output=majiq_output,
                        type_summary=type_summary, threshold=threshold, meta_postprocess=meta_postprocess,
                        pairwise_dir=pairwise)


def check_list_file(value):
    try:
        with open(value, 'r') as f:
            return [line.strip() for line in f]
    except IOError as e:
        print e.errno
        if e.errno == errno.ENOENT:
            raise argparse.ArgumentTypeError('cannot find file "{0}"'.format(value))
        else:
            raise


def psi_args():
    parser_single = argparse.ArgumentParser(add_help=False)
    parser_single.add_argument('--splice-graph',
                               help='path to splice graph file')
    return parser_single


def deltapsi_args():
    parser_delta = argparse.ArgumentParser(add_help=False)
    parser_delta.add_argument('--splice-graph', type=str,
                              help='Experiment 1 splice graph information file(s) or directory.')
    # Probability threshold used to sum the accumulative probability of inclusion/exclusion.
    parser_delta.add_argument('--threshold', type=float, default=0.2,
                              help='Filter out LSVs with no junction predicted to change over a certain value (in '
                                   'percentage).')
    parser_delta.add_argument('--show-all', dest='show_all', action='store_true', default=False,
                              help='Show all LSVs including those with no junction with significant change predicted')
    parser_delta.add_argument('--pairwise-dir', type=str, dest='pairwise',
                              help='Directory with the pairwise comparisons.')
    return parser_delta


def splicegraphs_args():
    parser_splicegraphs = argparse.ArgumentParser(add_help=False)
    parser_splicegraphs.add_argument('--max',
                                     type=int,
                                     default=20,
                                     help='Maximum number of splice graphs to show (*.splicegraph files may be quite '
                                          'large).')
    parser_splicegraphs.add_argument('--filter-genes',
                                     type=str,
                                     dest='gene_names',
                                     help='File with gene names to filter the results (one gene per line). Use - to '
                                          'type in the gene names.')
    return parser_splicegraphs


def thumbnails_args():
    parser_thumbs = argparse.ArgumentParser(add_help=False)
    parser_thumbs.add_argument('--collapsed', action='store_true', default=False,
                               help='Collapsed LSVs thumbnails in the HTML summary.')
    return parser_thumbs


def conditional_table_args():
    parser_conditional_table = argparse.ArgumentParser(add_help=False)
    parser_conditional_table.add_argument('--cond-pair', dest='cond_pair', required=True, nargs=2, metavar='M1 M2',
                                          help='Condition pair to compare.')

    parser_conditional_table.add_argument('--sample-files', dest='sample_files', required=True, nargs='+',
                                          metavar='M1_M2_sample1 [M1_M2_sample2 ...]',
                                          help='Samples Voila output files.')

    parser_conditional_table.add_argument('--sample-names', dest='sample_names', required=True, nargs='+',
                                          metavar='sample1 [sample2 ...]', help='Sample names.')

    parser_conditional_table.add_argument('--pairwise-dir', dest='pair_dir', required=True, type=str,
                                          help='Root directory where the pairwise delta psi VOILA summaries were created.')

    parser_conditional_table.add_argument('--thres-change', dest='thres_change', type=float, metavar='0.2',
                                          help='Threshold used to filter non-changing LSVs.')

    parser_conditional_table.add_argument('--filter-genes', type=str, dest='gene_names',
                                          help='File with gene names to filter the results (one gene per line). Use - to type '
                                               'in the gene names.')

    parser_conditional_table.add_argument('--filter-lsvs', type=str, dest='lsv_names',
                                          help='File with lsv names to filter the results (one gene per line). Use - to type '
                                               'in the gene names.')
    return parser_conditional_table


def base_args():
    parser_base = argparse.ArgumentParser(add_help=False)
    parser_required = parser_base.add_argument_group('required arguments')

    parser_required.add_argument('-o', '--output',
                                 dest='output_dir',
                                 required=True,
                                 help='path for output directory')

    parser_base.add_argument('-l', '--logger',
                             default=None,
                             help='path for log files')

    parser_base.add_argument('-s', '--silent',
                             action='store_true',
                             help='do not write logs to standard out')

    parser_base.add_argument('--no-html',
                             action='store_true',
                             help='do not write html files')

    parser_base.add_argument('--no-tsv',
                             action='store_true',
                             help='do not generate tab-separated values output file')

    parser_base.add_argument('--gtf',
                             action='store_true',
                             help='generate GTF (GFF2) files for LSVs')

    parser_base.add_argument('--gff',
                             action='store_true',
                             help='generate GFF3 files for LSVs')

    return parser_base


def common_ags():
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('majiq_quantifier',
                               help='location of majiq\'s quantifier file, sometimes referred to as the "voila" file')

    common_parser.add_argument('--lsv-types-file',
                               type=check_list_file,
                               dest='lsv_types',
                               help='location of file that contains a list of LSV types which should remain in the '
                                    'results, one type per line')

    common_parser.add_argument('--lsv-types',
                               nargs='*',
                               default=[],
                               help='LSV types which should remain in the results')

    common_parser.add_argument('--gene-names-file',
                               dest='gene_names',
                               type=check_list_file,
                               default=[],
                               help='location of file that contains a list of gene names which should remain in the '
                                    'results, one name per line')

    common_parser.add_argument('--gene-names',
                               nargs='*',
                               default=[],
                               help='gene names which should remain in the results')

    common_parser.add_argument('--lsv-ids-file',
                               type=check_list_file,
                               dest='lsv_ids',
                               help='location of file that contains a list of LSV IDs which should remain in the '
                                    'results, one ID per line')

    common_parser.add_argument('--lsv-ids',
                               nargs='*',
                               default=[],
                               help='LSV IDs which should remain in the results')

    return common_parser


def main():
    # Time execution time
    start_time = time.time()

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''
            VOILA is a visualization package for Alternative Local Splicing Events.
            -----------------------------------------------------------------------

            '''))
    parser.add_argument('-v', action='version', version=constants.VERSION)

    # Base script options
    base_parser = base_args()

    # Common script options
    common_parser = common_ags()

    # Subparser module to agglutinate all subparsers
    subparsers = parser.add_subparsers(dest='type_analysis')
    subparsers.required = True

    # Single LSV by Gene(s) of interest
    parser_single = psi_args()
    subparsers.add_parser(constants.ANALYSIS_PSI,
                          help='Single LSV analysis by gene(s) of interest.',
                          parents=[base_parser, common_parser, parser_single])

    # Delta LSV
    parser_delta = deltapsi_args()
    subparsers.add_parser(constants.ANALYSIS_DELTAPSI,
                          help='Delta LSV analysis by gene(s) of interest.',
                          parents=[base_parser, common_parser, parser_delta])

    # In-group out-group analysis option
    parser_conditional_table = conditional_table_args()
    subparsers.add_parser(constants.COND_TABLE,
                          help='Generate a HTML table with a list of LSVs changing between conditions in multiple '
                               'samples.',
                          parents=[base_parser, parser_conditional_table])

    # Splice graphs generation option
    parser_splicegraphs = splicegraphs_args()
    subparsers.add_parser(constants.SPLICE_GRAPHS,
                          help='Generate only splice graphs.',
                          parents=[base_parser, common_parser, parser_splicegraphs])

    # Thumbnails generation option
    parser_thumbs = thumbnails_args()
    subparsers.add_parser(constants.LSV_THUMBNAILS,
                          help='Generate LSV thumbnails.',
                          parents=[base_parser, common_parser, parser_thumbs])

    # get args
    args = parser.parse_args()

    # set up logging
    log_filename = 'voila.log'
    if args.logger:
        log_filename = os.path.join(args.logger, log_filename)
    log = voilaLog(filename=log_filename, silent=args.silent)

    if args.no_html and not args.gtf and args.no_tsv:
        log.warning('There are no available output formats.')

    # Parse input
    input_parsed = parse_input(args)

    # Generate outputs
    if not args.no_html:
        render_summary(
            input_parsed.output_dir, input_parsed.output_html, input_parsed.majiq_output,
            input_parsed.type_summary, input_parsed.threshold, input_parsed.meta_postprocess,
        )

    if not args.no_tsv:
        io_voila.write_tab_output(input_parsed)

    if args.gtf:
        io_voila.generic_feature_format_txt_files(input_parsed)

    if args.gff:
        io_voila.generic_feature_format_txt_files(input_parsed, out_gff3=True)

    if input_parsed:
        log.info("Voila! Summaries created in: %s" % input_parsed.output_dir)

    # Add elapsed time
    elapsed_str = utils_voila.secs2hms(time.time() - start_time)
    log.info("Execution time: {0}".format(elapsed_str))


if __name__ == '__main__':
    main()
