import sys
from collections import defaultdict
import json
import os
import textwrap
import collections as cc
import voila.io_voila as io_voila
import voila.module_locator as module_locator
import voila.utils.utils_voila as utils_voila
import voila.constants as constants

try:
    import cPickle as pkl
except ImportError:
    import pickle as pkl

EXEC_DIR = module_locator.module_path() + "/"
VERSION = '0.8.0.yeolab'


def table_marks_set(size):
    """
    Calculate the number of elements to show in LSV tables.

    :param size: total number of LSVs.
    :return: set of total number of elements to show.
    """
    # TODO: Customizable by script options
    ideal_set = (10, 20, 50, 100)
    index = 0
    for mark in ideal_set:
        if size < ideal_set[index]:
            break
        index += 1

    return ideal_set[0:index]


def render_summary(output_dir, output_html, majiq_output, type_summary, threshold=None, extra_args=None, logger=None):
    """Render a HTML summary using the Jinja2 template system in the output directory specified.

    :param output_dir: output directory for the summaries.
    :param output_html: name for the output html files.
    :param majiq_output: parsed data from majiq.
    :param type_summary: type of analysis performed.
    :param threshold: minimum change considered as significant (in deltapsi analysis).
    :param extra_args: additional arguments needed for certain summaries.
    :param logger: logger instance.
    :return: nothing.
    """
    logger.info("Creating the interactive HTML5 summary in %s ..." % output_dir)
    from jinja2 import Environment, FileSystemLoader, escape

    def to_json(value):
        return escape(json.dumps(value, cls=utils_voila.PickleEncoder))

    def to_json_especial(value):
        return escape(json.dumps(value, cls=utils_voila.LsvGraphicEncoder).replace('\"', '\''))

    def is_combined(spliceg_id):
        return spliceg_id.startswith(constants.COMBINED_PREFIX)


    env = Environment(extensions=["jinja2.ext.do"], loader=FileSystemLoader(os.path.join(EXEC_DIR, "templates/")))
    env.filters.update({'to_json': to_json,'to_json_especial': to_json_especial, 'debug': utils_voila.debug, 'is_combined': is_combined})
    template_file_name = type_summary.replace("-", "_") + "_summary_template.html"
    sum_template = env.get_template(template_file_name)

    if type_summary == constants.ANALYSIS_PSI:
        count_pages = 0
        gene_keys = sorted(majiq_output['genes_dict'].keys())

        logger.info("Number of genes detected in Voila: %d." % len(gene_keys))
        logger.info("Number of LSVs detected in Voila: %d." % sum([len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
        logger.info("Creating HTML5 with splice graphs summaries ...")
        links_dict = {}

        # Subfolder for summary pages
        summaries_subfolder = "%s/%s" % (output_dir, constants.SUMMARIES_SUBFOLDER)
        utils_voila.create_if_not_exists(summaries_subfolder)

        try:
            comb_spliceg_cond1 = [gg for gg in majiq_output['genes_exp'][0].keys() if gg.startswith(constants.COMBINED_PREFIX)][0]
        except IndexError:
            comb_spliceg_cond1 = majiq_output['genes_exp'][0].keys()[0]


        while count_pages*constants.MAX_GENES < len(gene_keys):
            prev_page = None
            next_page = None

            subset_keys = gene_keys[count_pages*constants.MAX_GENES: constants.MAX_GENES*(count_pages+1)]
            genes_dict = cc.OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)
            logger.info("Processing %d out of %d genes ..." % (min((count_pages+1)*constants.MAX_GENES, len(gene_keys)), len(majiq_output['genes_dict'])))
            if (count_pages+1)*constants.MAX_GENES < len(majiq_output['genes_dict']):
                next_page = str(count_pages+1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            full_path = "%s/%s" % (summaries_subfolder, name_page)
            voila_output = open(full_path, 'w')
            voila_output.write(sum_template.render(tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
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
                links_dict[glsv_list[0].get_gene_name()] = "%s/%s" % (constants.SUMMARIES_SUBFOLDER, name_page)
            count_pages += 1
        majiq_output['voila_links'] = links_dict

        # Generate index
        logger.info("Creating HTML5 index summary ...")
        sum_template = env.get_template("index_single_summary_template.html")
        voila_output = open(output_dir+"index.html", 'w')
        voila_output.write(sum_template.render( lsvList=majiq_output['lsv_list'],
                                                tableMarks=table_marks_set(len(majiq_output['lsv_list'])),
                                                lexps=majiq_output['meta_exps'],
                                                links_dict=links_dict,
                                                maxLsvs=constants.MAX_LSVS_PSI_INDEX
        ))
        voila_output.close()


    elif type_summary == constants.ANALYSIS_DELTAPSI:
        count_pages = 0

        gene_keys = sorted(majiq_output['genes_dict'].keys())
        logger.info("Number of genes detected in Voila: %d." % len(gene_keys))
        logger.info("Number of LSVs detected in Voila: %d." % sum([len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
        logger.info("Creating HTML5 with splice graphs summaries ...")
        links_dict = {}

        # Subfolder for summary pages
        summaries_subfolder = "%s/%s" % (output_dir, constants.SUMMARIES_SUBFOLDER)
        utils_voila.create_if_not_exists(summaries_subfolder)

        try:
            comb_spliceg_cond1 = [gg for gg in majiq_output['genes_exp'][0].keys() if gg.startswith(constants.COMBINED_PREFIX)][0]
        except IndexError:
            comb_spliceg_cond1 = majiq_output['genes_exp'][0].keys()[0]
        try:
            comb_spliceg_cond2 = [gg for gg in majiq_output['genes_exp'][1].keys() if gg.startswith(constants.COMBINED_PREFIX)][0]
        except IndexError:
            comb_spliceg_cond2 = majiq_output['genes_exp'][1].keys()[0]

        while count_pages*constants.MAX_GENES < len(gene_keys):
            prev_page = None
            next_page = None

            subset_keys = gene_keys[count_pages*constants.MAX_GENES: constants.MAX_GENES*(count_pages+1)]
            genes_dict = cc.OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)

            logger.info("Processing %d out of %d genes ..." % (min((count_pages+1)*constants.MAX_GENES, len(gene_keys)), len(majiq_output['genes_dict'])))
            if (count_pages+1)*constants.MAX_GENES < len(majiq_output['genes_dict']):
                next_page = str(count_pages+1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            full_path = "%s/%s" % (summaries_subfolder, name_page)
            voila_output = open(full_path, 'w')
            voila_output.write(sum_template.render( tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                                                    genes_dict=genes_dict,
                                                    genes_exps_list=majiq_output['genes_exp'],
                                                    prevPage = prev_page,
                                                    nextPage= next_page,
                                                    namePage= name_page,
                                                    threshold=threshold,
                                                    lexps=majiq_output['meta_exps'],
                                                    comb_spliceg_cond1=comb_spliceg_cond1,
                                                    comb_spliceg_cond2=comb_spliceg_cond2
            ))
            voila_output.close()
            for g_key, glsv_list in genes_dict.iteritems():
                links_dict[glsv_list[0]['lsv'].get_gene_name()] = "%s/%s" % (constants.SUMMARIES_SUBFOLDER, name_page)
            count_pages += 1
        majiq_output['voila_links'] = links_dict

        # Generate index
        logger.info("Creating HTML5 index summary ...")
        sum_template = env.get_template("index_delta_summary_template.html")
        voila_output = open(output_dir+"index.html", 'w')
        voila_output.write(sum_template.render( lsvList=majiq_output['lsv_list'],
                                                tableMarks=table_marks_set(len(majiq_output['lsv_list'])),
                                                threshold=threshold,
                                                lexps=majiq_output['meta_exps'],
                                                links_dict=links_dict,
                                                maxLsvs=constants.MAX_LSVS_DELTAPSI_INDEX
        ))
        voila_output.close()

    elif type_summary == constants.LSV_THUMBNAILS:
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render(lsvList=majiq_output,
                                               collapsed=int(extra_args['collapsed'])))

    elif type_summary == constants.SPLICE_GRAPHS:
        count_pages = 0
        genes = majiq_output['genes']

        logger.info("Number of genes detected in Voila: %d." % len(genes))
        while count_pages*constants.MAX_GENES < len(genes):
            prev_page = None
            next_page = None

            logger.info("Processing %d out of %d genes ..." % (min((count_pages+1)*constants.MAX_GENES, len(genes)), len(genes)))
            if (count_pages+1)*constants.MAX_GENES < len(genes):
                next_page = str(count_pages+1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            voila_output = open(output_dir+name_page, 'w')
            voila_output.write(sum_template.render(genes=genes[count_pages*constants.MAX_GENES:(count_pages+1)*constants.MAX_GENES],
                                                   prevPage=prev_page,
                                                   nextPage=next_page,
                                                   namePage=name_page
            ))
            voila_output.close()
            count_pages += 1

    elif type_summary == constants.COND_TABLE:
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render(lsvs=majiq_output['lsvs'],
                                               sample_names=majiq_output['sample_names'],
                                               tableMarks=table_marks_set(len(majiq_output['lsvs'].keys())),
                                               cond_pair=majiq_output['cond_pair'],
                                               thres=majiq_output['thres']
                           ))

    else:
        logger.error("summary type not recognized %s." % type_summary, exc_info=1)

    logger.info("Copying static files from Voila sources ...")
    utils_voila.copyanything(EXEC_DIR+"templates/static", output_dir+"static")
    utils_voila.copyanything(EXEC_DIR+"templates/static", "%s%s/static" % (output_dir, constants.SUMMARIES_SUBFOLDER))

    logger.info("HTML5 Summary successfully created in %s." % output_dir)


def combine_gg(gg_comb_dict, gg_new):
    gg = gg_comb_dict[gg_new.get_id()]
    if gg is None:
        gg_comb_dict[gg_new.get_id()] = gg_new
        return

    for i, eg in enumerate(gg.get_exons()):
        eg.type_exon = min(eg.type_exon, gg_new.get_exons()[i].type_exon)

    for j, jg in enumerate(gg.get_junctions()):
        jg.type_junction = min(jg.type_junction, gg_new.get_junctions()[j].type_junction)
        jg.num_reads += gg_new.get_junctions()[j].num_reads


def parse_gene_graphics(splicegraph_flist, gene_name_list, condition_names=('group1', 'group2'), logger=None):
    """
    Load and combine splice graph files.

    :param splicegraph_flist: list of splice graph files or directory containing splice graphs.
    :param gene_name_list: list of genes of interest.
    :param condition_names: ids for condition 1 [and condition 2, in deltapsi].
    :param logger: logger instance.
    :return: list of genes graphic per condition.
    """
    genes_exp1_exp2 = []
    logger.info("Parsing splice graph information files ...")
    for grp_i, gene_flist in enumerate(splicegraph_flist):
        genes_exp = defaultdict()
        splice_files = utils_voila.list_files_or_dir(gene_flist, suffix=constants.SUFFIX_SPLICEGRAPH)

        # Check that the folders have splicegraphs
        if not len(splice_files):
            logger.error("No file with extension .%s found in %s." % (constants.SUFFIX_SPLICEGRAPH, gene_flist))

        # Combined SpliceGraph data structures
        gg_combined = defaultdict(lambda: None)
        gg_combined_name = "%s%s" % (constants.COMBINED_PREFIX, condition_names[grp_i])

        for splice_graph_f in splice_files:
            logger.info("Loading %s." % splice_graph_f)
            genesG = pkl.load(open(splice_graph_f, 'r'))
            genes_graphic = defaultdict(list)
            genesG.sort()
            for gene_obj in genesG:
                if gene_obj.get_id() in gene_name_list or gene_obj.get_name().upper() in gene_name_list:
                    genes_graphic[gene_obj.get_id()].append(json.dumps(gene_obj, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'"))
                    genes_graphic[gene_obj.get_id()].append(gene_obj.get_strand())
                    genes_graphic[gene_obj.get_id()].append(gene_obj.get_coords())
                    genes_graphic[gene_obj.get_id()].append(gene_obj.get_chrom())
                    genes_graphic[gene_obj.get_id()].append(gene_obj.get_name())

                    # Combine genes from different Splice Graphs
                    combine_gg(gg_combined, gene_obj)

            ggenes_set = set(genes_graphic.keys())
            if not len(ggenes_set):
                logger.warning("No gene matching the splice graph file %s." % splice_graph_f)

            if len(gene_name_list) != len(ggenes_set):
                logger.warning("Different number of genes in splicegraph (%d) and majiq (%d) files! Hint: Are you sure "
                               "you are using bins and splicegraph files from the same execution?" % (len(ggenes_set), len(gene_name_list)))

            genes_exp[os.path.basename(splice_graph_f)] = genes_graphic

        # Add combined SpliceGraph (when more than one sample)
        if len(genes_exp.keys())>1:
            for gkey, gg_comb in gg_combined.iteritems():
                gg_combined[gkey] = [
                    json.dumps(gg_comb, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'"),
                    gg_comb.get_strand(),
                    gg_comb.get_coords(),
                    gg_comb.get_chrom(),
                    gg_comb.get_name()
                ]
            genes_exp[gg_combined_name] = gg_combined
        genes_exp1_exp2.append(cc.OrderedDict(sorted(genes_exp.items(), key=lambda t: t[0])))

    logger.info("Splice graph information files correctly loaded.")
    return genes_exp1_exp2




def create_gff3_txt_files(output_dir, majiq_output, logger):
    """
    Create GFF3 files for each LSV.

    :param output_dir: output directory for the file.
    :param majiq_output: parsed data from majiq.
    :param logger: logger instance.
    :return: nothing.
    """
    logger.info("Saving LSVs files in gff3 format ...")
    if 'genes_dict' not in majiq_output or len(majiq_output['genes_dict'])<1:
        logger.warning("No gene information provided. Genes files are needed to calculate the gff3 files.")
        return

    header = "##gff-version 3"

    odir = output_dir+"/static/doc/lsvs"
    utils_voila.create_if_not_exists(odir)
    for gkey, gvalue in majiq_output['genes_dict'].iteritems():
        for lsv_dict in gvalue:
            lsv = lsv_dict
            if type(lsv_dict) == dict:
                lsv = lsv_dict['lsv']
            lsv_file_basename = "%s/%s" % (odir, lsv.get_id())
            gff_file = "%s.gff3" % (lsv_file_basename)
            with open(gff_file, 'w') as ofile:
                ofile.write(header+"\n")
                ofile.write(lsv.get_gff3()+"\n")
            try:
                utils_voila.gff2gtf(gff_file, "%s.gtf" % lsv_file_basename)
            except UnboundLocalError, e:
                logger.warning("problem generating GTF file for %s" % lsv.get_id())
                logger.error(e.message)

    logger.info("Files saved in %s" % odir)


def create_summary(args):
    """This method generates an html summary from a majiq output file and the rest of the arguments."""

    type_summary    = args.type_analysis
    output_dir      = args.output_dir

    if not output_dir.endswith('/'):
        output_dir += '/'
    utils_voila.create_if_not_exists(output_dir)

    if args.logger is None:
        args.logger = output_dir
    utils_voila.create_if_not_exists(args.logger)

    logger = utils_voila.get_logger("%svoila.log" % args.logger, silent=args.silent)
    logger.info("Execution line: %s" % repr(args))
    logger.info("Processing %s summary." % type_summary)

    threshold   = None
    pairwise    = None

    try:
        voila_file = args.majiq_bins
        output_html = os.path.splitext(os.path.split(voila_file)[1])[0] + "_" + type_summary.replace("-", "_") + '.html'
    except Exception:
        pass

    majiq_output = None
    meta_postprocess = {}

    if type_summary == constants.ANALYSIS_PSI:

        lsv_types = args.lsv_types

        import fileinput
        gene_name_list = []
        if args.gene_names:
            for gene_name in fileinput.input(args.gene_names):
                gene_name_list.append(gene_name.rstrip().upper())
        else:
            gene_name_list = []
        voila_input = io_voila.load_voila_input(voila_file, logger=logger)
        majiq_output = utils_voila.get_lsv_single_exp_data(voila_input, gene_name_list=gene_name_list, lsv_types=lsv_types, logger=logger)

        if not gene_name_list:
            gene_name_list = majiq_output['genes_dict'].keys()

        if not gene_name_list:
            logger.warning("Number of LSVs detected in Voila: 0.")
            logger.info("End of Voila execution.")
            return

        # Get gene info
        majiq_output['genes_exp'] = parse_gene_graphics([args.genes_files], gene_name_list,
                                                        condition_names=[majiq_output['meta_exps'][0]['group'], None],
                                                        logger=logger)
        majiq_output['lsv_list'] = [ll for g in majiq_output['genes_dict'].viewvalues() for ll in g]

    if type_summary == constants.ANALYSIS_DELTAPSI:
        threshold   = args.threshold
        pairwise    = args.pairwise
        gene_name_list = []

        import fileinput
        if args.gene_names:
            for gene_name in fileinput.input(args.gene_names):
                gene_name_list.append(gene_name.rstrip().upper())

        voila_input = io_voila.load_voila_input(voila_file, logger=logger)
        majiq_output = utils_voila.get_lsv_delta_exp_data(voila_input, args.confidence, args.threshold, args.show_all, gene_name_list=gene_name_list, logger=logger)

        if not gene_name_list:
            gene_name_list = majiq_output['genes_dict'].keys()

        if not gene_name_list:
            logger.warning("Number of LSVs detected in Voila with E(Delta(PSI)) > %.2f: None." % threshold)
            logger.warning("End of Voila execution.")
            return

        # Get gene info
        majiq_output['genes_exp'] = parse_gene_graphics([args.genesf_exp1, args.genesf_exp2], gene_name_list,
                                                        condition_names=[majiq_output['meta_exps'][0][0]['group'],
                                                                majiq_output['meta_exps'][1][0]['group']], logger=logger)
        majiq_output['lsv_list'] = [ll['lsv'] for g in majiq_output['genes_dict'].viewvalues() for ll in g]

    if type_summary == constants.LSV_THUMBNAILS:
        try:
            majiq_output = []
            with open(voila_file, 'r') as types_file:
                for line in types_file:
                    majiq_output.append(line.rstrip())
        except IOError, e:
            logger.error(e.message, exc_info=1)

        meta_postprocess['collapsed'] = args.collapsed
        render_summary(output_dir, output_html, majiq_output, type_summary, threshold, meta_postprocess, logger=logger)
        return

    if type_summary == constants.SPLICE_GRAPHS:
        logger.info("Loading %s." % voila_file)
        genesG = pkl.load(open(voila_file, 'r'))[:args.max]
        majiq_output = {'genes': sorted(genesG, key=lambda t: t.id)}
        render_summary(output_dir, output_html, majiq_output, type_summary, logger=logger)
        return

    if type_summary == constants.COND_TABLE:
        cond_pair = args.cond_pair
        sample_files = args.sample_files
        sample_names = args.sample_names
        thres_change = args.thres_change

        output_html = "%s_%s_comp_table_%.2f.html" % (cond_pair[0], cond_pair[1], thres_change)
        lsvs_dict = io_voila.load_dpsi_tab(sample_files, sample_names, thres_change=thres_change)
        logger.info("LSVs added to the table: %d" % len(lsvs_dict.keys()))
        render_summary(output_dir, output_html, {'lsvs': lsvs_dict, 'sample_names': sample_names, 'cond_pair': cond_pair, 'thres': thres_change}, type_summary, logger=logger)
        return

    render_summary(output_dir, output_html, majiq_output, type_summary, threshold, meta_postprocess, logger=logger)
    io_voila.write_tab_output(output_dir, output_html, majiq_output, type_summary, logger=logger, pairwise_dir=pairwise, threshold=threshold)
    create_gff3_txt_files(output_dir, majiq_output, logger=logger)

    logger.info("Voila! Summaries created in: %s" % output_dir)
    return


def main():

    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''
            VOILA is a visualization package for Alternative Local Splicing Events.
            -----------------------------------------------------------------------

            '''))
    parser.add_argument('-v', action='version', version=VERSION)

    # Base script options
    base_parser = argparse.ArgumentParser(add_help=False)
    base_parser.add_argument('-o', '--output', metavar='output_dir', dest='output_dir', type=str, required=True, help='Output directory where the files will be placed.')
    base_parser.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    base_parser.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')

    # Common script options
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('majiq_bins', metavar='majiq_output.pickle', type=str, help='Pickle file with the bins produced by Majiq.')
    # common_parser.add_argument('-o', '--output', metavar='output_dir', dest='output_dir', type=str, required=True, help='Output directory where the files will be placed.')
    common_parser.add_argument('-c', '--confidence', metavar=0.95, dest='confidence', type=float, default=0.95, help='Percentage of confidence required (by default, 0.95).')
    # common_parser.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    # common_parser.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')

    # Subparser module to agglutinate all subparsers
    subparsers = parser.add_subparsers(dest='type_analysis')
    subparsers.required = True

    # Single LSV by Gene(s) of interest
    parser_single = argparse.ArgumentParser(add_help=False)
    parser_single.add_argument('--genes-exp1', nargs='+', required=True, dest='genes_files', metavar='Hippocampus1.splicegraph [Hippocampus2.splicegraph ...]', type=str, help='Splice graph information file(s) or directory with *.splicegraph file(s).')
    parser_single.add_argument('--lsv-types', nargs='*', default=[], type=str, dest='lsv_types', help='LSV type to filter the results. (If no gene list is provided, this option will display only genes containing LSVs of the specified type).')
    parser_single.add_argument('--gene-names-file', type=str, dest='gene_names', help='File with gene names to filter the results (one gene per line). Use - to type in the gene names.')
    subparsers.add_parser(constants.ANALYSIS_PSI, help='Single LSV analysis by gene(s) of interest.', parents=[base_parser, common_parser, parser_single])

    # Delta LSV
    parser_delta = argparse.ArgumentParser(add_help=False)
    parser_delta.add_argument('--threshold', type=float, default=0.2, help='Filter out LSVs with no junction predicted to change over a certain value (in percentage).')  # Probability threshold used to sum the accumulative probability of inclusion/exclusion.
    parser_delta.add_argument('--show-all', dest='show_all', action='store_true', default=False, help='Show all LSVs including those with no junction with significant change predicted')
    parser_delta.add_argument('--pairwise-dir', type=str, dest='pairwise', help='Directory with the pairwise comparisons.')
    parser_delta.add_argument('--genes-exp1', required=True, nargs='+', dest='genesf_exp1', metavar='Hippocampus1.splicegraph [Hippocampus2.splicegraph ...]', type=str, help='Experiment 1 splice graph information file(s) or directory.')
    parser_delta.add_argument('--genes-exp2', required=True, nargs='+', dest='genesf_exp2', metavar='Liver1.splicegraph [Liver2.splicegraph ...]', type=str, help='Experiment 2 splice graph information file(s) or directory.')
    parser_delta.add_argument('--gene-names-file', type=str, dest='gene_names', help='File with gene names to filter the results (one gene per line). Use - to type in the gene names.')
    parser_delta.add_argument('--lsv-types', nargs='+', default=[], type=str, dest='lsv_types', help='LSV type(s) used to filter the results. (If no gene list is provided, this option will display only genes containing LSVs of the specified type).')
    subparsers.add_parser(constants.ANALYSIS_DELTAPSI, help='Delta LSV analysis by gene(s) of interest.', parents=[base_parser, common_parser, parser_delta])

    # Thumbnails generation option (dev) TODO: Delete??
    parser_thumbs = argparse.ArgumentParser(add_help=False)
    parser_thumbs.add_argument('--collapsed', action='store_true', default=False, help='Collapsed LSVs thumbnails in the HTML summary.')
    subparsers.add_parser(constants.LSV_THUMBNAILS, help='Generate LSV thumbnails [DEBUGING!].', parents=[base_parser, common_parser, parser_thumbs])

    # Splice graphs generation option (dev) TODO: Delete??
    parser_splice_graphs = argparse.ArgumentParser(add_help=False)
    parser_splice_graphs.add_argument('--max', type=int, default=20, help='Maximum number of splice graphs to show (*.splicegraph files may be quite large).')
    subparsers.add_parser(constants.SPLICE_GRAPHS, help='Generate only splice graphs [DEBUGING!].', parents=[base_parser, common_parser, parser_splice_graphs])

    # In-group out-group analysis option (dev) TODO: Delete??
    parser_comptable = argparse.ArgumentParser(add_help=False)
    parser_comptable.add_argument('-cond-pair', dest='cond_pair', required=True, nargs=2, metavar='M1 M2', help='Condition pair to compare.')
    parser_comptable.add_argument('-sample-files', dest='sample_files', required=True, nargs='+', metavar='M1_M2_sample1 [M1_M2_sample2 ...]', help='Samples Voila output files.')
    parser_comptable.add_argument('-sample-names', dest='sample_names', required=True, nargs='+', metavar='sample1 [sample2 ...]', help='Sample names.')
    parser_comptable.add_argument('--thres-change',  dest='thres_change', type=float, metavar='0.2', help='Threshold used to filter non-changing LSVs.')
    subparsers.add_parser(constants.COND_TABLE, help='Generate a HTML table with a list of LSVs changing between conditions in multiple samples [DEBUGING!].', parents=[base_parser, parser_comptable])

    args = parser.parse_args()
    create_summary(args)
    sys.exit(1)

if __name__ == '__main__':
    main()

