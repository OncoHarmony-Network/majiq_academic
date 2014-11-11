import sys
from collections import defaultdict
import json
import os
import textwrap
import collections as cc
import voila.module_locator as module_locator
import voila.utils.utils_voila as utils_voila
import constants
try:
    import cPickle as pkl
except ImportError:
    import pickle as pkl

EXEC_DIR = module_locator.module_path() + "/"
VERSION = '0.1.0'

def table_marks_set(size):
    # TODO: Customizable by script options
    ideal_set = (10, 20, 50, 100)
    index = 0
    for mark in ideal_set:
        if size < ideal_set[index]:
            break
        index += 1

    return ideal_set[0:index]


def debug(text):
    print text
    return ''


def render_summary(output_dir, output_html, majiq_output, type_summary, threshold=None, post_process_info=None, logger=None):
    """
    Rendering the summary template to create the HTML file.

    @param output_dir: directory where the HTML will be created
    @param majiq_output: event list from MAJIQ output used as input in Voila
    @param type_summary: defines summary template used
    """
    logger.info("Creating the interactive HTML5 summary in %s ..." % output_dir)
    from jinja2 import Environment, FileSystemLoader, escape

    def to_json(value):
        return escape(json.dumps(value, cls=utils_voila.PickleEncoder))

    env = Environment(extensions=["jinja2.ext.do"], loader=FileSystemLoader(os.path.join(EXEC_DIR, "templates/")))
    env.filters.update({'to_json': to_json, 'debug': debug})
    sum_template = env.get_template(type_summary.replace("-", "_") + "_summary_template.html")

    if type_summary == constants.ANALYSIS_PSI:
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render(lsvList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               metadata=majiq_output['metadata'],
                                               lexps=majiq_output['meta_exps']
                                               ))
        voila_output.close()
    elif type_summary == constants.ANALYSIS_PSI_GENE:
        # Max. 10 genes per page, create as many HTMLs as needed

        MAX_GENES = 10
        count_pages = 0
        gene_keys = sorted(majiq_output['genes_dict'].keys())

        logger.info("Number of genes detected in Voila: %d." % len(gene_keys))
        logger.info("Number of LSVs detected in Voila: %d." % sum([len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
        while count_pages*MAX_GENES < len(gene_keys):
            prev_page = None
            next_page = None

            subset_keys = gene_keys[count_pages*MAX_GENES: MAX_GENES*(count_pages+1)]
            genes_dict = cc.OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)
            logger.info("Processing %d out of %d genes ..." % (min((count_pages+1)*MAX_GENES, len(gene_keys)), len(majiq_output['genes_dict'])))
            if (count_pages+1)*MAX_GENES < len(majiq_output['genes_dict']):
                next_page = str(count_pages+1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            voila_output = open(output_dir+name_page, 'w')
            voila_output.write(sum_template.render(tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                                                   # genes_json=genes_json_dict,
                                                   genes_dict=genes_dict,
                                                   prevPage = prev_page,
                                                   nextPage= next_page,
                                                   namePage= name_page,
                                                   lexps=majiq_output['meta_exps'],
                                                   genes_exps_list=majiq_output['genes_exp']
            ))
            voila_output.close()
            count_pages += 1

    elif type_summary == constants.ANALYSIS_DELTAPSI_GENE:
        # Max. 10 genes per page, create as many HTMLs as needed

        MAX_GENES = 10
        count_pages = 0

        gene_keys = sorted(majiq_output['genes_dict'].keys())
        logger.info("Number of genes detected in Voila: %d." % len(gene_keys))
        logger.info("Number of LSVs detected in Voila: %d." % sum([len(majiq_output['genes_dict'][g]) for g in majiq_output['genes_dict']]))
        while count_pages*MAX_GENES < len(gene_keys):
            prev_page = None
            next_page = None

            subset_keys = gene_keys[count_pages*MAX_GENES: MAX_GENES*(count_pages+1)]
            genes_dict = cc.OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)

            logger.info("Processing %d out of %d genes ..." % (min((count_pages+1)*MAX_GENES, len(gene_keys)), len(majiq_output['genes_dict'])))
            if (count_pages+1)*MAX_GENES < len(majiq_output['genes_dict']):
                next_page = str(count_pages+1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            voila_output = open(output_dir+name_page, 'w')
            voila_output.write(sum_template.render( tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                                                    genes_dict=genes_dict,
                                                    genes_exps_list=majiq_output['genes_exp'],
                                                    prevPage = prev_page,
                                                    nextPage= next_page,
                                                    namePage= name_page,
                                                    threshold=threshold,
                                                    lexps=majiq_output['meta_exps']
            ))
            voila_output.close()
            count_pages += 1

    elif type_summary == constants.ANALYSIS_DELTAPSI:
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render( lsvList=majiq_output['event_list'],
                                                tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                                metadata=majiq_output['metadata'],
                                                threshold=threshold,
                                                lexps=majiq_output['meta_exps']
        ))
        voila_output.close()

    else:
        logger.error("summary type not recognized %s." % type_summary, exc_info=1)

    logger.info("Copying static files from Voila sources to %sstatic/ ..." % output_dir)
    utils_voila.copyanything(EXEC_DIR+"templates/static", output_dir+"static")

    logger.info("HTML5 Summary successfully created in %s." % output_dir)


class ParseError(Exception):
    def __init__(self, msg, logger=None):
        self.msg = msg
        self.logger = logger

    def __repr__(self):
        return self.msg


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


def parse_gene_graphics(gene_exps_flist, gene_name_list, groups=('group1', 'group2'), logger=None):
    genes_exp1_exp2 = []
    logger.info("Parsing splice graph information files ...")
    for grp_i, gene_flist in enumerate(gene_exps_flist):
        genes_exp = defaultdict()
        splice_files = utils_voila.list_files_or_dir(gene_flist, suffix=constants.SUFFIX_SPLICEGRAPH)

        # Check that the folders have splicegraphs
        if not len(splice_files):
            raise ParseError("No file with extension .%s found in %s." % (constants.SUFFIX_SPLICEGRAPH, gene_flist), logger=logger)

        # Combined SpliceGraph data structures
        gg_combined = defaultdict(lambda: None)
        gg_combined_name = "ALL_%s" % groups[grp_i]

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
                raise ParseError("Different number of genes in splicegraph (%d) and majiq (%d) files." % (len(ggenes_set), len(gene_name_list)), logger=logger)

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


def render_tab_output(output_dir, output_html, majiq_output, type_summary, logger=None):

    ofile_str = "%s%s.%s" % (output_dir, output_html.rsplit('.html', 1)[0], constants.EXTENSION)
    tlb_categx = {'A5SS': 'prime5', 'A3SS': 'prime3', 'Num. Junctions': 'njuncs', 'Num. Exons': 'nexons', 'ES': 'ES'}

    logger.info("Creating Tab-delimited output file in %s..." % ofile_str)
    with open(ofile_str, 'w+') as ofile:
        headers = ['Gene name', 'LSV ID', 'E(PSI) per LSV junction', 'Var(E(PSI)) per LSV junction', 'LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'chr', 'strand', 'Junctions coords', 'Exons coords', 'Exons Alternative Start', 'Exons Alternative End']
        if 'delta' in type_summary:
            headers[2] = 'E(Delta(PSI)) per LSV junction'
            headers[3] = 'P(Delta(PSI)>%s) per LSV junction' % .1

        ofile.write(constants.DELIMITER.join(headers))
        ofile.write('\n')

        for gene in majiq_output['genes_dict']:
            for llsv in majiq_output['genes_dict'][gene]:
                lline = []
                lline.extend([gene, llsv[0].get_id()])
                lexpected = []
                lconfidence = []
                for i, bins in enumerate(llsv[0].get_bins()):
                    if 'delta' in type_summary:
                        lexpected.append(str(sum(llsv[0].get_excl_incl()[i])))
                        lconfidence.append(str(utils_voila.get_prob_delta_psi_greater_v(bins, float(lexpected[-1]), 0.1)))
                    else:
                        lexpected.append(repr(llsv[0].get_means()[i]))
                        lconfidence.append(repr(llsv[0].get_variances()[i]))

                lline.append('; '.join(lexpected))
                lline.append('; '.join(lconfidence))

                lline.append(llsv[0].get_type())
                lline.append(repr(llsv[0].get_categories()[tlb_categx['A5SS']]))
                lline.append(repr(llsv[0].get_categories()[tlb_categx['A3SS']]))
                lline.append(repr(llsv[0].get_categories()[tlb_categx['ES']]))
                lline.append(repr(llsv[0].get_categories()[tlb_categx['Num. Junctions']]))
                lline.append(repr(llsv[0].get_categories()[tlb_categx['Num. Exons']]))

                lline.append(llsv[1][4].get_chrom())
                lline.append(llsv[1][4].get_strand())

                lline.append('; '.join(['-'.join(str(c) for c in junc.get_coords()) for junc in llsv[1][4].get_junctions()]))
                lline.append('; '.join(['-'.join(str(c) for c in exon.get_coords()) for exon in llsv[1][4].get_exons()]))

                try:
                    lline.append('; '.join([' '.join([str(c) for c in exon.get_alt_starts()]) for exon in llsv[1][4].get_exons()]))
                    lline.append('; '.join([' '.join([str(c) for c in exon.get_alt_ends()]) for exon in llsv[1][4].get_exons()]))
                except TypeError:
                    pass

                ofile.write(constants.DELIMITER.join(lline))
                ofile.write('\n')

    logger.info("Delimited output file successfully created in: %s" % ofile_str)


def create_gff3_txt_files(output_dir, majiq_output, logger):
    logger.info("Saving LSVs files in gff3 format ...")
    if 'genes_dict' not in majiq_output or len(majiq_output['genes_dict'])<1:
        logger.warning("No gene information provided. Genes files are needed to calculate the gff3 files.")
        return

    header = "##gff-version 3"

    odir = output_dir+"/static/doc/lsvs"
    utils_voila.create_if_not_exists(odir)
    for gkey, gvalue in majiq_output['genes_dict'].iteritems():
        for lsv in gvalue:
            lsv_file_basename = "%s/%s" % (odir, lsv[0].get_id())
            gff_file = "%s.gff3" % (lsv_file_basename)
            with open(gff_file, 'w') as ofile:
                ofile.write(header+"\n")
                ofile.write(lsv[0].get_gff3()+"\n")
            utils_voila.gff2gtf(gff_file, "%s.gtf" % lsv_file_basename)

    logger.info("Files saved in %s" % odir)

def create_summary(args):
    """This method generates an html summary from a majiq output file"""

    type_summary    = args.type_analysis
    majiq_bins_file = args.majiq_bins
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

    # meta_preprocess = args.meta_preprocess
    threshold       = None

    output_html = os.path.splitext(os.path.split(majiq_bins_file)[1])[0] + "_" + type_summary.replace("-", "_") + '.html'
    majiq_output = None
    meta_postprocess = {}

    if type_summary == constants.ANALYSIS_PSI:
        majiq_output = utils_voila.get_lsv_single_exp_data(majiq_bins_file, args.confidence)

    if type_summary == constants.ANALYSIS_PSI_GENE:

        lsv_types = args.lsv_types

        import fileinput
        gene_name_list = []
        if args.gene_names:
            for gene_name in fileinput.input(args.gene_names):
                gene_name_list.append(gene_name.rstrip().upper())
        else:
            gene_name_list = []
        majiq_output = utils_voila.get_lsv_single_exp_data(majiq_bins_file, args.confidence, gene_name_list=gene_name_list, lsv_types=lsv_types, logger=logger)

        if not gene_name_list:
            gene_name_list = majiq_output['genes_dict'].keys()

        if not gene_name_list:
            logger.warning("Number of LSVs detected in Voila: 0.")
            logger.info("End of Voila execution.")
            return

        # Get gene info
        majiq_output['genes_exp'] = parse_gene_graphics([args.genes_files], gene_name_list,
                                                        groups=[majiq_output['meta_exps'][0]['group'], None],
                                                        logger=logger)

    if type_summary == constants.ANALYSIS_DELTAPSI:
        threshold = args.threshold

        majiq_output = utils_voila.get_lsv_delta_exp_data(majiq_bins_file, args.confidence, args.threshold, args.show_all, logger=logger)
        majiq_output['event_list'] = []
        majiq_output['metadata'] = []
        for elem_list in majiq_output['genes_dict'].values():
            for elem in elem_list:
                majiq_output['event_list'].append(elem[0])
                majiq_output['metadata'].append(elem[1])
        # del majiq_output['genes_dict']

    if type_summary == constants.ANALYSIS_DELTAPSI_GENE:
        threshold = args.threshold

        import fileinput
        gene_name_list = []

        if args.gene_names:
            for gene_name in fileinput.input(args.gene_names):
                gene_name_list.append(gene_name.rstrip().upper())

        majiq_output = utils_voila.get_lsv_delta_exp_data(majiq_bins_file, args.confidence, args.threshold, args.show_all, gene_name_list=gene_name_list, logger=logger)

        if not gene_name_list:
            gene_name_list = majiq_output['genes_dict'].keys()

        if not gene_name_list:
            logger.warning("Number of LSVs detected in Voila with E(Delta(PSI)) > %.2f: None." % threshold)
            logger.warning("End of Voila execution.")
            return

        # Get gene info
        majiq_output['genes_exp'] = parse_gene_graphics([args.genesf_exp1, args.genesf_exp2], gene_name_list,
                                                        groups=[majiq_output['meta_exps'][0][0]['group'],
                                                                majiq_output['meta_exps'][1][0]['group']], logger=logger)

    if type_summary == constants.LSV_THUMBNAILS:
        try:
            majiq_output = []
            with open(majiq_bins_file, 'r') as types_file:
                for line in types_file:
                    majiq_output.append(line.rstrip())
        except IOError, e:
            logger.error(e.message, exc_info=1)

        meta_postprocess['collapsed'] = args.collapsed

    render_summary(output_dir, output_html, majiq_output, type_summary, threshold, meta_postprocess, logger=logger)
    create_gff3_txt_files(output_dir, majiq_output, logger=logger)
    render_tab_output(output_dir, output_html, majiq_output, type_summary, logger=logger)

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

    # Common script options
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('-b', '--bins', metavar='majiq_output.pickle', dest='majiq_bins', type=str, required=True, help='Pickle file with the bins produced by Majiq.')
    common_parser.add_argument('-o', '--output', metavar='output_dir', dest='output_dir', type=str, required=True, help='Output directory where the files will be placed.')
    common_parser.add_argument('--event-names', metavar='event_names.majiq', dest='event_names', type=str, help='Event names.')
    common_parser.add_argument('-c', '--confidence', metavar=0.95, dest='confidence', type=float, default=0.95, help='Percentage of confidence required (by default, 0.95).')
    common_parser.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    common_parser.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')

    # Subparser module to agglutinate all subparsers
    subparsers = parser.add_subparsers(dest='type_analysis')
    subparsers.required = True

    # Single LSV
    parser_single = argparse.ArgumentParser(add_help=False)
    parser_single.add_argument('--key-plots', metavar='keysplots.pickle', dest='keys_plots', type=str, help='Heatmap plots.')
    subparsers.add_parser(constants.ANALYSIS_PSI, help='Single LSV analysis.', parents=[common_parser, parser_single])

    # Delta LSV
    parser_delta = argparse.ArgumentParser(add_help=False)
    parser_delta.add_argument('--threshold', type=float, default=0.2, help='Filter out LSVs with no junction predicted to change over a certain value (in percentage).')  # Probability threshold used to sum the accumulative probability of inclusion/exclusion.
    parser_delta.add_argument('--show-all', dest='show_all', action='store_true', default=False, help='Show all LSVs including those with no junction with significant change predicted')
    subparsers.add_parser(constants.ANALYSIS_DELTAPSI, help='Delta LSV analysis.', parents=[common_parser, parser_delta])

    # Single LSV by Gene(s) of interest
    parser_single_gene = argparse.ArgumentParser(add_help=False)
    parser_single_gene.add_argument('--genes-files', nargs='+', required=True, dest='genes_files', metavar='Hippocampus1.splicegraph [Hippocampus2.splicegraph ...]', type=str, help='Splice graph information file(s) or directory with *.splicegraph file(s).')
    parser_single_gene.add_argument('--gene-names', type=str, dest='gene_names', help='Gene names to filter the results.')
    parser_single_gene.add_argument('--lsv-types', nargs='*', default=[], type=str, dest='lsv_types', help='LSV type to filter the results. (If no gene list is provided, this option will display only genes containing LSVs of the specified type).')
    subparsers.add_parser(constants.ANALYSIS_PSI_GENE, help='Single LSV analysis by gene(s) of interest.', parents=[common_parser, parser_single_gene])

    # Delta LSV by Gene(s) of interest
    parser_delta_gene = argparse.ArgumentParser(add_help=False)
    parser_delta_gene.add_argument('--genes-exp1', required=True, nargs='+', dest='genesf_exp1', metavar='Hippocampus1.splicegraph', type=str, help='Experiment 1 splice graph information file(s) or directory.')
    parser_delta_gene.add_argument('--genes-exp2', required=True, nargs='+', dest='genesf_exp2', metavar='Liver1.splicegraph', type=str, help='Experiment 2 splice graph information file(s) or directory.')
    parser_delta_gene.add_argument('--gene-names', type=str, dest='gene_names', help='Gene names to filter the results.')
    parser_delta_gene.add_argument('--lsv-types', nargs='*', default=[], type=str, dest='lsv_types', help='LSV type to filter the results. (If no gene list is provided, this option will display only genes containing LSVs of the specified type).')
    subparsers.add_parser(constants.ANALYSIS_DELTAPSI_GENE, help='Delta LSV analysis by gene(s) of interest.', parents=[common_parser, parser_delta, parser_delta_gene])

    # Thumbnails generation option (dev) TODO: Delete
    # parser_thumbs = argparse.ArgumentParser(add_help=False)
    # parser_thumbs.add_argument('--collapsed', type=bool, default=False, help='Collapsed LSVs thumbnails in the HTML summary.')
    # subparsers.add_parser(constanst.LSV_THUMBNAILS, help='Generate LSV thumbnails [DEBUGING!].', parents=[common_parser, parser_thumbs])  # TODO: GET RID OF THESE OR GIVE IT A BETTER SHAPE!!!

    args = parser.parse_args()
    try:
        create_summary(args)
    except ParseError, e:
        if e.logger:
            e.logger.error(repr(e), exc_info=constants.DEBUG)
        else:
            sys.stdout.write(repr(e))
    sys.exit(1)

if __name__ == '__main__':
    main()

