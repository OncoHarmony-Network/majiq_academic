#!/usr/bin/python
from collections import defaultdict
import json
import os
import sys
import utils_voila
import collections as cc


__author__ = 'abarrera'

EXEC_DIR = os.path.dirname(os.path.realpath(__file__)) + "/"
VOILA_ANALYSIS_TYPES = ['single', 'delta', 'lsv_single', 'lsv_delta', 'lsv_thumbnails', 'lsv_gene']
VERSION = 'alpha'

def table_marks_set(size):
    # TODO: Customizable by script options
    ideal_set = (10, 20, 50, 100)
    index = 0
    for mark in ideal_set:
        if size < ideal_set[index]:
            break
        index += 1

    return ideal_set[0:index]


def mocked_get_array_bins():
    pass

def debug(text):
    print text
    return ''

def _render_template(output_dir, output_html, majiq_output, type_summary, threshold=None, post_process_info=None):
    """
    Rendering the summary template to create the HTML file.

    @param output_dir: directory where the HTML will be created
    @param majiq_output: event list from MAJIQ output used as input in Voila
    @param type_summary: defines summary template used
    """
    from jinja2 import Environment, FileSystemLoader, escape
    def to_json(value):
        return escape(json.dumps(value, cls=utils_voila.PickleEncoder))

    env = Environment(loader=FileSystemLoader(EXEC_DIR + "templates/"))
    env.filters.update({'to_json': to_json, 'debug': debug})
    sum_template = env.get_template(type_summary + "_summary_template.html")

    if not output_dir.endswith('/'):
        output_dir += '/'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if type_summary == 'single':
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render(eventList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               metadata=majiq_output['metadata_pre']))
        voila_output.close()
    elif type_summary == 'delta':
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render(eventList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               expList=majiq_output['experiments_info'],
                                               threshold=threshold
                                               ))
        voila_output.close()

    elif type_summary == 'lsv_single':
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render(lsvList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               metadata=majiq_output['metadata']
                                               ))
        voila_output.close()
    elif type_summary == 'lsv_thumbnails':
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render(lsvList=majiq_output,
                                               collapsed=post_process_info['collapsed']))
        voila_output.close()

    elif type_summary == 'lsv_single_gene':
        # Max. 10 genes per page, create as many HTMLs as needed

        MAX_GENES = 10
        count_pages = 0
        gene_keys = sorted(majiq_output['genes_dict'].keys())
        while count_pages*MAX_GENES < len(majiq_output['genes_json']):
            prev_page = None
            next_page = None

            subset_keys = gene_keys[count_pages*MAX_GENES: MAX_GENES*(count_pages+1)]

            genes_dict = cc.OrderedDict(dict((k, majiq_output['genes_dict'][k]) for k in subset_keys).items())
            genes_json_dict = cc.OrderedDict(dict((k, majiq_output['genes_json'][k]) for k in subset_keys).items())
            if (count_pages+1)*MAX_GENES < len(majiq_output['genes_dict']):
                print (count_pages+1)*MAX_GENES, len(majiq_output['genes_dict'])
                next_page = str(count_pages+1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            voila_output = open(output_dir+str(count_pages) + "_" + output_html, 'w')
            voila_output.write(sum_template.render( tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                                                    # metadata=majiq_output['metadata'],
                                                    genes_json=genes_json_dict,
                                                    genes_dict=genes_dict,
                                                    prevPage = prev_page,
                                                    nextPage= next_page
            ))
            voila_output.close()
            count_pages += 1

    elif type_summary == 'lsv_delta':
        if 'genes_dict' not in majiq_output:
            voila_output = open(output_dir+output_html, 'w')
            voila_output.write(sum_template.render( lsvList=majiq_output['event_list'],
                                                    tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                                    metadata=majiq_output['metadata'],
                                                    threshold=threshold
            ))
            voila_output.close()

    else:
        print "summary type not recognized %s." % type_summary
        import sys
        sys.exit(1)

    # Copy static files to the output directory
    utils_voila.copyanything(EXEC_DIR+"templates/static", output_dir+"static")


def create_summary(args):
    """This method generates an html summary from a majiq output file"""

    type_summary    = args.type_analysis.replace('-', '_')  # Notation preference
    majiq_bins_file = args.majiq_bins
    output_dir      = args.output_dir
    if not str(output_dir).endswith('/'):
        output_dir += '/'
    meta_preprocess = args.meta_preprocess
    threshold       = None

    output_html = os.path.splitext(os.path.split(majiq_bins_file)[1])[0] + "_" + type_summary + '.html'
    majiq_output = None
    meta_postprocess = {}

    # if type_summary == 'single':
    #     majiq_output = utils_voila.get_single_exp_data(majiq_bins_file, meta_preprocess, meta_postprocess, confidence)
    # elif type_summary == 'delta':
    #     majiq_output = utils_voila.get_delta_exp_data(majiq_bins_file, meta_postprocess, confidence, threshold)
    # el
    if type_summary == 'lsv_single':
        majiq_output = utils_voila.get_lsv_single_exp_data(majiq_bins_file, args.confidence)
    if type_summary == 'lsv_thumbnails':
        try:
            majiq_output = []
            with open(majiq_bins_file, 'r') as types_file:
                for line in types_file:
                    majiq_output.append(line.rstrip())
        except IOError, e:
            print e.message
            sys.exit(1)
        meta_postprocess['collapsed'] = args.collapsed

    if type_summary == 'lsv_single_gene':
        if not args.genes_file:
            print "[ERROR] :: parameter genes-file-info needed for filtering results by gene name."
            sys.exit(1)
        import pickle as pkl
        if args.genes_file:
            genes_file = pkl.load(open(args.genes_file, 'r'))

            import fileinput
            gene_name_list = []
            if args.gene_names:
                for gene_name in fileinput.input(args.gene_names):
                    gene_name_list.append(gene_name.rstrip().split(":")[0])
        else:
            gene_name_list = []
        majiq_output = utils_voila.get_lsv_single_exp_data(majiq_bins_file, args.confidence, gene_name_list=gene_name_list, lsv_type=args.lsv_type )

        # Get gene info
        try:
            genes_graphic = defaultdict(list)
            for gene_obj in genes_file:
                if gene_obj.get_name() in gene_name_list:
                    genes_graphic[gene_obj.get_name()].append(json.dumps(gene_obj, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'"))
                    genes_graphic[gene_obj.get_name()].append(gene_obj.get_strand())
                    # majiq_output['gene_json'] = json.dumps(gene_obj, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'")
                    # majiq_output['gene'] = gene_obj

            if not len(genes_graphic.keys()): raise Exception("[ERROR] :: No gene matching the visual information file.")
            majiq_output['genes_json'] = genes_graphic
            # majiq_output['genes'] = genes_graphic[1]

        except Exception, e:
            print e.message
            raise e


    if type_summary == 'lsv_delta':
        threshold = args.threshold
        if 'gene_names' in args:
            if 'genes_file' not in args :
                print "[ERROR] :: parameter genes-file-info needed for filtering results by gene name."
                sys.exit(1)
            import pickle as pkl
            genes_file = pkl.load(open(args.genes_file, 'r'))

            import fileinput
            gene_name_list = []
            for gene_name in fileinput.input(args.gene_names):
                gene_name_list.append(gene_name.rstrip())

            majiq_output = utils_voila.get_lsv_delta_exp_data(majiq_bins_file, args.confidence, args.threshold, args.show_all)

            # Get gene info
            try:
                genes_graphic = defaultdict(list)
                for gene_obj in genes_file:
                    if gene_obj.get_name() in gene_name_list:
                        genes_graphic[gene_obj.get_name()].append(json.dumps(gene_obj, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'"))
                        genes_graphic[gene_obj.get_name()].append(gene_obj.get_strand())

                if not len(genes_graphic.keys()): raise Exception("[ERROR] :: No gene matching the splice graph information file.")
                majiq_output['genes_json'] = genes_graphic

            except Exception, e:
                print e.message
                sys.exit(1)
        else:
            majiq_output = utils_voila.get_lsv_delta_exp_data(majiq_bins_file, args.confidence, args.threshold)
            majiq_output['event_list'] = []
            majiq_output['metadata'] = []
            for elem_list in majiq_output['genes_dict'].values():
                for elem in elem_list:
                    majiq_output['event_list'].append(elem[0])
                    majiq_output['metadata'].append(elem[1])
            del majiq_output['genes_dict']

    _render_template(output_dir, output_html, majiq_output, type_summary, threshold, meta_postprocess)
    print "Summary created in:\n%s" % output_dir
    return


def main():

    import argparse
    parser = argparse.ArgumentParser(description="VOILA is a visualization package for Alternative Splicing Events and Quantifications.")
    parser.add_argument('-v', action='version', version=VERSION)

    # Common script options
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('-b', '--bins', metavar='majiq_output.pickle', dest='majiq_bins', type=str, required=True, help='Pickle file with the bins produced by Majiq.')
    common_parser.add_argument('-o', '--output', metavar='output_dir', dest='output_dir', type=str, required=True, help='Output directory where the files will be placed.')
    common_parser.add_argument('--metadata_builder', metavar='metadata_pre.majiq', dest='meta_preprocess', type=str, help='Metadata from MAJIQ builder.')
    common_parser.add_argument('--event-names', metavar='event_names.majiq', dest='event_names', type=str, help='Event names.')
    common_parser.add_argument('-c', '--confidence', metavar=0.95, dest='confidence', type=float, default=0.95, help='Percentage of confidence required (by default, 0.95).')

    # Subparser module to agglutinate all subparsers
    subparsers = parser.add_subparsers(dest='type_analysis')
    subparsers.required = True

    # Single LSV
    parser_single = argparse.ArgumentParser(add_help=False)
    parser_single.add_argument('--key-plots', metavar='keysplots.pickle', dest='keys_plots', type=str, help='Heatmap plots.')
    subparsers.add_parser('lsv-single', help='Single LSV analysis.', parents=[common_parser, parser_single])

    # Delta LSV
    parser_delta = argparse.ArgumentParser(add_help=False)
    parser_delta.add_argument('--threshold-significant', type=float, dest='threshold', default=0.2, help='Filter out LSVs with no junction predicted to change over a certain value (in percentage).')  # Probability threshold used to sum the accumulative probability of inclusion/exclusion.
    parser_delta.add_argument('--show-all', dest='show_all', action='store_true', default=False, help='Show all LSVs including those with no junction with significant change predicted')
    subparsers.add_parser('lsv-delta', help='Delta LSV analysis.', parents=[common_parser, parser_delta])

    # Single LSV by Gene(s) of interest
    parser_single_gene = argparse.ArgumentParser(add_help=False)
    parser_single_gene.add_argument('--genes-file-info', dest='genes_file', metavar='visual_LSE.majiq', type=str, help='Pickle file with gene coords info.')
    parser_single_gene.add_argument('--gene-names', type=str, dest='gene_names', help='Gene names to filter the results.')
    parser_single_gene.add_argument('--lsv-type', type=str, dest='lsv_type', help='LSV type to filter the results. (If no gene list is provided, this option will display only genes containing LSVs of the specified type).')
    subparsers.add_parser('lsv-single-gene', help='Single LSV analysis by gene(s) of interest.', parents=[common_parser, parser_single_gene])

    # Thumbnails generation option (dev) TODO: Delete
    parser_thumbs = argparse.ArgumentParser(add_help=False)
    parser_thumbs.add_argument('--collapsed', type=bool, default=False, help='Collapsed LSVs thumbnails in the HTML summary.')
    subparsers.add_parser('lsv-thumbnails', help='Generate LSV thumbnails [DEBUGING!].', parents=[common_parser, parser_thumbs])  # TODO: GET RID OF THESE OR GIVE IT A BETTER SHAPE!!!


    args = parser.parse_args()
    print args

    create_summary(args)



if __name__ == '__main__':
    main()

