#!/usr/bin/python
from collections import defaultdict
import json
import os
import sys
import utils_voila


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

def _render_template(output_dir, output_html, majiq_output, type_summary, threshold, post_process_info=None):
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

    voila_output = open(output_dir+output_html, 'w+')

    if type_summary == 'single':
        voila_output.write(sum_template.render(eventList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               metadata=majiq_output['metadata_pre']))
    elif type_summary == 'delta':
        voila_output.write(sum_template.render(eventList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               expList=majiq_output['experiments_info'],
                                               threshold=threshold
                                               ))

    elif type_summary == 'lsv_single':
        voila_output.write(sum_template.render(lsvList=majiq_output['event_list'],
                                               tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                               metadata=majiq_output['metadata']
                                               ))
    elif type_summary == 'lsv_thumbnails':
        voila_output.write(sum_template.render(lsvList=majiq_output,
                                               collapsed=post_process_info['collapsed']))

    elif type_summary == 'lsv_gene':
        voila_output.write(sum_template.render( lsvList=majiq_output['event_list'],
                                                tableMarks=[table_marks_set(len(gene_set)) for gene_set in majiq_output['genes_dict']],
                                                metadata=majiq_output['metadata'],
                                                genes_json=majiq_output['genes_json'],
                                                genes_dict=majiq_output['genes_dict']
        ))

    elif type_summary == 'lsv_delta':
        if 'genes_dict' not in majiq_output:
            voila_output.write(sum_template.render( lsvList=majiq_output['event_list'],
                                                    tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                                    metadata=majiq_output['metadata'],
                                                    threshold=threshold
            ))

    else:
        print "summary type not recognized %s." % type_summary
        import sys
        sys.exit(1)

    # Copy static files to the output directory
    utils_voila.copyanything(EXEC_DIR+"templates/static", output_dir+"static")
    voila_output.close()


def create_summary(args):
    """This method generates an html summary from a majiq output file"""

    type_summary    = args.type_analysis.replace('-', '_')  # Notation preference
    majiq_bins_file = args.majiq_bins
    output_dir      = args.output_dir
    meta_preprocess = args.meta_preprocess

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
        genes_file = pkl.load(open(args.genes_file, 'r'))

        import fileinput
        gene_name_list = []
        for gene_name in fileinput.input(args.gene_names):
            gene_name_list.append(gene_name.rstrip())

        majiq_output = utils_voila.get_lsv_single_exp_data(majiq_bins_file, args.confidence, gene_name_list=gene_name_list)

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
            sys.exit(1)

    if type_summary == 'lsv_delta':
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

    _render_template(output_dir, output_html, majiq_output, type_summary, args.threshold, meta_postprocess)
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
    parser_single_gene.add_argument('--gene-names', type=str, dest='gene_names', help='Gene names to filter the results. ')
    subparsers.add_parser('lsv-single-gene', description='Single LSV analysis by gene(s) of interest.', parents=[common_parser, parser_single_gene])

    # Thumbnails generation option (dev) TODO: Delete
    parser_thumbs = argparse.ArgumentParser(add_help=False)
    parser_thumbs.add_argument('--collapsed', type=bool, default=False, help='Collapsed LSVs thumbnails in the HTML summary.')
    subparsers.add_parser('lsv-thumbnails', help='Generate LSV thumbnails [DEBUGING!].', parents=[common_parser, parser_thumbs])  # TODO: GET RID OF THESE OR GIVE IT A BETTER SHAPE!!!


    args = parser.parse_args()
    print args

    create_summary(args)



if __name__ == '__main__':
    main()

