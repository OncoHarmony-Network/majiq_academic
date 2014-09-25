#!/usr/bin/python
from collections import defaultdict
import json
import os
import subprocess
import sys
import utils_voila
import collections as cc
try:
    import cPickle as pkl
except:
    import pickle as pkl
from pdb import set_trace

__author__ = 'abarrera'

EXEC_DIR = os.path.dirname(os.path.realpath(__file__)) + "/"
VERSION = 'beta'

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

    env = Environment(extensions=["jinja2.ext.do"], loader=FileSystemLoader(EXEC_DIR + "templates/"))
    env.filters.update({'to_json': to_json, 'debug': debug})
    sum_template = env.get_template(type_summary + "_summary_template.html")

    if not output_dir.endswith('/'):
        output_dir += '/'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

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
            genes_dict = cc.OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)
            genes_json_dict = cc.OrderedDict((k, majiq_output['genes_json'][k]) for k in subset_keys)
            if (count_pages+1)*MAX_GENES < len(majiq_output['genes_dict']):
                print (count_pages+1)*MAX_GENES, len(majiq_output['genes_dict'])
                next_page = str(count_pages+1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            voila_output = open(output_dir+name_page, 'w')
            voila_output.write(sum_template.render( tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                                                    # metadata=majiq_output['metadata'],
                                                    genes_json=genes_json_dict,
                                                    genes_dict=genes_dict,
                                                    prevPage = prev_page,
                                                    nextPage= next_page,
                                                    namePage= name_page
            ))
            voila_output.close()
            count_pages += 1

    elif type_summary == 'lsv_delta_gene':
        # Max. 10 genes per page, create as many HTMLs as needed

        MAX_GENES = 10
        count_pages = 0

        # set_trace()
        gene_keys = sorted(majiq_output['genes_dict'].keys())
        while count_pages*MAX_GENES < len(gene_keys):
            prev_page = None
            next_page = None

            subset_keys = gene_keys[count_pages*MAX_GENES: MAX_GENES*(count_pages+1)]
            genes_dict = cc.OrderedDict((k, majiq_output['genes_dict'][k]) for k in subset_keys)

            if (count_pages+1)*MAX_GENES < len(majiq_output['genes_dict']):
                print (count_pages+1)*MAX_GENES, len(majiq_output['genes_dict'])
                next_page = str(count_pages+1) + "_" + output_html
            if not count_pages == 0:
                prev_page = str(count_pages - 1) + "_" + output_html

            name_page = str(count_pages) + "_" + output_html
            voila_output = open(output_dir+name_page, 'w')
            voila_output.write(sum_template.render( tableMarks=[table_marks_set(len(gene_set)) for gene_set in genes_dict],
                                                    # metadata=majiq_output['metadata'],
                                                    genes_dict=genes_dict,
                                                    genes_exps_list=majiq_output['genes_exp'],
                                                    expList=majiq_output['experiments_info'],
                                                    prevPage = prev_page,
                                                    nextPage= next_page,
                                                    namePage= name_page
            ))
            voila_output.close()
            count_pages += 1

    elif type_summary == 'lsv_delta':
        voila_output = open(output_dir+output_html, 'w')
        voila_output.write(sum_template.render( lsvList=majiq_output['event_list'],
                                                tableMarks=table_marks_set(len(majiq_output['event_list'])),
                                                metadata=majiq_output['metadata'],
                                                expList=majiq_output['experiments_info'],
                                                threshold=threshold
        ))
        voila_output.close()

    else:
        print "summary type not recognized %s." % type_summary
        import sys
        sys.exit(1)

    # Copy static files to the output directory
    utils_voila.copyanything(EXEC_DIR+"templates/static", output_dir+"static")


def parse_gene_graphics(gene_exps_flist, gene_name_list):
    genes_exp1_exp2 = []
    for gene_flist in gene_exps_flist:
        genes_exp = defaultdict()
        for splice_graph_f in gene_flist:
            genes_file = pkl.load(open(splice_graph_f, 'r'))
            genes_graphic = defaultdict(list)
            genes_file.sort()
            for gene_obj in genes_file:
                if gene_obj.get_name() in gene_name_list:
                    genes_graphic[gene_obj.get_name()].append(json.dumps(gene_obj, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'"))
                    genes_graphic[gene_obj.get_name()].append(gene_obj.get_strand())
                    genes_graphic[gene_obj.get_name()].append(gene_obj.get_coords())
                    genes_graphic[gene_obj.get_name()].append("1") #genes_graphic[gene_obj.get_name()].append(gene_obj.get_chrom())

            if not len(genes_graphic.keys()):
                raise Exception("[ERROR] :: No gene matching the splice graph file %s." % splice_graph_f)
            genes_exp[os.path.basename(splice_graph_f)] = genes_graphic
        genes_exp1_exp2.append(genes_exp)

    return genes_exp1_exp2


def create_tab_output(output_dir, output_html, majiq_output, type_summary, threshold, meta_postprocess):
    delimiter = '\t'
    ofile_str = output_dir + output_html.rsplit('.html', 1)[0] + '.csv'
    tlb_categx = {'A5SS': 'prime5', 'A3SS': 'prime3', 'Num. Junctions': 'njuncs', 'Num. Exons': 'nexons'}
    with open(ofile_str, 'w') as ofile:
        headers = ['Gene name', 'LSV ID', 'E(Delta(PSI)) per LSV way', 'P(Delta(PSI)>%s) per LSV way' % .1, 'LSV Type', 'A5SS', 'A3SS', 'ES', 'Num. Junctions', 'Num. Exons', 'chr', 'strand', 'Junctions coords', 'Exons coords', 'Exons Alternative Start', 'Exons Aletrnative End']
        ofile.write(delimiter.join(headers))
        for gene in majiq_output['genes_dict']:
            lline = []
            for llsv in majiq_output['genes_dict'][gene]:
                lline.extend([gene, llsv[0].get_id()])
                lexpected = []
                lconfidence = []
                for i, bins in enumerate(llsv[0].get_bins()):
                    if 'delta' in type_summary:
                        lexpected.append(sum(llsv[0].get_excl_incl()[i]))
                        lconfidence.append(utils_voila.get_prob_delta_psi_greater_v(bins, lexpected, 0.1))
                    else:
                        lexpected.append(llsv[0].get_means_psi()[i])
                        lconfidence.append(llsv[0].get_variances()[i])

                lline.append('; '.join(lexpected))
                lline.append('; '.join(lconfidence))

                lline.append(llsv[0].get_type())
                lline.append(llsv[0].get_categories()[tlb_categx['A5SS']])
                lline.append(llsv[0].get_categories()[tlb_categx['A3SS']])
                lline.append(llsv[0].get_categories()[tlb_categx['ES']])
                lline.append(llsv[0].get_categories()[tlb_categx['Num. Junctions']])
                lline.append(llsv[0].get_categories()[tlb_categx['Num. Exons']])

                lline.append(llsv[1].get_chrom())
                lline.append(llsv[1].get_strand())

                lline.append('; '.join(['-'.join(str(c) for c in junc.get_coords()) for junc in llsv[1].get_junctions()]))
                lline.append('; '.join(['-'.join(str(c) for c in exon.get_coords()) for exon in llsv[1].get_exons()]))

                try:
                    lline.append('; '.join([' '.join([str(c) for c in exon.get_alt_starts()]) for exon in llsv[1].get_exons()]))
                    lline.append('; '.join([' '.join([str(c) for c in exon.get_alt_ends()]) for exon in llsv[1].get_exons()]))

                except Exception:
                    pass

            ofile.write(delimiter.join(lline))
    print "Delimited output file successfully created in:\n%s" % ofile_str


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

    if type_summary == 'lsv_single':
        majiq_output = utils_voila.get_lsv_single_exp_data(majiq_bins_file, args.confidence)

    if type_summary == 'lsv_single_gene':

        lsv_types = args.lsv_types
        bed_file  = args.bed_file

        genes_file = pkl.load(open(args.genes_file, 'r'))

        import fileinput
        gene_name_list = []
        if args.gene_names:
            for gene_name in fileinput.input(args.gene_names):
                gene_name_list.append(gene_name.rstrip().split(":")[0])
        else:
            gene_name_list = []
        majiq_output = utils_voila.get_lsv_single_exp_data(majiq_bins_file, args.confidence, gene_name_list=gene_name_list, lsv_types=lsv_types)  #, bed_file=bed_file

        # Get gene info
        # try:
        genes_graphic = defaultdict(list)
        for gene_obj in genes_file:
            if gene_obj.get_name() in majiq_output['genes_dict']:
                genes_graphic[gene_obj.get_name()].append(json.dumps(gene_obj, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'"))
                genes_graphic[gene_obj.get_name()].append(gene_obj.get_strand())
                genes_graphic[gene_obj.get_name()].append(gene_obj.get_coords())
                # genes_graphic[gene_obj.get_name()].append(gene_obj.get_chrom())
                # majiq_output['gene_json'] = json.dumps(gene_obj, cls=utils_voila.LsvGraphicEncoder).replace("\"", "'")
                # majiq_output['gene'] = gene_obj


                for lsv_data in majiq_output['genes_dict'][gene_obj.get_name()]:
                    # Find which is the ending coordinate of the LSV
                    lsv_data[1].append(gene_obj.get_exons()[-1].get_coords()[1])

                    # Calculate extension of LSV
                    lsv_data[0].set_coords(lsv_data[1][0])
                    lsv_data[0].set_extension(lsv_data[1][4], lsv_data[1][2])

        if not len(genes_graphic.keys()): raise Exception("[ERROR] :: No gene matching the visual information file.")
        majiq_output['genes_json'] = genes_graphic
        # majiq_output['genes'] = genes_graphic[1]

        # except Exception, e:
        #     print e.message
        #     raise e

    if type_summary == 'lsv_delta':
        threshold = args.threshold

        majiq_output = utils_voila.get_lsv_delta_exp_data(majiq_bins_file, args.confidence, args.threshold, args.show_all)
        majiq_output['event_list'] = []
        majiq_output['metadata'] = []
        for elem_list in majiq_output['genes_dict'].values():
            for elem in elem_list:
                majiq_output['event_list'].append(elem[0])
                majiq_output['metadata'].append(elem[1])
        del majiq_output['genes_dict']

    if type_summary == 'lsv_delta_gene':  # For a Delta PSI gene summary, there might be more than one splicegraph per experiment
        threshold = args.threshold

        import fileinput
        gene_name_list = []

        if args.gene_names:
            for gene_name in fileinput.input(args.gene_names):
                gene_name_list.append(gene_name.rstrip())

        majiq_output = utils_voila.get_lsv_delta_exp_data(majiq_bins_file, args.confidence, args.threshold, args.show_all, gene_name_list=gene_name_list)

        if not gene_name_list:
            gene_name_list = majiq_output['genes_dict'].keys()

        # Get gene info
        majiq_output['genes_exp'] = parse_gene_graphics([args.genesf_exp1, args.genesf_exp2], gene_name_list)

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

    _render_template(output_dir, output_html, majiq_output, type_summary, threshold, meta_postprocess)
    create_tab_output(output_dir, output_html, majiq_output, type_summary, threshold, meta_postprocess)

    print "Summaries created in:\n%s" % output_dir
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
    parser_single_gene.add_argument('--genes-file-info', required=True, dest='genes_file', metavar='Hippocampus1.splicegraph', type=str, help='Splice graph information file.')
    parser_single_gene.add_argument('--gene-names', type=str, dest='gene_names', help='Gene names to filter the results.')
    parser_single_gene.add_argument('--bed-file', type=str, dest='bed_file', help='Bed file with coordinates and genes mapping.')
    parser_single_gene.add_argument('--lsv-types', nargs='*', default=[], type=str, dest='lsv_types', help='LSV type to filter the results. (If no gene list is provided, this option will display only genes containing LSVs of the specified type).')
    subparsers.add_parser('lsv-single-gene', help='Single LSV analysis by gene(s) of interest.', parents=[common_parser, parser_single_gene])

    # Delta LSV by Gene(s) of interest
    parser_delta_gene = argparse.ArgumentParser(add_help=False)
    parser_delta_gene.add_argument('--genes-exp1', required=True, nargs='+', dest='genesf_exp1', metavar='Hippocampus1.splicegraph', type=str, help='Splice graph information file(s) of experiment 1.')
    parser_delta_gene.add_argument('--genes-exp2', required=True, nargs='+', dest='genesf_exp2', metavar='Liver1.splicegraph', type=str, help='Splice graph information file(s) of experiment 2.')
    parser_delta_gene.add_argument('--gene-names', type=str, dest='gene_names', help='Gene names to filter the results.')
    parser_delta_gene.add_argument('--bed-file', type=str, dest='bed_file', help='Bed file with coordinates and genes mapping.')
    parser_delta_gene.add_argument('--lsv-types', nargs='*', default=[], type=str, dest='lsv_types', help='LSV type to filter the results. (If no gene list is provided, this option will display only genes containing LSVs of the specified type).')
    subparsers.add_parser('lsv-delta-gene', help='Delta LSV analysis by gene(s) of interest.', parents=[common_parser, parser_delta, parser_delta_gene])

    # Thumbnails generation option (dev) TODO: Delete
    # parser_thumbs = argparse.ArgumentParser(add_help=False)
    # parser_thumbs.add_argument('--collapsed', type=bool, default=False, help='Collapsed LSVs thumbnails in the HTML summary.')
    # subparsers.add_parser('lsv-thumbnails', help='Generate LSV thumbnails [DEBUGING!].', parents=[common_parser, parser_thumbs])  # TODO: GET RID OF THESE OR GIVE IT A BETTER SHAPE!!!


    args = parser.parse_args()
    print repr(args)

    create_summary(args)



if __name__ == '__main__':
    main()

