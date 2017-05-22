import json
from collections import OrderedDict
from distutils.dir_util import copy_tree
from itertools import zip_longest

from math import ceil
from os import path

from jinja2 import Environment, FileSystemLoader
from markupsafe import escape

from voila import constants
from voila import module_locator
from voila.splice_graphics import SpliceGraph
from voila.utils.voila_log import voila_log

EXEC_DIR = module_locator.module_path()


def parse_gene_graphics(splice_graph_file, metainfo, gene_ids_list):
    """
    Load and combine splice graph files.

    :param splicegraph_flist: list of splice graph files or directory containing splice graphs.
    :param gene_names_list: list of genes of interest.
    :param condition_names: ids for condition 1 [and condition 2, in deltapsi].
    :return: list of genes graphic per condition.
    """
    log = voila_log()
    log.info("Parsing splice graph information files ...")

    genes_exp1_exp2 = []

    with SpliceGraph(splice_graph_file, 'r') as sg:
        genes = sg.get_genes_list(gene_ids_list)
        gene_experiments_list = sg.get_experiments_list()

    genes.sort()

    for experiments in [metainfo['experiments1'], metainfo.get('experiments2', [])]:
        genes_exp = {}
        combined_genes_exp = {}

        for experiment in experiments:
            genes_exp[experiment] = {}

            for gene in genes:
                # map the metainfo experiment name to the experiment index in the splice graph file.
                experiment_index = gene_experiments_list.index(experiment)

                # get the data needed to render the html
                genes_exp[experiment][gene.gene_id] = gene.get_experiment(experiment_index)

                # record all genes and combine their experiment data
                combined_genes_exp[gene.gene_id] = gene.combine(experiment_index,
                                                                combined_genes_exp.get(gene.gene_id, None))

        # if there are more then 1 experiments, then record the combined data
        if len(experiments) > 1:
            genes_exp['Combined'] = {gene_id: combined_genes_exp[gene_id] for gene_id in combined_genes_exp}

        genes_exp1_exp2.append(OrderedDict(sorted(genes_exp.items(), key=lambda t: t[0])))

    log.info("Splice graph information files correctly loaded.")

    return genes_exp1_exp2


class VoilaNoLSVsException(Exception):
    def __init__(self):
        """
        Throw when there are no LSVs in voila file.
        """
        message = "There are no LSVs detected.  It could be the threshold is too high or the filters are " \
                  "incorrectly set."
        voila_log().error(message)
        super(VoilaNoLSVsException, self).__init__(message)


def table_marks_set(size):
    """
    Calculate the number of elements to show in LSV tables.

    :param size: total number of LSVs.
    :return: set of total number of elements to show.
    """
    ideal_set = (10, 20, 50, 100)
    for index, value in enumerate(ideal_set):
        if size < value:
            return ideal_set[0:index]
    return ideal_set


def get_env():
    """
    Get environment variable for creating html files from jinja2 templates.
    :return: env variable
    """

    def to_json(value):
        return json.dumps(value.to_dict()).replace('"', '\'')

    def to_json_especial(value):
        return escape(json.dumps(value).replace('\"', '\''))

    env = Environment(extensions=["jinja2.ext.do"], loader=FileSystemLoader(path.join(EXEC_DIR, "templates/")))
    env.filters.update({'to_json': to_json, 'to_json_especial': to_json_especial})
    return env


def get_summary_template(args, env):
    """
    Get summary template for specific type analysis.
    :param args: command line arguments
    :param env: environment variable
    :return: summary template
    """
    template_file_name = args.type_analysis.replace("-", "_") + "_summary_template.html"
    return env.get_template(template_file_name)


def get_output_html(args, file_name=None):
    """
    Get output html file name.
    :param args: command line arguments
    :param file_name: input file name
    :return:
    """

    if file_name:
        return '{0}_{1}.html'.format(path.splitext(path.split(file_name)[1])[0], args.type_analysis.replace("-", "_"))
    else:
        return '{0}.html'.format(args.type_analysis.replace("-", "_"))


def grouper(iterable, n, fillvalue=None):
    """
    Pulled from python api... https://docs.python.org/2/library/itertools.html#recipes

    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    :param iterable:
    :param n:
    :param fillvalue:
    :return:
    """

    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def copy_static(args):
    """
    Copy static files to output directory.
    :param args: command line arguments
    :return: None
    """
    log = voila_log()
    log.info("Copying static files from Voila sources ...")
    static_dir = path.join(EXEC_DIR, 'templates/static')
    copy_tree(static_dir, path.join(args.output, 'static'))
    copy_tree(static_dir, path.join(args.output, constants.SUMMARIES_SUBFOLDER, 'static'))


def get_prev_next_pages(page_number, genes_count, output_html, limit=None):
    """
    Get location of prev and next pages for rendering html files.
    :param page_number: current page number
    :param genes_count: number of genes
    :param output_html: html file name
    :param limit: limit for genes
    :return: prev page, next page
    """
    if limit:
        genes_count = min(genes_count, limit)

    last_page = ceil(genes_count / float(constants.MAX_GENES)) - 1

    next_page = None
    prev_page = None

    if page_number != last_page:
        next_page = '{0}_{1}'.format(page_number + 1, output_html)

    if page_number != 0:
        prev_page = '{0}_{1}'.format(page_number - 1, output_html)

    return prev_page, next_page
