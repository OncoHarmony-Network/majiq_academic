import os
import uuid
from collections import OrderedDict

from voila import constants
from voila.api.view_splice_graph import HtmlSpliceGraph
from voila.utils import utils_voila
from voila.utils.run_voila_utils import get_env, get_output_html


class Html(object):
    def __init__(self, args):
        self.args = args
        self.voila_links = {}
        self.env = get_env()
        self._database_name = None

    def gene_experiments(self, experiments, genes, gene_experiments):

        genes_exp = {}
        gene_experiments_tuple = tuple(gene_experiments)

        with HtmlSpliceGraph(self.args.splice_graph) as sg:
            for experiment in experiments:

                genes_exp[experiment] = {}

                for gene_id in genes:
                    gene = sg.gene(gene_id)

                    # map the metainfo experiment name to the experiment index in the splice graph file.
                    experiment_index = gene_experiments_tuple.index(experiment)

                    # get the data needed to render the html
                    genes_exp[experiment][gene_id] = gene.get.get_experiment(experiment_index)

            # if there are more then 1 experiments, then record the combined data
            if len(experiments) > 1:
                genes_exp['Combined'] = sg.combine_genes(experiments, genes)

            return OrderedDict(sorted(genes_exp.items(), key=lambda t: t[0]))

    def add_to_voila_links(self, lsv_dict, page_name):
        for gene_id in lsv_dict.keys():
            self.voila_links[gene_id] = '{0}#{1}'.format(os.path.join(constants.SUMMARIES_SUBFOLDER, page_name),
                                                         gene_id)

    def database_name(self):
        if self._database_name is None:
            self._database_name = 'voila_{}'.format(uuid.uuid4().hex)

        return self._database_name

    @staticmethod
    def get_summaries_subfolder(args):
        summaries_subfolder = os.path.join(args.output, constants.SUMMARIES_SUBFOLDER)
        utils_voila.create_if_not_exists(summaries_subfolder)
        return summaries_subfolder

    @staticmethod
    def voila_links(lsv_dict, page_name):
        for gene_id in lsv_dict.keys():
            yield gene_id, '{0}#{1}'.format(os.path.join(constants.SUMMARIES_SUBFOLDER, page_name), gene_id)

    @staticmethod
    def get_page_name(args, index):
        try:
            output_html = get_output_html(args, args.voila_file)
        except AttributeError:
            output_html = get_output_html(args, args.splice_graph)
        return '{0}_{1}'.format(index, output_html)

    @classmethod
    def get_next_page(cls, args, index, page_count):
        if index + 1 == page_count:
            return None
        else:
            return cls.get_page_name(args, index + 1)

    @classmethod
    def get_prev_page(cls, args, index):
        if index - 1 < 0:
            return None
        else:
            return cls.get_page_name(args, index - 1)

    @staticmethod
    def chunkify(lst, n):
        for i in range(n):
            yield lst[i::n]