import os
from collections import OrderedDict

from voila import constants
from voila.utils import utils_voila
from voila.utils.run_voila_utils import get_env, get_output_html


class Html(object):
    def __init__(self, args):
        self.args = args
        self.voila_links = {}
        self.env = get_env()

    def get_summaries_subfolder(self):
        summaries_subfolder = os.path.join(self.args.output, constants.SUMMARIES_SUBFOLDER)
        utils_voila.create_if_not_exists(summaries_subfolder)
        return summaries_subfolder

    @staticmethod
    def gene_experiments(experiments, genes, gene_experiments):

        genes_exp = {}
        combined_genes_exp = {}
        gene_experiments_tuple = tuple(gene_experiments)

        for experiment in experiments:
            genes_exp[experiment] = {}

            for gene in genes:
                # map the metainfo experiment name to the experiment index in the splice graph file.
                experiment_index = gene_experiments_tuple.index(experiment)

                # get the data needed to render the html
                genes_exp[experiment][gene.gene_id] = gene.get_experiment(experiment_index)

                # record all genes and combine their experiment data
                combined_genes_exp[gene.gene_id] = gene.combine(experiment_index,
                                                                combined_genes_exp.get(gene.gene_id, None))

        # if there are more then 1 experiments, then record the combined data
        if len(experiments) > 1:
            genes_exp['Combined'] = {gene_id: combined_genes_exp[gene_id] for gene_id in combined_genes_exp}

        return OrderedDict(sorted(genes_exp.items(), key=lambda t: t[0]))

    def add_to_voila_links(self, lsv_dict, page_name):
        for lsvs in lsv_dict.values():
            for lsv in lsvs:
                self.voila_links[lsv.name] = '{0}#{1}'.format(os.path.join(constants.SUMMARIES_SUBFOLDER, page_name),
                                                              lsv.name)

    def get_page_name(self, index):
        args = self.args
        try:
            output_html = get_output_html(args, args.voila_file)
        except AttributeError:
            output_html = get_output_html(args, args.splice_graph)
        return '{0}_{1}'.format(index, output_html)

    def get_next_page(self, index, page_count):
        if index + 1 == page_count:
            return None
        else:
            return self.get_page_name(index + 1)
