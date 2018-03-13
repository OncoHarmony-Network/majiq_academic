import math
from itertools import zip_longest

from voila import constants
from voila.api import splice_graph_model
from voila.api import view_splice_graph_model as model
from voila.api.sql import SQL


class ViewSpliceGraph(SQL):
    def __init__(self, args):
        super().__init__(args.splice_graph, model)
        self.args = args

    @property
    def genome(self):
        return self.session.query(splice_graph_model.Genome.name).one()[0]

    @property
    def experiment_names(self):
        return (e for e, in self.session.query(splice_graph_model.Experiment.name).all())

    def genes(self, gene_ids):
        for gene_id in gene_ids:
            yield self.gene(gene_id)

    def gene(self, gene_id):
        return self.session.query(model.ViewGene).get(gene_id)

    def get_paginated_genes(self):
        def grouper(iterable, n, fillvalue=None):
            a = [iter(iterable)] * n
            return zip_longest(*a, fillvalue=fillvalue)

        for page in grouper(self.get_genes(), constants.MAX_GENES):
            yield tuple(p for p in page if p is not None)

    def get_genes(self):
        args = self.args
        if args.gene_ids:
            for gene_id in args.gene_ids:
                yield self.gene(gene_id)
        else:
            return self.session.query(model.ViewGene).all()

    def get_page_count(self):
        """
        This page count is only for the splice graph view. Not to be used with anything containing quantified data.
        :param args: arg parse
        :return: number of pages
        """
        gene_count = len(tuple(self.get_genes()))
        return math.ceil(gene_count / constants.MAX_GENES)
