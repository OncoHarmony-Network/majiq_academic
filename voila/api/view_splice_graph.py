import math
from itertools import zip_longest

from voila import constants
from voila.api import SpliceGraph
from voila.api import splice_graph_model as model


class ViewSpliceGraph(SpliceGraph):
    def get_page_count(self, args):

        gene_count = self.session.query(model.Gene).count()
        if args.limit < gene_count:
            gene_count = args.limit

        return int(math.ceil(gene_count / float(constants.MAX_GENES)))

    def get_paginated_genes(self, args):
        def grouper(iterable, n, fillvalue=None):
            args = [iter(iterable)] * n
            return zip_longest(*args, fillvalue=fillvalue)

        for page in grouper(self.get_gene_ids(args), constants.MAX_GENES):
            yield page

    def get_gene_ids(self, args=None):
        if args and args.gene_ids:
            return args.gene_ids

        if args and hasattr(args, 'lsv_ids') and args.lsv_ids:
            return (lsv_id.split(':')[0] for lsv_id in args.lsv_ids)

        for gene_id, in self.session.query(model.Gene.id).all():
            yield gene_id
