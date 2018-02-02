import math
from itertools import zip_longest

from voila import constants
from voila.api import SpliceGraph
from voila.api import splice_graph_model as model
from voila.api.matrix_hdf5 import lsv_id_to_gene_id


class ViewSpliceGraph(SpliceGraph):
    def get_paginated_genes(self, args):
        def grouper(iterable, n, fillvalue=None):
            args = [iter(iterable)] * n
            return zip_longest(*args, fillvalue=fillvalue)

        for page in grouper(self.get_gene_ids(args), constants.MAX_GENES):
            yield tuple(p for p in page if p is not None)

    def get_gene_ids(self, args):
        if args.gene_ids:
            yield from args.gene_ids
        elif hasattr(args, 'lsv_ids') and args.lsv_ids:
            for lsv_id in args.lsv_ids:
                yield lsv_id_to_gene_id(lsv_id)
        else:
            for gene_id, in self.session.query(model.Gene.id).all():
                yield gene_id


class SpliceGraphViewSpliceGraph(ViewSpliceGraph):
    def get_page_count(self, args):
        """
        This page count is only for the splice graph view. Not to be used with anything containing quantified data.
        :param args: arg parse
        :return: number of pages
        """
        gene_count = len(tuple(self.get_gene_ids(args)))
        return math.ceil(gene_count / constants.MAX_GENES)
