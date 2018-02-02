import math
from itertools import zip_longest

from voila import constants
from voila.api import SpliceGraph
from voila.api import splice_graph_model as model
from voila.api.view_matrix import ViewDeltaPsi, ViewPsi
from voila.utils.exceptions import GeneIdNotFoundInVoilaFile
from voila.utils.voila_log import voila_log


class ViewSpliceGraph(SpliceGraph):
    def get_page_count(self, args):
        gene_count = self.session.query(model.Gene).count()
        return int(math.ceil(gene_count / float(constants.MAX_GENES)))

    def get_paginated_genes(self, args):
        def grouper(iterable, n, fillvalue=None):
            args = [iter(iterable)] * n
            return zip_longest(*args, fillvalue=fillvalue)

        for page in grouper(self.get_gene_ids(args), constants.MAX_GENES):
            yield tuple(p for p in page if p is not None)

    def get_gene_ids(self, args=None):
        if args and args.gene_ids:
            return args.gene_ids

        if args and hasattr(args, 'lsv_ids') and args.lsv_ids:
            return (lsv_id.split(':')[0] for lsv_id in args.lsv_ids)

        for gene_id, in self.session.query(model.Gene.id).all():
            yield gene_id


class HtmlSpliceGraph(ViewSpliceGraph):
    def combine_genes(self, experiments, genes):
        gene_dict = {}

        for gene_id in genes:
            gene = self.gene(gene_id).get

            comb_gene = None

            for exp in experiments:
                if comb_gene is None:
                    comb_gene = gene.get_experiment(exp)
                else:
                    new_gene = gene.get_experiment(exp)
                    for x, y in zip(comb_gene['junctions'], new_gene['junctions']):
                        x['reads'] += y['reads']

            gene_dict[gene.id] = comb_gene

        return gene_dict


class PsiSpliceGraph(ViewSpliceGraph):
    def get_page_count(self, args=None):
        gene_count = self.session.query(model.Gene).count()
        return int(math.ceil(gene_count / float(constants.MAX_GENES)))

    def get_paginated_genes_with_lsvs(self, args):
        log = voila_log()
        log.debug('Getting paginated genes with LSVs')

        gene_list = []
        lsv_dict = {}

        with ViewPsi(args.voila_file) as m:
            for gene_id in self.get_gene_ids(args):
                lsvs = tuple(m.get_lsvs(args, gene_id=gene_id))

                if lsvs:
                    lsv_dict[gene_id] = tuple(dict(m.psi(lsv_id).get()) for lsv_id in lsvs)
                    gene_list.append(gene_id)

                if len(gene_list) == constants.MAX_GENES:
                    yield lsv_dict, gene_list
                    gene_list = []
                    lsv_dict = {}

            if gene_list:
                yield lsv_dict, gene_list


class DeltaPsiSpliceGraph(PsiSpliceGraph):
    def get_paginated_genes_with_lsvs(self, args):
        log = voila_log()
        log.debug('Getting paginated genes with LSVs')

        gene_list = []
        lsv_dict = {}

        with ViewDeltaPsi(args.voila_file) as m:
            for gene_id in self.get_gene_ids(args):
                try:
                    lsvs = tuple(m.get_lsv_ids(args, gene_id=gene_id))
                except GeneIdNotFoundInVoilaFile:
                    lsvs = None

                if lsvs:
                    lsv_dict[gene_id] = tuple(dict(m.delta_psi(lsv_id).get()) for lsv_id in lsvs)
                    gene_list.append(gene_id)

                if len(gene_list) == constants.MAX_GENES:
                    yield lsv_dict, gene_list
                    gene_list = []
                    lsv_dict = {}

            if gene_list:
                yield lsv_dict, gene_list
