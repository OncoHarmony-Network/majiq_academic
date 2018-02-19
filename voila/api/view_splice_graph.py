import math
from itertools import zip_longest

from sqlalchemy import inspect

from voila import constants
from voila.api import SpliceGraph
from voila.api import splice_graph_model as model
from voila.api.matrix_hdf5 import lsv_id_to_gene_id
from voila.api.splice_graph_model import Reads


class ViewSpliceGraph(SpliceGraph):
    def get_paginated_genes(self, args):
        def grouper(iterable, n, fillvalue=None):
            args = [iter(iterable)] * n
            return zip_longest(*args, fillvalue=fillvalue)

        for page in grouper(self.get_genes(args), constants.MAX_GENES):
            yield tuple(p for p in page if p is not None)

    def get_genes(self, args):
        if args.gene_ids:
            for gene_id in args.gene_ids:
                yield self.session.query(model.Gene).get(gene_id)
        elif hasattr(args, 'lsv_ids') and args.lsv_ids:
            for lsv_id in args.lsv_ids:
                yield self.session.query(model.Gene).get(lsv_id_to_gene_id(lsv_id))
        else:
            yield from self.session.query(model.Gene).all()


class SpliceGraphViewSpliceGraph(ViewSpliceGraph):
    def get_page_count(self, args):
        """
        This page count is only for the splice graph view. Not to be used with anything containing quantified data.
        :param args: arg parse
        :return: number of pages
        """
        gene_count = len(tuple(self.get_genes(args)))
        return math.ceil(gene_count / constants.MAX_GENES)


class ViewJunction:
    def __init__(self, junction):
        self.junction = junction

    def junction_types(self, experiment_names):
        for experiment_name in experiment_names:
            reads = self.reads(experiment_name)

            junc_type = constants.JUNCTION_TYPE_RNASEQ

            if reads == 0 and self.junction.annotated:
                if self.reads_sum(experiment_names) - reads > 0:
                    junc_type = constants.JUNCTION_TYPE_DB_OTHER_RNASEQ
                else:
                    junc_type = constants.JUNCTION_TYPE_DB

            if self.junction.annotated and reads > 0:
                junc_type = constants.JUNCTION_TYPE_DB_RNASEQ

            if not self.junction.annotated and reads > 0:
                junc_type = constants.JUNCTION_TYPE_RNASEQ

            yield experiment_name, junc_type

    def reads_sum(self, experiment_names):
        return sum(self.reads(e) for e in experiment_names)

    def reads(self, experiment_name):
        session = inspect(self.junction).session
        r = session.query(Reads).get((experiment_name, self.junction.gene_id, self.junction.start, self.junction.end))
        try:
            return r.reads
        except AttributeError:
            return 0

    def get_experiment(self):
        j = dict(self.junction)
        j['intron_retention'] = self.get_intron_retention_type()
        del j['gene_id']
        return j

    def get_intron_retention_type(self):
        if self.junction.intron_retention:
            for exon in filter(lambda e: not e.intron_retention, self.junction.gene.exons):
                if self.junction.start in exon:
                    return constants.IR_TYPE_START
                elif self.junction.end in exon:
                    return constants.IR_TYPE_END

        return self.junction.intron_retention


class ViewExon:
    def __init__(self, exon):
        self.exon = exon

    @property
    def view_start(self):
        if self.exon.start == -1:
            return self.exon.end - 10
        return self.exon.start

    @property
    def view_end(self):
        if self.exon.end == -1:
            return self.exon.start + 10
        return self.exon.end

    def has_reads(self, experiment_name):
        return any(ViewJunction(j).reads(experiment_name) > 0 for j in self.exon.a5) or any(
            ViewJunction(j).reads(experiment_name) > 0 for j in self.exon.a3)

    def get_exon_type(self, experiment_name):
        if self.exon.start == -1:
            return constants.EXON_TYPE_MISSING_START
        if self.exon.end == -1:
            return constants.EXON_TYPE_MISSING_END

        has_reads = self.has_reads(experiment_name)

        if self.exon.annotated and not has_reads:
            return constants.EXON_TYPE_DB

        if self.exon.annotated and has_reads:
            return constants.EXON_TYPE_DB_RNASEQ

        if not self.exon.annotated and has_reads:
            return constants.EXON_TYPE_RNASEQ

        return constants.EXON_TYPE_RNASEQ

    def get_experiment(self):
        exon = dict(self.exon)
        exon['start'] = self.view_start
        exon['end'] = self.view_end
        del exon['gene_id']
        return exon


class ViewGene:
    def __init__(self, gene):
        self.gene = gene

    def lsv_reference_exon(self, lsv_id):
        lsv_coords = lsv_id.split(':')[-1]
        if lsv_coords.startswith('-1'):
            coords_list = [-1, int(lsv_coords.split('-')[-1])]
        elif lsv_coords.endswith('-1'):
            coords_list = [int(lsv_coords.split('-')[0]), -1]
        else:
            coords_list = list(map(int, lsv_coords.split('-')))

        for exon in self.gene.exons:
            if coords_list == [exon.start, exon.end]:
                return exon

    def lsv_exons(self, lsv, lsv_junctions=None):
        lsv_id = lsv.lsv_id
        if lsv_junctions is None:
            lsv_junctions = self.lsv_junctions(lsv)
        is_target = lsv_id.split(':')[-2] == 't'

        def find_exons():
            yield self.lsv_reference_exon(lsv_id)

            for junc in lsv_junctions:
                for exon in self.gene.exons:
                    if is_target:
                        if self.gene.strand == '+':
                            if junc.start in exon:
                                yield exon
                        else:
                            if junc.end in exon:
                                yield exon
                    else:
                        if self.gene.strand == '+':
                            if junc.end in exon:
                                yield exon
                        else:
                            if junc.start in exon:
                                yield exon

        yield from sorted(find_exons(), key=lambda e: [e.start, e.end])

    def lsv_junctions(self, lsv):
        for start, end in lsv.junctions:
            for junc in self.gene.junctions:
                if [start, end] == [junc.start, junc.end]:
                    yield junc

    def lsv_ucsc_coordinates(self, lsv):
        exons = tuple(self.lsv_exons(lsv))
        start_exon = sorted(e.start for e in exons if e.start != -1)[0]
        end_exon = sorted(e.end for e in exons if e.end != -1)[-1]
        return {'start': start_exon, 'end': end_exon}

    def get_experiment(self, experiment_names_list):
        gene = dict(self.gene)
        gene['_id'] = self.gene.id

        # todo: exons should NOT be sorted.
        exons = tuple(ViewExon(e).get_experiment() for e in self.gene.exons)
        exons = sorted(exons, key=lambda e: (e['start'], e['end']))
        gene['exons'] = exons

        gene['junctions'] = tuple(ViewJunction(j).get_experiment() for j in self.gene.junctions)

        gene['start'] = self.gene.start
        gene['end'] = self.gene.end

        gene['exon_types'] = {}
        gene['junction_types'] = {}
        gene['reads'] = {}

        for experiment_names in experiment_names_list:
            combined_name = next((n for n in experiment_names if ' Combined' in n), '')
            experiment_names = experiment_names[experiment_names != combined_name]

            for name in experiment_names:
                gene['reads'][name] = {}
                gene['junction_types'][name] = {}
                if combined_name:
                    gene['reads'][combined_name] = {}
                    gene['junction_types'][combined_name] = {}

            for exon in self.gene.exons:
                view_exon = ViewExon(exon)
                for experiment_name in experiment_names:
                    try:
                        exon_start = gene['exon_types'][view_exon.view_start]
                    except KeyError:
                        gene['exon_types'][view_exon.view_start] = {}
                        exon_start = gene['exon_types'][view_exon.view_start]

                    try:
                        exon_end = exon_start[view_exon.view_end]
                    except KeyError:
                        exon_start[view_exon.view_end] = {}
                        exon_end = exon_start[view_exon.view_end]

                    exon_end[experiment_name] = view_exon.get_exon_type(experiment_name)

                if combined_name:
                    exon_end[combined_name] = min(exon_end[x] for x in experiment_names)

            for junc in self.gene.junctions:
                view_junc = ViewJunction(junc)

                for experiment_name in experiment_names:
                    try:
                        gene['reads'][experiment_name][junc.start][junc.end] = view_junc.reads(experiment_name)
                    except KeyError:
                        gene['reads'][experiment_name][junc.start] = {junc.end: view_junc.reads(experiment_name)}

                for name, junc_type in view_junc.junction_types(experiment_names):
                    try:
                        gene['junction_types'][name][junc.start][junc.end] = junc_type
                    except KeyError:
                        gene['junction_types'][name][junc.start] = {junc.end: junc_type}

                if combined_name:
                    summed_reads = sum(view_junc.reads(n) for n in experiment_names)
                    try:
                        gene['reads'][combined_name][junc.start][junc.end] = summed_reads
                    except KeyError:
                        gene['reads'][combined_name][junc.start] = {junc.end: summed_reads}

                    combined_type = min(gene['junction_types'][n][junc.start][junc.end] for n in experiment_names)
                    try:
                        gene['junction_types'][combined_name][junc.start][junc.end] = combined_type
                    except KeyError:
                        gene['junction_types'][combined_name][junc.start] = {junc.end: combined_type}

        return gene
