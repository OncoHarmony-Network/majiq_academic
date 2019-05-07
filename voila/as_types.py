from bisect import bisect_left, bisect_right
from itertools import combinations
from pathlib import Path

from voila import constants
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_utils import generate_means

import csv

import argparse

def check_file(value):
    """
    Check if file exists.
    :param value: file path
    :return:
    """
    value = Path(value)

    value = value.expanduser()
    value = value.absolute()

    if value.exists():
        return value
    else:
        raise Exception("Can't find file %s" % value)

parser = argparse.ArgumentParser(description='LSV Classification Script')
parser.add_argument('psi_file', type=check_file,
                        help='File path of input voila psi file')
parser.add_argument('splicegraph_file', type=check_file,
                        help='File path of input splicegraph file')
parser.add_argument('-f', '--file-name', required=True, help="Output TSV file path")
parser.add_argument('--hide-sub-complex', action='store_true',
                    help='If a module is detected as complex, do not show any other events detected in that module')
parser.add_argument('--psi-threshold', type=float, default=0.0,
                    help='Only detect junctions where target and source psi values pass threshold. '
                         '(0.0, the default, accepts everything)')

args = parser.parse_args()


PSI_THRESHOLD = args.psi_threshold
DPSI_THRESHOLD = None
HIDE_SUB_COMPLEX = args.hide_sub_complex


class UnsupportedVoilaFile(Exception):
    pass


class Graph:
    def __init__(self, gene_id, sg_file, voila_file):
        """
        This contains the edges and nodes used to find modules.

        :param gene_id: gene id
        :param sg_file: splice graph file name
        :param voila_file: voila file name (psi/delta psi)
        """

        self.nodes = []  # all nodes in the graph
        self.edges = []  # all edges in the graph
        self.gene_id = gene_id
        self.sg_file = Path(sg_file).expanduser().resolve()
        self.voila_file = Path(voila_file).expanduser().resolve()

        # populate the graph with data from the splice graph
        self._populate()

        # get psi/dpsi data from voila file and associate with edges
        with Matrix(self.voila_file) as m:
            analysis_type = m.analysis_type

        with SpliceGraph(self.sg_file) as sg:
            self.strand = sg.gene(self.gene_id)['strand']

        if analysis_type == constants.ANALYSIS_PSI:
            self._psi()
        elif analysis_type == constants.ANALYSIS_DELTAPSI:
            self._delta_psi()
        else:
            raise UnsupportedVoilaFile()

        # find connections between nodes
        self._find_connections()

    class Node:
        def __init__(self, exon):
            """
            Graph wrapper for exons.

            :param exon: exon dictionary from splice graph file.
            """

            self.edges = []  # all edge starting in this exon
            self.exon = exon  # exon dictionary

        def __eq__(self, other):
            """
            Exons are uniquely defined by gene id, start, and end.  Since Graphs work on one gene at a time, equality
            is determined by start and end.
            :param other: other exon
            :return: boolean
            """

            return self.start == other.start and self.end == other.end

        def __lt__(self, other):
            """
            To better sort exons, we use visual start and ends of exons.

            :param other: other exon
            :return: boolean
            """

            return self.view_start < other.view_start and self.view_end < other.view_start

        def __repr__(self):
            """
            A string representation of this exon including it's start and end.
            :return: string
            """

            return '<{} {},{}>'.format(self.__class__.__name__, self.start, self.end)

        @property
        def start(self):
            """
            Start of exon.
            :return: integer
            """

            return self.exon['start']

        @property
        def end(self):
            """
            End of exon.
            :return: integer
            """

            return self.exon['end']

        @property
        def view_start(self):
            """
            Alter start of exon if exon has no start. Visual start of exon.
            :return: integer
            """

            start = self.exon['start']
            if start == -1:
                start = self.end - 10
            return start

        @property
        def view_end(self):
            """
            Alter end of exon if exon has no edn. Visual end of exon.
            :return: integer
            """

            end = self.exon['end']
            if end == -1:
                end = self.start + 10
            return end

        def connects(self, node, filter=None):
            """
            Search through junctions for this exon to see if this exon has a junction that connects to supplied exon.
            :param node: the exon that this exon might connect to.
            :param filter: function to filter junctions
            :return: boolean
            """

            edges = self.edges
            if filter:
                edges = filter(edges)

            # print(node)
            # print([x.node for x in edges])
            # print([x.node for x in edges])

            return any(e.node == node for e in edges)

    class Edge:
        def __init__(self, junc):
            """
            Graph wrapper for junctions.
            :param junc: junction dictionary from the splice graph file.
            """

            self.junc = junc  # junction dictionary
            self.node = None  # node this junction connects to
            self.lsvs = {}

        def __lt__(self, other):
            """
            Junction are uniquely identified by gene_id, start, end.  Since graphs work on one gene at a time, we order
            junction by start and end.
            :param other: other junction
            :return: boolean
            """

            return self.start < other.start and self.end < other.start

        def __eq__(self, other):
            """
            Equality is determined by start and end of junction.
            :param other:
            :return:
            """

            return self.start == other.start and self.end == other.end

        def __repr__(self):
            """
            String representation of junction with start and end.
            :return: string
            """

            return '<{} {},{}>'.format(self.__class__.__name__, self.start, self.end)

        @property
        def start(self):
            """
            Junction start.
            :return: integer
            """

            return self.junc['start']

        @property
        def end(self):
            """
            Junction end.
            :return: interger
            """

            return self.junc['end']

        @property
        def view_start(self):
            """
            For compatibility with using bisect to find which exon this junction starts/stops in.
            :return: integer
            """

            return self.start

        @property
        def view_end(self):
            """
            For compatibility with using bisect to find which exon this junction starts/stops in.
            :return: integer
            """

            return self.end

    def start_node(self, edge):
        """
        Get exon where this junction starts.
        :param edge: supplied junction
        :return: node object
        """

        i = bisect_left(self.nodes, edge)
        return self.nodes[i]

    def end_node(self, edge):
        """
        Get exon where this junction ends.
        :param edge: supplied junction
        :return: node object
        """

        i = bisect_right(self.nodes, edge)
        return self.nodes[i - 1]

    def _add_junc(self, junc):
        """
        Add junction to graph as edge object. If junction starts and ends in the same exon, it is not added to graph.
        :param junc: junction dictionary from splice graph.
        :return: None
        """

        edge = self.Edge(junc)
        start_node = self.start_node(edge)
        end_node = self.end_node(edge)

        # Since majiq doesn't quantify junctions that start/stop in same exon, filter them.
        if start_node != end_node:
            self.edges.append(edge)
            edge.node = end_node

    def _add_exon(self, exon):
        """
        Added exon to graph as node object.
        :param exon: exon dictionary from splice graph
        :return: None
        """

        self.nodes.append(self.Node(exon))

    def _find_connections(self):
        """
        When this has completed, each exon should have a list of exons that start there.
        :return: None
        """

        for edge in self.edges:
            node = self.start_node(edge)
            node.edges.append(edge)

    def _populate(self):
        """
        Add all juctions and exons to graph and sort those lists.
        :return: None
        """

        with SpliceGraph(self.sg_file) as sg:
            for exon in sg.exons(self.gene_id):
                self._add_exon(exon)
            for junc in sg.junctions(self.gene_id):
                self._add_junc(junc)

        self.edges.sort()
        self.nodes.sort()

    def modules(self):
        """
        Search through edges to find where they don't cross.  At this point is where the previous module ends and the
        next module starts. There will be an over lapping exon.
        :return: Generator of modules
        """

        start_idx = 0
        j = 0
        for edge in self.edges:


            if not any(e.start < edge.end < e.end or (e.start > edge.start and e.end == edge.end) for e in self.edges):
                j += 1
                i = bisect_left(self.nodes, edge.node)
                if j == 2 or True:
                    yield self.Module(self.nodes[start_idx: i + 1], self.strand)
                start_idx = i

    def _psi(self):
        """
        When a psi voila file is supplied, this is where the psi data is added to the junctions.
        :return: None
        """

        lsv_store = {}

        with Matrix(self.voila_file) as m:
            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):
                lsv = m.psi(lsv_id)
                experiment_names = [m.experiment_names[0][0]]

                for (start, end), means in zip(lsv.junctions, lsv.get('means')):

                    key = str(start) + '-' + str(end)

                    if key not in lsv_store:
                        lsv_store[key] = {}

                    if lsv_id not in lsv_store[key]:
                        lsv_store[key][lsv_id] = {'psi': [], 'delta_psi': []}

                    lsv_store[key][lsv_id]['psi'].append(means)

                    if lsv.intron_retention:
                        with SpliceGraph(self.sg_file) as sg:
                            if [x for x in sg.intron_retention_reads_exp({'gene_id': self.gene_id,
                                                              'start': int(start)+1,
                                                              'end': int(end)-1}, experiment_names)]:
                                lsv_store[key]['IR'] = True


                    if lsv.a3ss:
                        lsv_store[key]['a3ss'] = True


                    if lsv.a5ss:
                        lsv_store[key]['a5ss'] = True



        for edge in self.edges:
            key = str(edge.start) + '-' + str(edge.end)
            if key in lsv_store:
                edge.lsvs = lsv_store[key]

    def _delta_psi(self):
        """
        When a delta psi voila file is supplied, this is where the psi/delta psi data is added to the junctions.
        :return: None
        """

        lsv_store = {}

        with Matrix(self.voila_file) as m:

            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):

                lsv = m.delta_psi(lsv_id)
                juncs = lsv.junctions

                for group_means in lsv.get('group_means'):
                    for (start, end), means in zip(juncs, group_means):
                        key = str(start) + '-' + str(end)

                        if key not in lsv_store:
                            lsv_store[key] = {}

                        if lsv_id not in lsv_store[key]:
                            lsv_store[key][lsv_id] = {'psi': [], 'delta_psi': []}

                        lsv_store[key][lsv_id]['psi'].append(means)

                for (start, end), means in zip(juncs, generate_means(lsv.get('bins'))):
                    key = str(start) + '-' + str(end)
                    lsv_store[key][lsv_id]['delta_psi'].append(means)

        for edge in self.edges:
            key = str(edge.start) + '-' + str(edge.end)
            edge.lsvs = lsv_store[key]

    class Module:
        def __init__(self, nodes, strand):
            """
            Module is subset of a gene.  The divide between modules is where junctions don't cross.
            :param nodes: list of nodes that belong to module
            """

            self.nodes = nodes  # subset of nodes for this module
            self.Filters.strand = strand

        def alternate_downstream(self):
            """
            Check if alternate downstream occurs in this module.
            :return: boolean
            """

            # for node in self.nodes[:-1]:
            #     for e1, e2 in combinations(node.edges, 2):
            #         if e1.start == e2.start and e1.end != e2.end:
            #             return True

        def alternate_upstream(self):
            """
            Check if alternate upstream occurs in this module.
            :return: boolean
            """

            # for node in self.nodes[:-1]:
            #     for e1, e2 in combinations(node.edges, 2):
            #         if e1.end == e2.end and e1.start != e2.start:
            #             return True

        def exon_skipping(self):
            """
            Check if exon skipping occurs in this module.

            """

            b = self.Filters.target_source_psi
            s = self.Filters.source_psi
            t = self.Filters.target_psi

            num_found = 0
            for n1, n2, n3 in combinations(self.nodes, 3):
                # print(n1, n2, n3)
                # print('--------------------')
                # print(n1.connects(n2, t))
                # print('--------------------')
                # print(n1.connects(n3, b))
                # print('--------------------')
                # print(n2.connects(n3, s))
                if n1.connects(n2, s) and n1.connects(n3, b) and n2.connects(n3, t):
                    #return True
                    num_found += 1
            return num_found

        def multi_exon_skipping(self):
            """
            Check if multi exon skipping occurs in this module.
            :return: boolean
            """

            if len(self.nodes) < 4:
                return False

            b = self.Filters.target_source_psi

            # exclude half exons
            full_exons = list(filter(lambda ex: ex.start != -1 and ex.end != -1, self.nodes))

            if len(full_exons) < 4:
                return False

            # find combinations of nodes with at least two nodes in between
            # we assume nodes are ordered by default
            # we don't care to find different ordering either
            for i, n1 in enumerate(full_exons):
                for j, n2 in enumerate(full_exons):
                    if j - i > 2 and n1.connects(n2, b):
                        return True



        def mutually_exclusive(self):
            """
            Check if mutually exclusive occurs in this module.
            :return: boolean
            """

            f = self.Filters.target_source_psi

            # for n1, n2, n3, n4 in combinations(self.nodes, 4):
            #     if n1.connects(n2, f) and n1.connects(n3, f) and n2.connects(n4, f) and n3.connects(n4, f):
            #         if not n2.connects(n3):
            #             return True

            for n1 in self.nodes[:-1]:
                for e1, e2, in combinations(f(n1.edges), 2):
                    if e1.node != e2.node and not (e1.node.connects(e2.node, f) or e2.node.connects(e1.node, f)):
                        for i1 in f(e1.node.edges):
                            for i2 in f(e2.node.edges):
                                if i1.node == i2.node:
                                    return True

        def intron_retention(self):
            """
            Check if intron retention occurs in this module.
            :return: boolean
            """

            b = self.Filters.target_source_psi

            for n1, n2 in combinations(self.nodes, 2):
                for edge in n1.edges:
                    if 'IR' in edge.lsvs and n1.connects(n2, b):
                        return True
                for edge in n2.edges:
                    if 'IR' in edge.lsvs and n2.connects(n1, b):
                        return True


        def altTranscStart(self):
            for node in self.nodes:
                if node.start == -1:
                    return True

        def altTranscEnd(self):
            for node in self.nodes:
                if node.end == -1:
                    return True

        def alt5ss(self):


            b = self.Filters.target_source_psi

            # print(self.nodes)
            # for node in self.nodes:
            #
            #     for e1, e2 in combinations(node.edges, 2):
            #         print('h')
            #         print(e1.start == e2.start and e1.end != e2.end)
            #         if e1.start == e2.start and e1.end != e2.end:
            #             return True

            for n1, n2 in combinations(self.nodes, 2):
                for edge in n1.edges:
                    if 'a5ss' in edge.lsvs and n1.connects(n2, b):
                        return True
                for edge in n2.edges:
                    if 'a5ss' in edge.lsvs and n2.connects(n1, b):
                        return True

        def alt3ss(self):

            b = self.Filters.target_source_psi

            for n1, n2 in combinations(self.nodes, 2):
                for edge in n1.edges:
                    if 'a3ss' in edge.lsvs and n1.connects(n2, b):
                        return True
                for edge in n2.edges:
                    if 'a3ss' in edge.lsvs and n2.connects(n1, b):
                        return True

        def as_types(self):
            """
            Helper function that returns a list of types found in this module.
            :return: list of AS types
            """

            as_type_dict = {
                # 'alt_downstream': self.alternate_downstream,
                # 'alt_upstream': self.alternate_upstream,
                'cassette_exon': self.exon_skipping,
                'mutually_exclusive': self.mutually_exclusive,
                'intron_retention': self.intron_retention,
                'alt3ss': self.alt3ss,
                'alt5ss': self.alt5ss,
                'altTranscStart': self.altTranscStart,
                'altTranscEnd': self.altTranscEnd,
                'multi_exon_skipping': self.multi_exon_skipping
            }
            ret = []
            for k, v in as_type_dict.items():
                res = v()
                if k == 'cassette_exon' and res and res > 1:
                    ret.append('complex')
                if res:
                    ret.append(k)
            if not 'complex' in ret and len(ret) > 1:
                ret.append('complex')
            if HIDE_SUB_COMPLEX and 'complex' in ret and len(ret) > 1:
                ret = ['complex']
            return ret

        class Filters:
            """
            Class to act as Namespace to hold filter methods.
            """
            rev_map = {':t:': ':s:', ':s:': ':t:'}

            @classmethod
            def strand_select(cls, lsv_str):
                """
                For some reason, the way it is designed, the sources and targets are reversed for - strands compared
                to + strands. I am guessing this may eventually be more complicated, where I would not be able to
                do a simple reversal like here, but currently it works correctly.
                :param lsv_str:
                :return:
                """
                if cls.strand == '+':
                    return lsv_str
                else:
                    return cls.rev_map[lsv_str]

            @classmethod
            def target_source_psi(cls, edges):
                """
                Returns edge if it contains junctions that have target and source psi values that pass threshold.
                :param edges: list of edges
                :return: filtered list of edges
                """

                if PSI_THRESHOLD:
                    return cls.target_psi(cls.source_psi(edges))
                else:
                    return edges

            @classmethod
            def target_psi(cls, edges):
                """
                Returns edge if it contains junctions that have target psi values that pass threshold.
                :param edges: list of edges
                :return: filtered list of edges
                """

                for edge in edges:
                    target_ids = (l for l in edge.lsvs if cls.strand_select(':t:') in l)
                    #print([x for x in target_ids])
                    if any(p >= PSI_THRESHOLD for l in target_ids for p in edge.lsvs[l]['psi']):
                        yield edge

            @classmethod
            def source_psi(cls, edges):
                """
                Returns edge if it contains junctions that have source psi values that pass threshold.
                :param edges: list of edges
                :return: filtered list of edges
                """

                for edge in edges:
                    source_ids = (l for l in edge.lsvs if cls.strand_select(':s:') in l)
                    if any(p >= PSI_THRESHOLD for l in source_ids for p in edge.lsvs[l]['psi']):
                        yield edge


def output_tsv(genes_modules):
    """
    Write the output file
    :param genes_modules: a list of (gene_id (str), gene_modules (obj)) tuples
    :return: NOTHING
    """
    with open(args.file_name, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

        writer.writerow(['gene_id', 'module_idx', 'complex', 'cassette_exon', 'mutually_exclusive',
                    'intron_retention', 'alt3ss', 'alt5ss', 'altTranscStart', 'altTranscEnd',
                    'multi_exon_skipping'])

        for gene_id, modules in genes_modules:
            for i, module in enumerate(modules, 1):
                types = module.as_types()
                writer.writerow((gene_id, i,
                                 1 if 'complex' in types else 0,
                                 1 if 'cassette_exon' in types else 0,
                                 1 if 'mutually_exclusive' in types else 0,
                                 1 if 'intron_retention' in types else 0,
                                 1 if 'alt3ss' in types else 0,
                                 1 if 'alt5ss' in types else 0,
                                 1 if 'altTranscStart' in types else 0,
                                 1 if 'altTranscEnd' in types else 0,
                                 1 if 'multi_exon_skipping' in types else 0,)
                                )


if __name__ == "__main__":




    sg_file = '/home/paul/PycharmProjects/majiq/test_cases/VOILA_DEV/ORIGINAL_BUILD/splicegraph.sql'
    psi_file = '/home/paul/PycharmProjects/majiq/test_cases/VOILA_DEV/ORIGINAL_PSI/test.psi.voila'
    #dpsi_file = '~/Development/small_test/majiq_deltapsi_all_v_all/Adr_Cer.deltapsi.voila'

    psi_file = args.psi_file
    sg_file = args.splicegraph_file


    # Find all gene ids in splice graph
    # with SpliceGraph(sg_file) as sg:
    #    gene_ids = list(g['id'] for g in sg.genes())

    # Find all gene ids in voila file
    with Matrix(Path(psi_file).expanduser()) as m:
        gene_ids = list(m.gene_ids)

    # for gene_id in gene_ids:
    # gene_id = 'ENSMUSG00000001419'
    genes_modules = []
    for gene_id in gene_ids:
        #if gene_id == "ENSMUSG00000021820":
        #if gene_id == "ENSMUSG00000026843":
        #if gene_id == "ENSMUSG00000024097":
        #if gene_id == "ENSMUSG00000006498":
        #if gene_id == "ENSMUSG00000027287":
        #if gene_id == "ENSMUSG00000001419":
        #     print(gene_id)
        #     graph = Graph(gene_id, sg_file, psi_file)
        #
        #     mod = [x for x in graph.modules()][-3]
        #     print(mod.nodes)
        #     print(mod.as_types())

        graph = Graph(gene_id, sg_file, psi_file)
        genes_modules.append((gene_id, graph.modules()))
        #for module in graph.modules():
            # t = timeit.Timer(module.as_types)
            # print(t.timeit(100), module.as_types())

            #print(module.as_types())


    output_tsv(genes_modules)

    # what to do with exon 19-20 event in http://localhost:5005/gene/ENSMUSG00000021820/


    ### TESTS ###

    # for gene_id in gene_ids:
    #     if gene_id == "ENSMUSG00000021820":
    #
    #         expected = {5: ['cassette_exon'], 14: ['cassette_exon']}
    #         graph = Graph(gene_id, sg_file, psi_file)
    #
    #         for i, module in enumerate(graph.modules()):
    #
    #             if i in expected:
    #                 assert module.as_types() == expected[i]
    #
    # print("Succcess!")