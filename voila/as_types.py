from bisect import bisect_left, bisect_right
from itertools import combinations
from pathlib import Path

from voila import constants
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_utils import generate_means
import argparse


PSI_THRESHOLD = 0.01
DPSI_THRESHOLD = None


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
        for edge in self.edges:
            if not any(e.start < edge.end < e.end or (e.start > edge.start and e.end == edge.end) for e in self.edges):
                i = bisect_left(self.nodes, edge.node)
                yield self.Module(self.nodes[start_idx: i + 1])
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

                for (start, end), means in zip(lsv.junctions, lsv.get('means')):

                    key = str(start) + '-' + str(end)

                    if key not in lsv_store:
                        lsv_store[key] = {}

                    if lsv_id not in lsv_store[key]:
                        lsv_store[key][lsv_id] = {'psi': [], 'delta_psi': []}

                    lsv_store[key][lsv_id]['psi'].append(means)

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
        def __init__(self, nodes):
            """
            Module is subset of a gene.  The divide between modules is where junctions don't cross.
            :param nodes: list of nodes that belong to module
            """

            self.nodes = nodes  # subset of nodes for this module

        def alternate_downstream(self):
            """
            Check if alternate downstream occurs in this module.
            :return: boolean
            """

            for node in self.nodes[:-1]:
                for e1, e2 in combinations(node.edges, 2):
                    if e1.start == e2.start and e1.end != e2.end:
                        return True

        def alternate_upstream(self):
            """
            Check if alternate upstream occurs in this module.
            :return: boolean
            """

            for node in self.nodes[:-1]:
                for e1, e2 in combinations(node.edges, 2):
                    if e1.end == e2.end and e1.start != e2.start:
                        return True

        def exon_skipping(self):
            """
            Check if exon skipping occurs in this module.
            :return: boolean
            """

            b = self.Filters.target_source_psi
            s = self.Filters.source_psi
            t = self.Filters.target_psi

            for n1, n2, n3 in combinations(self.nodes, 3):
                if n1.connects(n2, s) and n1.connects(n3, b) and n2.connects(n3, t):
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

        def as_types(self):
            """
            Helper function that returns a list of types found in this module.
            :return: list of AS types
            """

            as_type_dict = {
                'alt_downstream': self.alternate_downstream,
                'alt_upstream': self.alternate_upstream,
                'exon_skipping': self.exon_skipping,
                'mutually_exclusive': self.mutually_exclusive
            }
            return [k for k, v in as_type_dict.items() if v()]

        class Filters:
            """
            Class to act as Namespace to hold filter methods.
            """

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

            @staticmethod
            def target_psi(edges):
                """
                Returns edge if it contains junctions that have target psi values that pass threshold.
                :param edges: list of edges
                :return: filtered list of edges
                """

                for edge in edges:
                    target_ids = (l for l in edge.lsvs if ':t:' in l)
                    if any(p >= PSI_THRESHOLD for l in target_ids for p in edge.lsvs[l]['psi']):
                        yield edge

            @staticmethod
            def source_psi(edges):
                """
                Returns edge if it contains junctions that have source psi values that pass threshold.
                :param edges: list of edges
                :return: filtered list of edges
                """

                for edge in edges:
                    source_ids = (l for l in edge.lsvs if ':s:' in l)
                    if any(p >= PSI_THRESHOLD for l in source_ids for p in edge.lsvs[l]['psi']):
                        yield edge


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('sg_file', help='Splicegraph file that comes from the build execution')
    parser.add_argument('voila_file', help='voila file')


    args = parser.parse_args()

    # Find all gene ids in voila file
    with Matrix(Path(args.voila_file).expanduser()) as m:
        gene_ids = list(m.gene_ids)

    # for gene_id in gene_ids:
    # gene_id = 'ENSMUSG00000001419'
    for gene_id in gene_ids:
        print(gene_id)
        graph = Graph(gene_id, args.sg_file, args.voila_file)
        for module in graph.modules():

            # t = timeit.Timer(module.as_types)
            # print(t.timeit(100), module.as_types())
            print(module)
