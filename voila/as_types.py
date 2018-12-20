from bisect import bisect_left, bisect_right
from itertools import combinations, permutations
from pathlib import Path

from voila import constants
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_utils import generate_means

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

            self.source_psi = []  # list of source psi values
            self.target_psi = []  # list of target psi values
            self.delta_psi = []  # list of dpsi values
            self.lsv_ids = []  # list of lsv ids

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
        target = {}
        source = {}

        with Matrix(self.voila_file) as m:
            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):
                lsv = m.psi(lsv_id)
                for (start, end), means in zip(lsv.junctions, lsv.get('means')):
                    key = str(start) + '-' + str(end)
                    if lsv_id.split(':')[-2] == 's':
                        source[key] = means
                    else:
                        target[key] = means

        for edge in self.edges:
            key = str(edge.start) + '-' + str(edge.end)
            if key in source:
                edge.source_psi.append(source[key])
            if key in target:
                edge.target_psi.append(target[key])

    def _delta_psi(self):
        """
        When a delta psi voila file is supplied, this is where the psi/delta psi data is added to the junctions.
        :return: None
        """
        target_psi = {}
        source_psi = {}
        junc_lsv_ids = {}
        dpsi = {}

        with Matrix(self.voila_file) as m:
            group_names = m.group_names
            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):

                lsv = m.delta_psi(lsv_id)
                juncs = lsv.junctions

                for group, group_means in zip(group_names, lsv.get('group_means')):

                    for (start, end), means in zip(juncs, group_means):
                        key = str(start) + '-' + str(end)

                        if lsv_id.split(':')[-2] == 's':
                            psi_store = source_psi
                        else:
                            psi_store = target_psi

                        if key in psi_store:
                            psi_store[key].append(means)
                        else:
                            psi_store[key] = [means]

                        if key in junc_lsv_ids:
                            junc_lsv_ids[key].add(lsv_id)
                        else:
                            junc_lsv_ids[key] = {lsv_id}

                for (start, end), means in zip(juncs, generate_means(lsv.get('bins'))):
                    key = str(start) + '-' + str(end)
                    dpsi[key] = means

        for edge in self.edges:
            key = str(edge.start) + '-' + str(edge.end)
            if key in source_psi:
                edge.source_psi += source_psi[key]
            if key in target_psi:
                edge.target_psi += target_psi[key]
            if key in dpsi:
                edge._delta_psi = dpsi[key]
            if key in junc_lsv_ids:
                edge.lsv_ids += list(junc_lsv_ids[key])

    class Module:
        def __init__(self, nodes):
            self.nodes = nodes  # subset of nodes for this module

        def starts_equal(self, e1, e2):
            return e1.start == e2.start

        def ends_equal(self, e1, e2):
            return e1.end == e2.end

        def alternate_downstream(self):
            for node in self.nodes[:-1]:
                for e1, e2 in combinations(node.edges, 2):
                    if self.starts_equal(e1, e2) and not self.ends_equal(e1, e2):
                        return True

        def alternate_upstream(self):
            for node in self.nodes[:-1]:
                for e1, e2 in combinations(node.edges, 2):
                    if self.ends_equal(e1, e2) and not self.starts_equal(e1, e2):
                        return True

        def exon_skipping(self):
            b = self.Filters.target_source_psi
            s = self.Filters.source_psi
            t = self.Filters.target_psi

            for n1, n2, n3 in permutations(self.nodes, 3):
                if n1.connects(n2, s) and n1.connects(n3, b) and n2.connects(n3, t):
                    return True

        def mutually_exclusive(self):
            f = self.Filters.target_source_psi

            for n1, n2, n3, n4 in permutations(self.nodes, 4):
                if n1.connects(n2, f) and n1.connects(n3, f) and n2.connects(n4, f) and n3.connects(n4, f):
                    if not n2.connects(n3):
                        return True

        def as_types(self):
            as_type_dict = {
                'alt_downstream': self.alternate_downstream,
                'alt_upstream': self.alternate_upstream,
                'exon_skipping': self.exon_skipping,
                'mutually_exclusive': self.mutually_exclusive
            }
            return [k for k, v in as_type_dict.items() if v()]

        class Filters:

            @classmethod
            def target_source_psi(cls, edges):
                if PSI_THRESHOLD:
                    return cls.target_psi(cls.source_psi(edges))
                else:
                    return edges

            @staticmethod
            def target_psi(edges):
                if PSI_THRESHOLD:
                    return list(e for e in edges if any(p >= PSI_THRESHOLD for p in e.target_psi))
                else:
                    return edges

            @staticmethod
            def source_psi(edges):
                if PSI_THRESHOLD:
                    return list(e for e in edges if any(p >= PSI_THRESHOLD for p in e.source_psi))
                else:
                    return edges


if __name__ == "__main__":
    sg_file = '~/Development/small_test/majiq_build/splicegraph.sql'
    psi_file = '~/Development/small_test/majiq_psi_all/Adr.psi.voila'
    dpsi_file = '~/Development/small_test/majiq_deltapsi_all_v_all/Adr_Cer.deltapsi.voila'

    # for gene_id in gene_ids:
    gene_id = 'ENSMUSG00000001419'
    graph = Graph(gene_id, sg_file, psi_file)

    for module in graph.modules():
        for as_type in module.as_types():
            print(as_type)
