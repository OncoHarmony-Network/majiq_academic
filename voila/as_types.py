import itertools
from bisect import bisect_left
from pathlib import Path

from voila.api import SpliceGraph

sg_file = '~/Development/small_test/majiq_build/splicegraph.sql'
sg_file = Path(sg_file)
sg_file = sg_file.expanduser().resolve()


class Graph:
    def __init__(self, gene_id):
        self.nodes = []
        self.edges = []
        self.gene_id = gene_id

    class Module:
        def __init__(self, nodes):
            self.nodes = nodes

        def starts_equal(self, e1, e2):
            return e1.start == e2.start

        def ends_equal(self, e1, e2):
            return e1.end == e2.end

        def alternate_downstream(self):
            for node in self.nodes:
                for e1, e2 in itertools.combinations(node.edges, 2):
                    assert e1 != e2
                    if self.starts_equal(e1, e2) and not self.ends_equal(e1, e2):
                        print(e1)
                        print(e2)
                        return True

        def alternate_upstream(self):
            for node in self.nodes:
                for e1, e2 in itertools.combinations(node.edges, 2):
                    assert e1 != e2
                    if self.ends_equal(e1, e2) and not self.starts_equal(e1, e2):
                        return True

    class Node:
        def __init__(self, exon):
            self.edges = []
            self.exon = exon

        def __eq__(self, other):
            return self.start == other.start and self.end == other.end

        def __lt__(self, other):
            return self.start < other.start and self.end < other.start

        def __repr__(self):
            return '<{} {},{}>'.format(self.__class__.__name__, self.start, self.end)

        @property
        def start(self):
            return self.exon['start']

        @property
        def end(self):
            return self.exon['end']

    class Edge:
        def __init__(self, junc):
            self.junc = junc
            self.node = None

        def __lt__(self, other):
            return self.start < other.start and self.end < other.start

        def __eq__(self, other):
            return self.start == other.start and self.end == other.end

        def __repr__(self):
            return '<{} {},{}>'.format(self.__class__.__name__, self.start, self.end)

        @property
        def start(self):
            return self.junc['start']

        @property
        def end(self):
            return self.junc['end']

        def starts_in(self, node):
            if node.end == -1:
                return False
            if node.start == -1:
                return self.start == node.end
            else:
                return node.start <= self.start <= node.end

        def ends_in(self, nodes):
            # todo: binary search?
            for node in nodes:
                if (node.end == -1 and node.start == self.end) or (node.start <= self.end <= node.end):
                    self.node = node
                    return

    def add_junc(self, junc):
        self.edges.append(self.Edge(junc))

    def add_exon(self, exon):
        self.nodes.append(self.Node(exon))

    def find_connections(self):
        for node in self.nodes:
            for edge in self.edges:
                if edge.starts_in(node):
                    node.edges.append(edge)
                    edge.ends_in(self.nodes)

    def populate(self):
        with SpliceGraph(sg_file) as sg:
            for exon in sg.exons(self.gene_id):
                self.add_exon(exon)
            for junc in sg.junctions(self.gene_id):
                self.add_junc(junc)

        self.edges.sort()
        self.nodes.sort()

    def modules(self):
        start_idx = 0
        for edge in self.edges:
            if not any(e.start < edge.end < e.end or (e.start > edge.start and e.end == edge.end) for e in self.edges):
                i = bisect_left(self.nodes, edge.node)
                yield self.Module(self.nodes[start_idx: i + 1])
                start_idx = i


if __name__ == "__main__":
    graph = Graph('ENSMUSG00000039361')
    graph.populate()
    graph.find_connections()
    for module in graph.modules():
        if module.alternate_downstream():
            print('module:', module.nodes[0], module.nodes[-1])
