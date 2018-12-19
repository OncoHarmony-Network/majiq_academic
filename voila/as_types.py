from bisect import bisect_left, bisect_right
from itertools import combinations, chain
from pathlib import Path

from voila.api import SpliceGraph, Matrix
from voila.api.matrix_utils import generate_means

sg_file = '~/Development/small_test/majiq_build/splicegraph.sql'
psi_file = '~/Development/small_test/majiq_psi_all/Adr.psi.voila'
dpsi_file = '~/Development/small_test/majiq_deltapsi_all_v_all/Adr_Cer.deltapsi.voila'

files = {
    'sg': sg_file,
    'psi': psi_file,
    'dpsi': dpsi_file
}

for k in files:
    files[k] = Path(files[k])
    files[k] = files[k].expanduser().resolve()

PSI_THRESHOLD = None
DPSI_THRESHOLD = None


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
            for n1, n2, n3 in combinations(self.nodes, 3):
                if n1.connects(n2) and n1.connects(n3) and n2.connects(n3):
                    return True

        def mutually_exclusive(self):
            for n1, n2, n3, n4 in combinations(self.nodes, 4):
                if n1.connects(n2) and n1.connects(n3) and n2.connects(n4) and n3.connects(n4):
                    if not n2.connects(n3):
                        return True

        def as_type(self):
            as_type_dict = {
                'alt_downstream': self.alternate_downstream,
                'alt_upstream': self.alternate_upstream,
                'exon_skipping': self.exon_skipping,
                'mutually_exclusive': self.mutually_exclusive
            }
            return [k for k, v in as_type_dict.items() if v()]

    class Node:
        def __init__(self, exon):
            self.edges = []
            self.exon = exon

        def __eq__(self, other):
            return self.start == other.start and self.end == other.end

        def __lt__(self, other):
            return self.view_start < other.view_start and self.view_end < other.view_start

        def __repr__(self):
            return '<{} {},{}>'.format(self.__class__.__name__, self.start, self.end)

        def __hash__(self):
            return hash(str(self.start) + '-' + str(self.end))

        @property
        def start(self):
            return self.exon['start']

        @property
        def end(self):
            return self.exon['end']

        @property
        def view_start(self):
            start = self.exon['start']
            if start == -1:
                start = self.end - 10
            return start

        @property
        def view_end(self):
            end = self.exon['end']
            if end == -1:
                end = self.start + 10
            return end

        def connects(self, node):
            return any(e.node == node for e in self.edges)

    class Edge:
        def __init__(self, junc):
            self.junc = junc
            self.node = None
            self.source_psi = []
            self.target_psi = []
            self.delta_psi = None

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

        @property
        def view_start(self):
            return self.start

        @property
        def view_end(self):
            return self.end

    def start_node(self, edge):
        i = bisect_left(self.nodes, edge)
        return self.nodes[i]

    def end_node(self, edge):
        i = bisect_right(self.nodes, edge)
        return self.nodes[i - 1]

    def add_junc(self, junc):
        edge = self.Edge(junc)

        # we can use edge to search for left and right nodes in node list.
        start_node = self.start_node(edge)
        end_node = self.end_node(edge)

        if start_node != end_node:
            self.edges.append(edge)
            edge.node = end_node

    def add_exon(self, exon):
        self.nodes.append(self.Node(exon))

    def find_connections(self):
        for edge in self.edges:
            node = self.start_node(edge)
            node.edges.append(edge)

    def populate(self):
        with SpliceGraph(files['sg']) as sg:
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

    def filter_nodes(self):
        keep_nodes = set()
        for edge in self.edges:
            keep_nodes.add(self.start_node(edge))
            keep_nodes.add(self.end_node(edge))

        self.nodes = list(keep_nodes)
        self.nodes.sort()

    def filter_edges(self):
        if PSI_THRESHOLD is not None:
            self.edges = [e for e in self.edges if e.source_psi or e.target_psi]
            self.edges = [e for e in self.edges if any(p >= PSI_THRESHOLD for p in chain(e.source_psi, e.target_psi))]

        if DPSI_THRESHOLD is not None:
            self.edges = [e for e in self.edges if e.delta_psi is not None]
            self.edges = [e for e in self.edges if abs(e.delta_psi) >= DPSI_THRESHOLD]

    def psi(self):
        target = {}
        source = {}

        with Matrix(files['psi']) as m:
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

    def delta_psi(self):
        target = {}
        source = {}
        dpsi = {}

        with Matrix(files['dpsi']) as m:
            group_names = m.group_names
            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):

                lsv = m.delta_psi(lsv_id)
                juncs = lsv.junctions

                for group, group_means in zip(group_names, lsv.get('group_means')):

                    for (start, end), means in zip(juncs, group_means):
                        key = str(start) + '-' + str(end)
                        if lsv_id.split(':')[-2] == 's':

                            if key in source:
                                source[key].append(means)
                            else:
                                source[key] = [means]

                        else:

                            if key in target:
                                target[key].append(means)
                            else:
                                target[key] = [means]

                for (start, end), means in zip(juncs, generate_means(lsv.get('bins'))):
                    key = str(start) + '-' + str(end)
                    dpsi[key] = means

        for edge in self.edges:
            key = str(edge.start) + '-' + str(edge.end)
            if key in source:
                edge.source_psi += source[key]
            if key in target:
                edge.target_psi += target[key]
            if key in dpsi:
                edge.delta_psi = dpsi[key]


if __name__ == "__main__":

    with SpliceGraph(files['sg']) as sg:
        gene_ids = [g['id'] for g in sg.genes()]

    for gene_id in gene_ids:

        graph = Graph(gene_id)

        # populate the graph with data from the splice graph
        graph.populate()

        # get psi/dpsi data from voila file and associate with an edge
        graph.delta_psi()

        # remove edges that don't pass thresholds
        graph.filter_edges()

        # find connections between nodes
        graph.find_connections()

        # remove nodes that don't have edges
        graph.filter_nodes()

        for module in graph.modules():
            if module.exon_skipping():
                print(gene_id, module.nodes)
