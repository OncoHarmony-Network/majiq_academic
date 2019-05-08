from bisect import bisect_left, bisect_right
from itertools import combinations
from pathlib import Path

from voila import constants
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_utils import generate_means

from operator import itemgetter
import csv

import argparse
from voila.classifier.tsv_writer import TsvWriter
from voila.config import ClassifyConfig
from voila.exceptions import GeneIdNotFoundInVoilaFile, VoilaException

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



PSI_THRESHOLD = 0.0
DPSI_THRESHOLD = None
HIDE_SUB_COMPLEX = False




class UnsupportedVoilaFile(Exception):
    pass

class Printable_Event:

    def range_str(self):
        return '{}-{}'.format(self.start, self.end)


class Graph:
    def __init__(self, gene_id):
        """
        This contains the edges and nodes used to find modules.

        :param gene_id: gene id
        :param sg_file: splice graph file name
        :param voila_file: voila file name (psi/delta psi)
        """

        self.nodes = []  # all nodes in the graph
        self.edges = []  # all edges in the graph
        self.gene_id = gene_id
        self.config = ClassifyConfig()

        # populate the graph with data from the splice graph
        self._populate()

        # get psi/dpsi data from voila file and associate with edges
        with Matrix(self.config.voila_file) as m:
            analysis_type = m.analysis_type

        with SpliceGraph(self.config.splice_graph_file) as sg:
            gene_meta = sg.gene(self.gene_id)
            if not gene_meta:
                raise VoilaException("Gene ID not found in SpliceGraph File: %s" % self.gene_id)
            self.strand, self.gene_name, self.chromosome = itemgetter('strand', 'name', 'chromosome')(gene_meta)


        if analysis_type == constants.ANALYSIS_PSI:
            self._psi()
        elif analysis_type == constants.ANALYSIS_DELTAPSI:
            self._delta_psi()
        else:
            raise UnsupportedVoilaFile()

        # find connections between nodes
        self._find_connections()

    class Node(Printable_Event):
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

        @property
        def is_half_exon(self):
            return self.exon['end'] == -1 or self.exon['start'] == -1

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
            connected = []
            for edge in edges:
                if edge.node == node:
                    connected.append(edge)
            return connected

    class Edge(Printable_Event):
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

        with SpliceGraph(self.config.splice_graph_file) as sg:
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
        :return: List of modules
        """

        modules = []
        j = 0
        nextEndShift = 0
        start_idx = 0
        for edge in self.edges:


            if not any(e.start < edge.end < e.end or (e.start > edge.start and e.end == edge.end) for e in self.edges):
                j += 1
                i = bisect_left(self.nodes, edge.node)
                if j == 2 or True:
                    #print(edge.lsvs)
                    if((self.nodes[i].end == -1 or self.nodes[i].start == -1) and True):
                        # handling case like exon 19-20 in ENSMUSG00000021820
                        # we aim to make sure that the half exons are in the middle of the module
                        # so that we don't mark the next module as having that half exon
                        modules.append(self.Module(self.nodes[start_idx: i + 1 + 1], self.strand))
                        nextEndShift = 1
                    else:
                        modules.append(self.Module(self.nodes[start_idx + nextEndShift: i + 1], self.strand))
                        nextEndShift = 0

                start_idx = i

        if self.strand == '-':
            modules.reverse()
        for i, mod in enumerate(modules, 1):
            mod.set_idx(i)


        return modules

    def _psi(self):
        """
        When a psi voila file is supplied, this is where the psi data is added to the junctions.
        :return: None
        """

        lsv_store = {}

        with Matrix(self.config.voila_file) as m:
            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):
                lsv = m.psi(lsv_id)
                experiment_names = [m.experiment_names[0][0]]

                for (start, end), means in zip(lsv.junctions, lsv.get('means')):

                    key = str(start) + '-' + str(end)

                    if key not in lsv_store:
                        lsv_store[key] = {'matrix_events': {}}

                    if lsv_id not in lsv_store[key]:
                        lsv_store[key][lsv_id] = {'psi': [], 'delta_psi': []}

                    lsv_store[key][lsv_id]['psi'].append(means)

                    if lsv.intron_retention:
                        with SpliceGraph(self.config.splice_graph_file) as sg:
                            # print(start)
                            # print(end)
                            # print("true 1")

                            if [x for x in sg.intron_retention_reads_exp({'gene_id': self.gene_id,
                                                              'start': int(start)+1,
                                                              'end': int(end)-1}, experiment_names)]:
                                # print("true 2")
                                lsv_store[key]['matrix_events']['IR'] = {'start': int(start)+1,
                                                                         'end': int(end)-1}








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

        with Matrix(self.config.voila_file) as m:

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
            self.source_lsv_ids, self.target_lsv_ids = self.get_lsv_ids()

        def set_idx(self, idx):
            self.idx = idx


        def get_lsv_ids(self):
            """
            This is just used for outputting to tsv later
            :return:
            """
            sources = set()
            targets = set()
            for node in self.nodes:
                for edge in node.edges:
                    for lsv in edge.lsvs:
                        if lsv != 'matrix_events':
                            if ":s:" in lsv:
                                sources.add(lsv)
                            elif ":t:" in lsv:
                                targets.add(lsv)
            return sources, targets

        def strand_case(self, case_plus, case_minus):
            """

            """
            if self.Filters.strand == '+':
                return case_plus
            else:
                return case_minus

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

        def cassette_exon(self):
            """
            Check if exon skipping occurs in this module.

            """

            found = []
            b = self.Filters.target_source_psi
            s = self.Filters.source_psi
            t = self.Filters.target_psi

            for n1, n2, n3 in combinations(self.nodes, 3):
                # print(n1, n2, n3)
                # print('--------------------')
                # print(n1.connects(n2, t))
                # print('--------------------')
                # print(n1.connects(n3, b))
                # print('--------------------')
                # print(n2.connects(n3, s))
                include1 = n1.connects(n2, s)
                include2 = n2.connects(n3, t)
                skip = n1.connects(n3, b)
                if include1 and include2 and skip:
                    # assert len(include1) > 1
                    # assert len(include2) > 1
                    # assert len(skip) > 1
                    found.append({'event': 'cassette_exon', 'C1': n1, 'C2': n3,
                                  'A': n2, 'Include1': include1,
                                  'Include2': include2, 'Skip': skip})
                    #return True

            return found

        def multi_exon_skipping(self):
            """
            Check if multi exon skipping occurs in this module.
            :return: boolean
            """
            found = []

            if len(self.nodes) < 4:
                return []

            b = self.Filters.target_source_psi

            # exclude half exons
            full_exons = list(filter(lambda ex: ex.start != -1 and ex.end != -1, self.nodes))

            if len(full_exons) < 4:
                return []

            # find combinations of nodes with at least two nodes in between
            # we assume nodes are ordered by default
            # we don't care to find different ordering either
            for i, n1 in enumerate(full_exons):
                for j, n2 in enumerate(full_exons):
                    if j - i > 2:
                        skip = n1.connects(n2, b)
                        if skip:
                            include1 = n1.connects(self.nodes[i+1], b)
                            include2 = self.nodes[j].connects(n2, b)
                            includes = []
                            found.append({'event': 'multi_exon_skipping', 'C1': n1, 'C2': n2, 'As': self.nodes[i+1:j],
                                          'Skip': skip, 'Include1': include1, 'Include2': include2,
                                          'Includes': includes})
            return found

        def mutually_exclusive(self):
            """
            Check if mutually exclusive occurs in this module.
            :return: boolean
            """
            found = []
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
                                    found.append({'event': 'mutually_exclusive',
                                                  'C1': n1, 'C2': i1.node,
                                                  'A1': e1.node, 'A2': e2.node,
                                                  'Include1': e1, 'Include2': i2,
                                                  'SkipA1': e2, 'SkipA2': i1})
            return found

        def intron_retention(self):
            """
            Check if intron retention occurs in this module.
            :return: boolean
            """
            found = []

            #b = self.Filters.target_source_psi
            #s = self.Filters.source_psi
            #t = self.Filters.target_psi

            for n1, n2 in combinations(self.nodes, 2):
                for edge in n1.edges:
                    if 'IR' in edge.lsvs.get('matrix_events', {}) and (n1.connects(n2) or n2.connects(n1)):
                        found.append({'event': 'intron_retention', 'C1': n1, 'C2': n2,
                                      'Intron': Graph.Edge(edge.lsvs['matrix_events']['IR'])})

                for edge in n2.edges:
                    if 'IR' in edge.lsvs.get('matrix_events', {}) and (n1.connects(n2) or n2.connects(n1)):
                        found.append({'event': 'intron_retention', 'C1': n1, 'C2': n2,
                                      'Intron': Graph.Edge(edge.lsvs['matrix_events']['IR'])})


            return found




        def alt5ss(self):

            found = []

            for n1, n2 in combinations(self.nodes, 2):
                connections = n1.connects(n2)
                if connections and len(connections) > 1 and \
                        len(set((self.strand_case(x.end, x.start) for x in connections))) == 1:
                    closest_edge = connections[0]
                    pop_i = 0
                    for i, edge in enumerate(connections):
                        if edge.start > closest_edge.start:
                            closest_edge = edge
                            pop_i = i
                    proximal = connections.pop(pop_i)
                    for distal in connections:
                        found.append({'event': 'alt5ss', 'E1': n1, 'E2': n2, 'Proximal': proximal, 'Distal': distal})
            return found

        def alt3ss(self):

            found = []

            for n1, n2 in combinations(self.nodes, 2):
                connections = n1.connects(n2)
                if connections and len(connections) > 1 and \
                        len(set((self.strand_case(x.start, x.end) for x in connections))) == 1:
                    closest_edge = connections[0]
                    pop_i = 0
                    for i, edge in enumerate(connections):
                        if edge.start > closest_edge.start:
                            closest_edge = edge
                            pop_i = i
                    proximal = connections.pop(pop_i)
                    for distal in connections:
                        found.append({'event': 'alt3ss', 'E1': n1, 'E2': n2, 'Proximal': proximal, 'Distal': distal})
            return found

        # def alt3and5ss(self):
        #
        #     found = []
        #
        #     for n1, n2 in combinations(self.nodes, 2):
        #         connections = n1.connects(n2)
        #         if connections and len(connections) > 1 and len(set((x.start for x in connections))) == 1:
        #             for i, edge in enumerate(connections):
        #                 if edge.end == n2.start:
        #                     # this should be the 'Proximal' connection
        #                     proximal = connections.pop(i)
        #                     break
        #             else:
        #                 print("Warning: did not find proximal connection for alt3ss %s %s" % (n1, n2))
        #                 continue
        #             for distal in connections:
        #                 found.append({'event': 'alt3ss', 'E1': n1, 'E2': n2, 'Proximal': proximal, 'Distal': distal})
        #     return found


        def p_alt_last_exon(self):
            found = []
            for node in self.nodes:
                if node.start == -1:
                    found.append({'event': 'p_ale', 'A1': node, 'A2': None, 'C1': None,
                                      'SkipA2': None, 'SkipA1': None})
            return found

        def p_alt_first_exon(self):
            found = []
            for node in self.nodes:
                if node.end == -1:
                    found.append({'event': 'p_afe', 'A1': node, 'A2': None, 'C1': None,
                                      'SkipA2': None, 'SkipA1': None})
            return found

        # def alt_last_exon(self):
        #
        #     found = []
        #     t = self.Filters.target_psi
        #     s = self.Filters.source_psi
        #
        #     for n1, n2, n3 in combinations(self.nodes, 3):
        #         if len(n2.edges) == 1 and len(n3.edges) == 1 and n2.end < n3.start:
        #
        #             connections1 = n1.connects(n2, s)
        #             connections2 = n1.connects(n3, s)
        #             connectionsTest = n2.connects(n3)
        #             if len(connections1) == 1 and len(connections2) == 1 and not connectionsTest:
        #                 found.append({'event': 'ale', 'A1': n2, 'A2': n3, 'C1': n1,
        #                               'SkipA2': connections1, 'SkipA1': connections2})
        #     return found

        def alt_last_exon(self):

            found = []

            nodes = [x for x in reversed(self.nodes)]
            current_submodule = []
            for i, node in enumerate(nodes[:-1]):

                current_submodule.append(node)

                for fwd_node in nodes[:len(current_submodule) - 1]:
                    #print(node.connects(fwd_node))
                    if node.connects(fwd_node):
                        break
                else:

                    if nodes[-1].connects(node) and not nodes[i + 1].connects(node):
                        # end the submodule

                        found.append(
                            {'event': self.strand_case('ale', 'afe'), 'A1': current_submodule[0], 'A2': node, 'C1': nodes[-1],
                             'SkipA2': current_submodule[0].edges, 'SkipA1': node.edges})
                        current_submodule = []

            return found

        # def alt_first_exon(self):
        #
        #     found = []
        #     t = self.Filters.target_psi
        #
        #     for n1, n2, n3 in combinations(self.nodes, 3):
        #         # we care that there is some exon which is connected by skipping one or more exons
        #         # but NOT connected to the skipped exons
        #         # here, when looking for this second condition
        #         if n1.end < n2.start:
        #
        #             connections1 = n1.connects(n3, t)
        #             connections2 = n2.connects(n3, t)
        #             connectionsTest = n1.connects(n2)
        #             if len(connections1) == 1 and len(connections2) == 1 and not connectionsTest:
        #                 found.append({'event': 'afe', 'A1': n1, 'A2': n2, 'C1': n3,
        #                               'SkipA2': connections1, 'SkipA1': connections2})
        #
        #     return found



        def alt_first_exon(self):

            found = []

            current_submodule = []
            for i, node in enumerate(self.nodes[:-1]):
                #print(node)

                current_submodule.append(node)

                # assume these are ordered?
                # we care that there is some exon which is connected by skipping one or more exons
                # but NOT connected to the skipped exons
                # here, when looking for this second condition, we need to do iteration
                # we check if the first exon is connected to the next, etc, until finding a break
                # to find the break we check that the current group is only connected to the last node
                # not other nodes ahead of it.

                for back_node in self.nodes[:len(current_submodule)-1]:
                    if back_node.connects(node):
                        break
                else:
                    # did not find any connections backward

                    if node.connects(self.nodes[-1]) and not node.connects(self.nodes[i+1]):
                        # end the submodule

                        found.append({'event': self.strand_case('afe', 'ale'), 'A1': current_submodule[0], 'A2': node,
                                      'C1': self.nodes[-1], 'SkipA2': current_submodule[0].edges, 'SkipA1': node.edges})
                        current_submodule = []

            #
            # for n1, n2, n3 in combinations(self.nodes, 3):
            #     # we care that there is some exon which is connected by skipping one or more exons
            #     # but NOT connected to the skipped exons
            #     # here, when looking for this second condition
            #     if n1.end < n2.start:
            #
            #         connections1 = n1.connects(n3, t)
            #         connections2 = n2.connects(n3, t)
            #         connectionsTest = n1.connects(n2)
            #         if len(connections1) == 1 and len(connections2) == 1 and not connectionsTest:
            #             found.append({'event': 'afe', 'A1': n1, 'A2': n2, 'C1': n3,
            #                           'SkipA2': connections1, 'SkipA1': connections2})

            return found

        def as_types(self):
            """
            Helper function that returns a list of types found in this module.
            :return: list of AS types, flag is module is complex true or false
            """
            print('---------------------------', self.idx, '--------------------------------')
            as_type_dict = {
                # 'alt_downstream': self.alternate_downstream,
                # 'alt_upstream': self.alternate_upstream,
                'cassette_exon': self.cassette_exon,
                'mutually_exclusive': self.mutually_exclusive,
                'intron_retention': self.intron_retention,
                'alt3ss': self.alt3ss,
                'alt5ss': self.alt5ss,
                'p_alt_last_exon': self.p_alt_last_exon,
                'p_alt_first_exon': self.p_alt_first_exon,
                'alt_last_exon': self.alt_last_exon,
                'alt_first_exon': self.alt_first_exon,
                'multi_exon_skipping': self.multi_exon_skipping
            }
            ret = []
            complex = False
            for k, v in as_type_dict.items():
                res = v()
                ret += res
            if len(ret) > 1:
                complex = True
            if HIDE_SUB_COMPLEX and complex:
                ret = []
            return ret, complex



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



if __name__ == "__main__":

    import pprint

    parser = argparse.ArgumentParser(description='LSV Classification Script')
    parser.add_argument('psi_path', type=check_file,
                        help='Path to either psi file or directory containing psi files')
    parser.add_argument('splicegraph_file', type=check_file,
                        help='File path of input splicegraph file')
    parser.add_argument('-d', '--directory', required=True, help="Output Directory for TSV files")
    parser.add_argument('--hide-sub-complex', action='store_true',
                        help='If a module is detected as complex, do not show any other events detected in that module')
    parser.add_argument('--psi-threshold', type=float, default=0.0,
                        help='Only detect junctions where target and source psi values pass threshold. '
                             '(0.0, the default, accepts everything)')

    args = parser.parse_args()

    PSI_THRESHOLD = args.psi_threshold
    HIDE_SUB_COMPLEX = args.hide_sub_complex
    #dpsi_file = '~/Development/small_test/majiq_deltapsi_all_v_all/Adr_Cer.deltapsi.voila'

    psi_file = args.psi_path
    sg_file = args.splicegraph_file

    sg_file = '/home/paul/PycharmProjects/majiq/test_cases/VOILA_DEV/ORIGINAL_BUILD/splicegraph.sql'
    voila_files = ['/home/paul/PycharmProjects/majiq/test_cases/VOILA_DEV/ORIGINAL_MULTIPSI/test.psi.voila',
                 '/home/paul/PycharmProjects/majiq/test_cases/VOILA_DEV/ORIGINAL_MULTIPSI/test2.psi.voila']

    sg_file = Path(sg_file).expanduser().resolve()
    voila_files = [Path(x).expanduser().resolve() for x in voila_files]

    # Find all gene ids in splice graph
    with SpliceGraph(sg_file) as sg:
       gene_ids = list(g['id'] for g in sg.genes())

    # Find all gene ids in voila file
    # with Matrix(Path(psi_file).expanduser()) as m:
    #     gene_ids = list(m.gene_ids)



    # for gene_id in gene_ids:
    # gene_id = 'ENSMUSG00000001419'
    genes_modules = []
    # genes_events = {'complex', 'cassette_exon', 'mutually_exclusive',
    #                 'intron_retention', 'alt3ss', 'alt5ss', 'altTranscStart', 'altTranscEnd',
    #                 'multi_exon_skipping'}
    out_rows = []
    TsvWriter.delete_tsvs(args.directory)
    for gene_id in gene_ids:
        #if gene_id == "ENSMUSG00000021820":
        #if gene_id == "ENSMUSG00000026843":
        #if gene_id == "ENSMUSG00000024097":
        #if gene_id == "ENSMUSG00000006498":
        #if gene_id == "ENSMUSG00000026843":
        if gene_id == "ENSMUSG00000026843":
        #if gene_id == "ENSMUSG00000049550":
        #     print(gene_id)
        #     graph = Graph(gene_id, sg_file, psi_file)
        #
        #     mod = [x for x in graph.modules()][-3]
        #     print(mod.nodes)
        #     print(mod.as_types())

            graph = Graph(gene_id, sg_file, voila_files)


            writer = TsvWriter(args.directory, graph, gene_id, voila_files, )
            writer.cassette()
            writer.alt3prime()
            writer.alt5prime()
            writer.alt3and5prime()
            writer.mutually_exclusive()
            writer.intron_retention()
            writer.summary()



            # break
            # # genes_modules.append((gene_id, graph.modules()))
            # for i, module in enumerate(graph.modules()):
            #     print(i+1)
            #     for as_type in module.as_types():
            #         # row = {'module_id': i, 'lsv_ids': semicolon(module.lsv_ids), 'gene_id': gene_id,
            #         #        'gene_name': graph.gene_name, 'chr': graph.chromosome, 'strand': graph.strand}
            #
            #         print(as_type)
            #         #pprint.pprint(as_type)
            #     #break
            #         #out_rows.append()
            #     # t = timeit.Timer(module.as_types)
            #     # print(t.timeit(100), module.as_types())
            #
            #     #print(module.as_types())



    #pprint.pprint(out_rows)

    #output_tsv(genes_modules)

    # what to do with exon 19-20 event in http://localhost:5005/gene/ENSMUSG00000021820/


