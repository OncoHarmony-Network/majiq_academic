from itertools import combinations, permutations, product
from pathlib import Path
from operator import itemgetter
from rna_voila.api import _SpliceGraphSQL as SpliceGraph
from bisect import bisect_left, bisect_right

from collections import namedtuple

exon = namedtuple('exon', 'start end')

class UnsupportedVoilaFile(Exception):
    pass

class Printable_Event:

    def _ranges_to_string(self, start, end):
        return f'{"na" if start == -1 else start}-{"na" if end == -1 else end}'

    def untrimmed_range_str(self):
        return self._ranges_to_string(getattr(self, 'untrimmed_start', self.start),
                                      getattr(self, 'untrimmed_end', self.end))

    def range_str(self):
        if ClassifyConfig().untrimmed_exons:
            return self.untrimmed_range_str()
        else:
            return self._ranges_to_string(self.start, self.end)

class Graph:
    def __init__(self, gene_id, experiment_names, splice_graph_file):
        """
        This contains the edges and nodes used to find modules.

        :param gene_id: gene id
        :param sg_file: splice graph file name
        :param voila_file: voila file name (psi/delta psi)
        """

        self.nodes = []  # all nodes in the graph
        self.edges = []  # all edges in the graph
        self.gene_id = gene_id

        self.experiment_names = experiment_names
        self.splice_graph_file = splice_graph_file




        # populate the graph with data from the splice graph
        self._populate()



        with SpliceGraph(self.splice_graph_file) as sg:
            gene_meta = sg.gene(self.gene_id)
            if not gene_meta:
                raise Exception("Gene ID not found in SpliceGraph File: %s" % self.gene_id)
            self.strand, self.gene_name, self.chromosome = itemgetter('strand', 'name', 'chromosome')(gene_meta)

        self.priors = {}

        # find connections between nodes
        self._find_connections()

        class C:
            keep_no_lsvs_modules = True
            keep_constitutive = False
            show_all = True
        self.config = C

    class Node(Printable_Event):
        def __init__(self, exon):
            """
            Graph wrapper for exons.

            :param exon: exon dictionary from splice graph file.
            """

            self.edges = []  # all edge starting in this exon
            self.back_edges = []
            self.exon = exon  # exon dictionary
            self.graph = None


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

            return '<{} {} ({}),{} ({})>'.format(self.__class__.__name__, self.start,
                                                 getattr(self, 'untrimmed_start', self.start),
                                                 self.end, getattr(self, 'untrimmed_end', self.end))




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
        def absolute_start(self):
            return self.start

        @property
        def absolute_end(self):
            return self.end

        @property
        def is_half_exon(self):
            return self.exon['end'] == -1 or self.exon['start'] == -1

        @property
        def short_name(self):
            # 'h' for half-exon
            if self.is_half_exon:
                return "h"
            # 'e' for Exon
            return "e"

        def is_de_novo(self):
            return False if self.exon['annotated'] else True

        def get_exitrons(self):
            """
            From a node, return exitron str coordinates in a list. Returns empty list if no exitrons found.
            :return: [<exitron coord>, <exitron coord>, <etc>]f
            """
            exitrons = []
            for edge in self.edges:
                if self.start < edge.start < self.end and self.start < edge.end < self.end:
                    exitrons.append(edge.range_str())
            return exitrons

        @property
        def edges_no_exitrons(self):
            """
            get self.edges, excluding exitrons
            """
            return [x for x in self.edges if not x.is_exitron(self)]

        @property
        def back_edges_no_exitrons(self):
            """
            get self.edges, excluding exitrons
            """
            return [x for x in self.back_edges if not x.is_exitron(self)]

        def get_constant_region(self):
            exitrons = self.get_exitrons()
            if len(exitrons) > 0:
                return ""
            # TBD do something about exitrons...
            # find first start from non-exitron edges that falls within exon bounds
            # if no edges, use node end (i.e. last exon)
            if len(self.edges) == 0:
                first_start = self.end
            else:
                first_start = float("Inf")
                for edge in self.edges:
                    if edge in exitrons:
                        continue
                    if edge.start < first_start and self.start <= edge.start and self.end >= edge.start:
                        first_start = edge.start

            # find last end from non-exitron edges that falls within exon bounds
            # if no back edges, use node start (i.e. first exon)
            if len(self.back_edges) == 0:
                last_end = self.start
            else:
                last_end = float("-Inf")
                for edge in self.back_edges:
                    if edge in exitrons:
                        continue
                    if edge.end > last_end and self.start <= edge.end and self.end >= edge.end:
                        last_end = edge.end
            if last_end >= first_start:
                return "No Constant Region"
            return ("%s-%s" % (last_end, first_start))


        def connects(self, node, filter=None, ir=False, only_ir=False):
            """
            Search through junctions for this exon to see if this exon has a junction that connects to supplied exon.
            :param node: the exon that this exon might connect to.
            :param filter: function to filter junctions
            :param ir: include IR edges and standard junctions if true
            :param only_ir: include ONLY IR edges and NOT standard junctions if true
            :return: list of edges connecting self node and other node, or empty list
            """
            edges = self.edges
            if filter:
                edges = filter(edges)

            if only_ir is True:
                edges = [edge for edge in edges if edge.ir]
            elif ir is False:
                edges = [edge for edge in edges if not edge.ir]

            # print(node)
            # print([x.node for x in edges])
            # print([x.node for x in edges])
            connected = []
            for edge in edges:
                if edge.node == node:
                    connected.append(edge)
            return connected


    class Edge(Printable_Event):
        def __init__(self, junc, ir=False):
            """
            Graph wrapper for junctions.
            :param junc: junction dictionary from the splice graph file.
            """

            self.ir = ir
            self.junc = junc  # junction dictionary
            self.de_novo = True if junc.get('annotated', 1) == 0 else False
            self.node = None  # node this junction connects to
            self.lsvs = {}
            self.is_constitutive = self.junc.get('is_constitutive', False)

        def __lt__(self, other):
            """
            Junction are uniquely identified by gene_id, start, end.  Since graphs work on one gene at a time, we order
            junction by start and end.
            :param other: other junction
            :return: boolean
            """
            #return self.start < other.start and self.end < other.start
            return self.view_start < other.view_start and self.view_end < other.view_start

        def __hash__(self):
            return hash(str(self))

        def __eq__(self, other):
            """
            Equality is determined by start and end of junction.
            :param other:
            :return:
            """

            return self.start == other.start and self.end == other.end and self.ir == other.ir

        def __repr__(self):
            """
            String representation of junction with start and end.
            :return: string
            """

            return '<{} {},{}>'.format(self.__class__.__name__, self.start, self.end)

        def __len__(self):
            return abs(self.start-self.end)

        def range_str(self):
            if ClassifyConfig().untrimmed_exons:
                return self.untrimmed_range_str()
            else:
                return '{}-{}'.format(self.absolute_start, self.absolute_end)

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

        @property
        def absolute_start(self):
            if not self.ir:
                return self.start
            return self.start + 1

        @property
        def absolute_end(self):
            if not self.ir:
                return self.end
            return self.end - 1

        @property
        def short_name(self):
            # 'i' for Intron
            if self.ir:
                return "i"
            else:
                return "j"

        def is_de_novo(self):
            return True if self.de_novo else False

        def is_exitron(self, node):
            if node.is_half_exon or self.ir:
                return False
            return self.start >= node.start and self.end <= node.end


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
        assert i > 0  # this should never be negative / loop to the last / first node
        return self.nodes[i - 1]

    def _add_junc(self, junc, ir=False):
        """
        Add junction to graph as edge object. If junction starts and ends in the same exon, it is not added to graph.
        This function follows decomplexify rules accordingly. Voila files will be read for this junction. If the
        thresholds are not passed, the munction is not added.
        :param junc: junction dictionary from splice graph.
        :param ir: junction is intron retention
        :return: None
        """

        if ir:
            junc['start'] -= 1
            junc['end'] += 1

        edge = self.Edge(junc, ir)

        # Since majiq doesn't quantify junctions that start/stop in same exon, filter them.
        #if start_node != end_node:

        self.edges.append(edge)

    def _add_exon(self, exon):
        """
        Added exon to graph as node object.
        :param exon: exon dictionary from splice graph
        :return: None
        """
        node = self.Node(exon)
        node.idx = "%d_%d" % (exon['start'], exon['end'])
        self.nodes.append(node)

    def _find_connections(self):
        """
        When this has completed, each exon should have a list of junctions that start/end there.
        :return: None
        """

        for edge in self.edges:
            node = self.start_node(edge)
            node.edges.append(edge)
            node = self.end_node(edge)
            node.back_edges.append(edge)


    def _decomplexify(self):
        """
        Pre Module-building discarding of junctions

        :return:
        """
        for voila_file in self.config.voila_files:
            self._add_matrix_values(voila_file)

        #assert False

        for i in range(len(self.edges) - 1, -1, -1):

            if self.edges[i].lsvs:

                group_means_psi = []
                for v in self.edges[i].lsvs.values():
                    for group_vals in v['group_psi'].values():
                        group_means_psi.append(sum(group_vals) / len(group_vals))

                try:
                    max_single_psi = max(map(max, (v['psi'] for v in self.edges[i].lsvs.values())), default=None)
                except ValueError:
                    max_single_psi = None

                max_group_psi = max(group_means_psi) if group_means_psi else None

                if max_single_psi and max_group_psi:
                    psi = max(max_single_psi, max_group_psi)
                elif max_single_psi:
                    psi = max_single_psi
                elif max_group_psi:
                    psi = max_group_psi
                else:
                    psi = None
                delta_psi = max((abs(y) for v in self.edges[i].lsvs.values() for y in v['delta_psi']), default=None)

                # We need both psi and deltapsi to pass threshold to keep
                assert psi is not None or delta_psi is not None

                if psi is not None and psi < self.config.decomplexify_psi_threshold:
                    del self.edges[i]
                elif delta_psi is not None and delta_psi < self.config.decomplexify_deltapsi_threshold:
                    del self.edges[i]

            else:
                if not self.config.keep_no_lsvs_junctions:
                    # if there are no lsvs, and it is not flagged constitutive by Majiq, delete it
                    if not self.edges[i].is_constitutive:
                        del self.edges[i]


    def _remove_empty_exons(self):
        """
        Remove exons / nodes with no junctions
        Also remove nodes that only have either:
            -a forward junction at the very beginning
            or
            -a backward junction at the very end
        These are not supposed to exist but come up in the case of collided exons sometimes.
        """
        new_nodes = []
        for node in self.nodes:
            for i, edge in enumerate(self.edges):
                end_in = self.in_exon(node, edge.end)
                start_in = self.in_exon(node, edge.start)
                # if (exitron junction) or (start in or end in but not start == beginning of exon or end == end of exon)
                if (start_in and end_in) or ((self.in_exon(node, edge.end) or self.in_exon(node, edge.start)) and
                                             not (edge.start == node.start or edge.end == node.end)):
                    new_nodes.append(node)
                    break
            else:
                #print("Removed", node)
                pass

        self.nodes[:] = new_nodes



    def in_exon(self,
                exon,
                coordinate):
        """
        Check if the coordinate falls inside the exon
        :param exon:
        :param coordinate:
        :return:
        """
        if exon.start == -1:
            return coordinate == exon.end
        elif exon.end == -1:
            return coordinate == exon.start
        return exon.start <= coordinate <= exon.end

    def _enough_reads(self, reads):
        """
        We check all the experiments in this group, and find if enough of them are above the read threshold
        (or the % of experiments threshold)
        :param junc: splicegraph object representing the junc
        :return: True of False
        """
        for exp in reads:
            # temporarialy changed to just one experiment overcoming reads threshold is acceptable
            if exp['reads'] >= self.config.decomplexify_reads_threshold:
                return True
        return False

    def _populate(self):
        """
        Add all juctions and exons to graph and sort those lists.
        :return: None
        """

        with SpliceGraph(self.splice_graph_file) as sg:
            for exon in sg.exons(self.gene_id):
                self._add_exon(exon)
            for junc in sg.junctions(self.gene_id, omit_simplified=True):
                # if self.config.decomplexify_reads_threshold == 0 or self._enough_reads(
                #         sg.junction_reads_exp(junc, self.experiment_names)):
                self._add_junc(junc)

            for ir in sg.intron_retentions(self.gene_id, omit_simplified=True):
                # if self.config.decomplexify_reads_threshold == 0 or self._enough_reads(
                #         sg.intron_retention_reads_exp(ir, self.experiment_names)):
                self._add_junc(ir, ir=True)

        # remove exons that don't have any junctions
        # this is done by looking at the start and end of each junction and seeing if any of those ends
        # fall inside of each node

        #self._decomplexify()
        self._remove_empty_exons()
        self._trim_exons()

        self.edges.sort()
        self.nodes.sort()

        # setting end_node wont work properly until sorting is finished
        for edge in self.edges:
            edge.node = self.end_node(edge)



    def _clear_connections(self):
        for node in self.nodes:
            node.edges = []
            node.back_edges = []

    def _trim_exons(self):
        """
        modify self.nodes to remove parts of exons which exist outside of any junction connection
        """
        for i, node in enumerate(self.nodes):
            # find conditions where we should not trim! ---

            # not half exon
            if node.is_half_exon:
                continue

            # first find exitrons, we will need them later
            exitrons = []
            for edge in self.edges:
                # this is different then using in_exon() because it is exclusive instead of inclusive
                # this is how we differentiate exitrons from junctions in overlapping exons
                if node.start < edge.start < node.end and node.start < edge.end < node.end:
                    exitrons.append(edge)

            coords_in_exon = []

            trim_end = False
            for edge in self.edges:
                # look through all edges
                # if we can't find any going ahead (other end greater value than exon) (excluding IR), don't trim end
                if self.in_exon(node, edge.start) and edge.end >= node.end:
                    # intron retention ahead immediately disqualified trimming ahead
                    if edge.ir:
                        trim_end = False
                        break
                    # check that the edge allowing trimming fwd is completely ahead of exitrons, otherwise
                    # if does not count
                    for exitron in exitrons:
                        if edge.start <= exitron.end:
                            break
                    else:
                        coords_in_exon.append(edge.start)
                        trim_end = True

            trim_start = False
            for edge in self.edges:
                # similar for backwards
                if self.in_exon(node, edge.end) and edge.start <= node.start:
                    if edge.ir:
                        trim_start = False
                        break
                    for exitron in exitrons:
                        if edge.end >= exitron.start:
                            break
                    else:
                        coords_in_exon.append(edge.end)
                        trim_start = True

            # need to check for the special case that there are is only one coordinate on the exon where
            # all junctions are connected. In this case we should not trim
            if coords_in_exon and all(x == coords_in_exon[0] for x in coords_in_exon):
                trim_start = False
                trim_end = False

            # end find conditions part ---

            node.untrimmed_start = node.start
            node.untrimmed_end = node.end
            global_min = float('inf')
            global_max = float('-inf')
            if trim_end or trim_start:
                edges_starting = []
                edges_ending = []
                for _e in self.edges:
                    if self.in_exon(node, _e.start) and not _e.start == node.start:
                        global_max = max(_e.start, global_max)
                        global_min = min(_e.start, global_min)
                        edges_starting.append(_e)
                for _e in self.edges:
                    if self.in_exon(node, _e.end) and not _e.end == node.end:
                        global_max = max(_e.end, global_max)
                        global_min = min(_e.end, global_min)
                        edges_ending.append(_e)

            if trim_start:
                node.exon['start'] = global_min

            if trim_end:
                node.exon['end'] = global_max

        # after all the regular trimming is done, there still may be some collision cases due to AFE/ALE that are
        # collided. We look for any remaining collisions, and trim the exon without junctions by one unit to
        # resolve the collision
        for i, node in enumerate(self.nodes[:-1]):
            if node.end == self.nodes[i + 1].start:
                #if not node.edges: # I don't think this is possible at this point? No edges yet?
                if not node.end in [edge.start for edge in self.edges]:
                    node.untrimmed_end = node.end
                    node.exon['end'] -= 1
                #elif not self.nodes[i + 1].back_edges: # I don't think this is possible at this point? No backedges yet?
                elif not node.start in [edge.end for edge in self.edges]:
                    self.nodes[i + 1].untrimmed_start = self.nodes[i + 1].start
                    self.nodes[i + 1].exon['start'] += 1
                else:
                    print(
                        "Found two exons in gene %s which are collided and both have junctions! Can not trim!" % self.gene_id)


    def _module_is_valid(self, module):
        """
        Make sure that module passes checks pertaining to current settings, before being added to module list
        :return:
        """

        # removing modules with no lsv ids
        if not self.config.keep_no_lsvs_modules:
            if not module.source_lsv_ids and not module.target_lsv_ids:
                return False

        # removing modules with only one junction
        if not self.config.keep_constitutive:
            if not module.get_num_edges(ir=True) > 1:
                return False

        # make sure module is changing, unless 'show all' is enabled
        if not self.config.show_all:
            if not self._confidence_changing(module):
                return False

        return True


    def modules(self):
        """
        Search through edges to find where they don't cross.  At this point is where the previous module ends and the
        next module starts. There will be an over lapping exon.
        :return: List of modules
        """

        modules = []
        nextEndShift = 0
        start_idx = 0
        #edges = [x for x in self.edges if not x.ir]  # we exclude IR from module creation for now
        edges = self.edges
        for edge in edges:

            # there is a chance that there are simply no junctions at all, so this algorithm will misplace
            # slightly. It always overlaps by one exon as stated, so we will end up getting two exons
            # smack against each other with no junctions between in a module, which gives superflous
            # definitions.
            # to work around this, we need to check the first exon in a module actually connects to
            # ANY other node alead, and if not remove this node from the module.
            # I just need to think of some efficient way to do this...

            if not any(e.start < edge.end < e.end or (e.start > edge.start and e.end == edge.end) for e in edges):
                i = bisect_left(self.nodes, edge.node)

                if (i - start_idx) > 0:

                    #print(edge.lsvs)
                    if(i < len(self.nodes) and (self.nodes[i].end == -1 or self.nodes[i].start == -1) and True):
                        # handling case like exon 19-20 in ENSMUSG00000021820
                        # we aim to make sure that the half exons are in the middle of the module
                        # so that we don't mark the next module as having that half exon
                        module = self.Module(self.nodes[start_idx: i + 1 + 1], self)
                        module._global_node_start_idx = start_idx
                        module._global_node_end_idx = i + 1 + 1
                        if self._module_is_valid(module):
                            modules.append(module)
                        nextEndShift = 1
                    else:
                        module = self.Module(self.nodes[start_idx + nextEndShift: i + 1], self)
                        module._global_node_start_idx = start_idx + nextEndShift
                        module._global_node_end_idx = i + 1
                        if self._module_is_valid(module):
                            modules.append(module)
                        nextEndShift = 0

                    start_idx = i

        if self.strand == '-':
            modules.reverse()

        # removing beginning node of module if it does not have any forward junctions
        # this can happen when there are complete breaks in the gene
        # this section also serves to calculate the broken gene region "events" for later
        last_break_idx = 0
        p_multi_gene_regions = []
        num_regions_found = 0
        last_region = None
        for i, mod in enumerate(modules, 1):
            if mod.get_num_edges(ir=True) > 0 and not mod.nodes[0].edges:

                if self.strand == '+':
                    # if there are prior entries, we need to update the last one instead of ending
                    # on the exon at the end of the gene, to end on the exon at the end of
                    # the last found region

                    if p_multi_gene_regions:
                        p_multi_gene_regions[-1]['ExonEnd'] = modules[last_break_idx-2].nodes[-1]
                    p_multi_gene_regions.append({'ExonStart': modules[last_break_idx-1].nodes[0] if p_multi_gene_regions else modules[0].nodes[0],
                                                 'ExonEnd': mod.nodes[0],
                                                 'idx': num_regions_found + 1})
                    last_region = {'ExonStart': mod.nodes[1],
                                   'ExonEnd': modules[-1].nodes[-1]}

                else:
                    if p_multi_gene_regions:
                        p_multi_gene_regions[-1]['ExonEnd'] = modules[last_break_idx-1].nodes[0]
                    p_multi_gene_regions.append({'ExonStart': modules[last_break_idx].nodes[-1] if p_multi_gene_regions else modules[0].nodes[-1],
                                                 'ExonEnd': mod.nodes[1],
                                                 'idx': num_regions_found + 1})
                    last_region = {'ExonStart': mod.nodes[0],
                                   'ExonEnd': modules[-1].nodes[0]}

                last_break_idx = i
                num_regions_found += 1

                del mod.nodes[0]  # this line actually removes the problem exon for the module
            mod.set_idx(i)
        if num_regions_found > 0:
            last_region['idx'] = num_regions_found + 1
            p_multi_gene_regions.append(last_region)

        if modules:
            modules[0].p_multi_gene_regions = p_multi_gene_regions

        return modules

    class Module:
        def __init__(self, nodes, graph):
            """
            Module is subset of a gene.  The divide between modules is where junctions don't cross.

            :param nodes: list of nodes that belong to module
            """

            self.nodes = nodes  # subset of nodes for this module
            self.graph = graph
            self.source_lsv_ids, self.target_lsv_ids = self.get_lsv_ids()
            # junctions that have been classified in this module
            self.classified_lsvs = list()
            self.classified_junctions = list()
            # LSV IDs to classify
            self.all_lsvs = self.source_lsv_ids | self.target_lsv_ids
            self.p_multi_gene_regions = []  # meta event

        def set_idx(self, idx):
            self.idx = idx

        def get_all_edges(self, ir=False):
            edges = []
            for node, n2 in permutations(self.nodes, 2):
                connections = node.connects(n2, ir=ir)
                if connections:
                    edges += connections

            return edges

        def get_num_edges(self, ir=False):
            num_edges = 0
            for node, n2 in permutations(self.nodes, 2):
                num_edges += len(node.connects(n2, ir=ir))
            return num_edges

        def get_lsv_ids(self, nodes_in=None):
            """
            This is just used for outputting to tsv later
            :return:
            """
            sources = set()
            targets = set()
            if nodes_in is None:
                nodes_in = self.nodes

            for node, n2 in permutations(nodes_in, 2):
                connections = node.connects(n2)
                if connections:
                    for edge in connections:
                        for lsv in edge.lsvs:
                            if ":s:" in lsv:
                                sources.add(lsv)
                            elif ":t:" in lsv:
                                targets.add(lsv)
            return sources, targets

        def strand_case(self, case_plus, case_minus):
            """

            """
            if self.graph.strand == '+':
                return case_plus
            else:
                return case_minus

