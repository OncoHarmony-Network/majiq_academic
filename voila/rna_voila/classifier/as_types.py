from bisect import bisect_left, bisect_right
from itertools import combinations, permutations, product
from pathlib import Path
from rna_voila.voila_log import voila_log
from rna_voila import constants
from rna_voila.api import SpliceGraph, Matrix
from rna_voila.api.matrix_utils import generate_means, unpack_bins, generate_high_probability_non_changing, \
    generate_bins_prior_removed, generate_prior_removed_expected_dpsi
from rna_voila.api import view_matrix
from rna_voila.exceptions import LsvIdNotFoundInVoilaFile

from operator import itemgetter
import csv

import argparse
from rna_voila.classifier.tsv_writer import TsvWriter
from rna_voila.config import ClassifyConfig
from rna_voila.exceptions import GeneIdNotFoundInVoilaFile, VoilaException
from rna_voila.vlsv import get_expected_psi, matrix_area, get_expected_dpsi
import numpy as np

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

#np.seterr(all='raise')


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
    def __init__(self, gene_id, experiment_names):
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
        self.experiment_names = experiment_names




        # populate the graph with data from the splice graph
        self._populate()



        with SpliceGraph(self.config.splice_graph_file) as sg:
            gene_meta = sg.gene(self.gene_id)
            if not gene_meta:
                raise VoilaException("Gene ID not found in SpliceGraph File: %s" % self.gene_id)
            self.strand, self.gene_name, self.chromosome = itemgetter('strand', 'name', 'chromosome')(gene_meta)
            if self.config.junc_gene_dist_column:
                self.gene_start, self.gene_end = sg.gene_extent(self.gene_id)

        self.priors = {}

        # find connections between nodes
        self._find_connections()





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

    def _confidence_changing(self, module):
        """
        Post Module-building discarding of modules -- this does not actually remove modules in this function,
        but it calculates if the current group of nodes will pass the confidence filter (otherwise, no module is
        even made, to save on efficiency)

        at least N juncs have:
            any(max(Prob(|E(dPSI)|>=changing-thresh)))>=probability-changing-thresh
        """

        all_edge_lsv_ids = set((y for x in module.get_all_edges(ir=True) for y in x.lsvs.keys()))

        for voila_file in self.config.voila_files:
            with Matrix(voila_file) as m1:
                analysis_type = m1.analysis_type

            if analysis_type == constants.ANALYSIS_HETEROGEN:
                with view_matrix.ViewHeterogen(voila_file) as m:

                    for lsv_id in all_edge_lsv_ids:
                        try:
                            lsv = m.lsv(lsv_id)
                            if lsv.changing(self.config.changing_pvalue_threshold,
                                            self.config.changing_between_group_dpsi).any():
                                return True
                        except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                            continue

            elif analysis_type == constants.ANALYSIS_DELTAPSI:
                with view_matrix.ViewDeltaPsi(voila_file) as m:

                    for lsv_id in all_edge_lsv_ids:
                        try:
                            lsv = m.lsv(lsv_id)
                            if lsv.changing(self.config.changing_between_group_dpsi,
                                            self.config.probability_changing_threshold):
                                return True
                        except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
                            continue

        return False


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

        # if self.config.min_experiments < 1:
        #     # we need to find based on percentage of experiments
        #     num_acceptable = 0.0
        #     percent_acceptable = 0.0
        #     exps = [x for x in reads]
        #
        #     for exp in exps:
        #         if exp['reads'] >= self.config.decomplexify_reads_threshold:
        #             num_acceptable += 1
        #             percent_acceptable = num_acceptable / float(len(exps))
        #
        #         if percent_acceptable >= self.config.min_experiments:
        #             return True
        # else:
        #     # we need to find based on raw number of acceptable experiments
        #     num_acceptable = 0
        #     for exp in reads:
        #         if exp['reads'] >= self.config.decomplexify_reads_threshold:
        #             num_acceptable += 1
        #         if num_acceptable >= self.config.min_experiments:
        #             return True
        #
        # return False

    def _populate(self):
        """
        Add all juctions and exons to graph and sort those lists.
        :return: None
        """

        with SpliceGraph(self.config.splice_graph_file) as sg:
            for exon in sg.exons(self.gene_id):
                self._add_exon(exon)
            for junc in sg.junctions(self.gene_id, omit_simplified=True):
                if self.config.decomplexify_reads_threshold == 0 or self._enough_reads(
                        sg.junction_reads_exp(junc, self.experiment_names)):
                    self._add_junc(junc)

            for ir in sg.intron_retentions(self.gene_id, omit_simplified=True):
                if self.config.decomplexify_reads_threshold == 0 or self._enough_reads(
                        sg.intron_retention_reads_exp(ir, self.experiment_names)):
                    self._add_junc(ir, ir=True)

        # remove exons that don't have any junctions
        # this is done by looking at the start and end of each junction and seeing if any of those ends
        # fall inside of each node

        self._decomplexify()
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

                    voila_log().warning(
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
                            if self._module_is_valid(module):
                                modules.append(module)
                            nextEndShift = 1
                        else:
                            module = self.Module(self.nodes[start_idx + nextEndShift: i + 1], self)
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

    def _add_matrix_values(self, voila_file):
        with Matrix(voila_file) as m:
            analysis_type = m.analysis_type
        if analysis_type == constants.ANALYSIS_PSI:
            self._psi(voila_file)
        elif analysis_type == constants.ANALYSIS_DELTAPSI:
            self._delta_psi(voila_file)
        elif analysis_type == constants.ANALYSIS_HETEROGEN:
            self._heterogen(voila_file)
        else:
            raise UnsupportedVoilaFile()

    def _add_lsvs_to_edges(self, lsv_store):

        for edge in self.edges:


            if edge.ir:
                key = str(edge.start + 1) + '-' + str(edge.end - 1)
            else:
                key = str(edge.start) + '-' + str(edge.end)
            # if we found values for this edge in the voila file

            if key in lsv_store:
                # if the edge already has values found from another voila file
                if edge.lsvs:
                    # for each lsv in the voila file just read
                    for lsv_id in lsv_store[key]:
                        if lsv_id in edge.lsvs:
                            for means in lsv_store[key][lsv_id]['psi']:
                                edge.lsvs[lsv_id]['psi'].add(means)
                            for means in lsv_store[key][lsv_id]['delta_psi']:
                                edge.lsvs[lsv_id]['delta_psi'].add(means)
                            for group_name in lsv_store[key][lsv_id]['group_psi']:
                                if not group_name in edge.lsvs[lsv_id]['group_psi']:
                                    edge.lsvs[lsv_id]['group_psi'][group_name] = set()
                                edge.lsvs[lsv_id]['group_psi'][group_name].update(lsv_store[key][lsv_id]['group_psi'][group_name])
                        else:
                            edge.lsvs[lsv_id] = lsv_store[key][lsv_id]
                else:
                    edge.lsvs = lsv_store[key]

    def _psi(self, voila_file):
        """
        When a psi voila file is supplied, this is where the psi data is added to the junctions.
        :return: None
        """

        lsv_store = {}


        with Matrix(voila_file) as m:
            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):
                lsv = m.psi(lsv_id)



                for (start, end), means in zip(lsv.junctions, lsv.get('means')):

                    key = str(start) + '-' + str(end)



                    if key not in lsv_store:
                        lsv_store[key] = {}

                    if lsv_id not in lsv_store[key]:
                        lsv_store[key][lsv_id] = {'psi': set(), 'delta_psi': set(), 'voila_file': voila_file,
                                                  'group_psi': {}}

                    lsv_store[key][lsv_id]['psi'].add(means)



        self._add_lsvs_to_edges(lsv_store)

    def _delta_psi(self, voila_file):
        """
        When a delta psi voila file is supplied, this is where the psi/delta psi data is added to the junctions.
        :return: None
        """

        lsv_store = {}

        with Matrix(voila_file) as m:

            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):

                lsv = m.delta_psi(lsv_id)
                juncs = lsv.junctions

                for group_means in lsv.get('group_means'):
                    for (start, end), means in zip(juncs, group_means):
                        key = str(start) + '-' + str(end)

                        if key not in lsv_store:
                            lsv_store[key] = {}

                        if lsv_id not in lsv_store[key]:
                            lsv_store[key][lsv_id] = {'psi': set(), 'delta_psi': set(), 'voila_file': voila_file,
                                                      'group_psi': {}}

                        lsv_store[key][lsv_id]['psi'].add(means)


                _bins = lsv.get('bins')

                for (start, end), means, bins in zip(juncs, generate_means(_bins), _bins):
                    key = str(start) + '-' + str(end)
                    lsv_store[key][lsv_id]['delta_psi'].add(means)


        self._add_lsvs_to_edges(lsv_store)


    def _heterogen(self, voila_file):
        """
        When a psi voila file is supplied, this is where the psi data is added to the junctions.
        :return: None
        """

        lsv_store = {}

        with view_matrix.ViewHeterogen(voila_file) as m:

            for lsv_id in m.lsv_ids(gene_ids=[self.gene_id]):
                lsv = m.lsv(lsv_id)

                for junc_i, junc in enumerate(lsv.junctions):

                    key = str(junc[0]) + '-' + str(junc[1])

                    if key not in lsv_store:
                        lsv_store[key] = {}

                    if lsv_id not in lsv_store[key]:
                        lsv_store[key][lsv_id] = {'psi': set(), 'delta_psi': set(), 'voila_file': voila_file,
                                                  'group_psi': {}}

                    junc_median_psi = lsv.median_psi()[junc_i]

                    for i, group_name in enumerate(m.group_names):
                        if not group_name in lsv_store[key][lsv_id]['group_psi']:
                            lsv_store[key][lsv_id]['group_psi'][group_name] = set()

                        lsv_store[key][lsv_id]['group_psi'][group_name].add(junc_median_psi[i])

                    dpsi_val = lsv.median_psi()[junc_i][0] - lsv.median_psi()[junc_i][1]

                    lsv_store[key][lsv_id]['delta_psi'].add(dpsi_val)

        self._add_lsvs_to_edges(lsv_store)

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

            
        def other_event(self, edges_of_interest, lsvs):
            """
            Catch-all for modules lacking events that match classic AS type definitions
            
            
            :return: catch-all event
            """
            nodes_other = []
            edges_other = []
            if len(edges_of_interest) > 0:
                for n1, n2 in permutations(self.nodes, 2):
                    connections = n1.connects(n2)
                    for conn in connections:
                        # other.tsv only ouptuts junction coordinates
                        # that were *not* already accounted for in another event in the module
                        if conn not in edges_of_interest:
                            continue
                        edges_other.append(conn)
                        for n in [n1, n2]:
                            if n not in nodes_other:
                                nodes_other.append(n)


            other_event = {
                'event': 'other_event',
                'nodes': nodes_other,
                'edges': edges_other,
                'lsvs': lsvs
            }
            return [other_event]
        

        def cassette_exon(self):
            """
            Check if exon skipping occurs in this module.

            """

            found = []
            i = 0
            # b = self.Filters.target_source_psi
            # s = self.Filters.source_psi
            # t = self.Filters.target_psi
            # print(self.nodes)
            for n1, n2, n3 in combinations(self.nodes, 3):
                # print(n1, n2, n3)
                # print('--------------------')
                # print(n1.connects(n2))
                # print('--------------------')
                # print(n1.connects(n3))
                # print('--------------------')
                # print(n2.connects(n3))
                # list of edges connecting n1 and n2
                include1s = n1.connects(n2)
                # list of edges connecting n2 and n3
                include2s = n2.connects(n3)
                # list of edges connecting n1 and n3
                skips = n1.connects(n3)
                if include1s and include2s and skips:
                    assert len(include1s) >= 1
                    assert len(include2s) >= 1
                    assert len(skips) >= 1
                    for include1, include2, skip in product(include1s, include2s, skips):
                        # ensure skipping junction shares start with include 1
                        # and skipping junction shares end with include 2
                        if include1.start == skip.start and include2.end == skip.end:
                            # update Module's seen junctions
                            # will only see LSVs associated with skip junction
                            # include junctions might have extra LSV with ref exon as alternative exon
                            # this LSV won't ever be associated with the cassette (appropriately)
                            for edge in [skip]:
                                if len(edge.lsvs) > 0:
                                    self.classified_lsvs.extend(edge.lsvs)
                                self.classified_junctions.extend([include1, include2, edge])

                            event_constitutive = []
                            if self.graph.config.cassettes_constitutive_column:
                                indices = [self.nodes.index(x) for x in (n1, n2, n3)]
                                if list(range(min(indices), max(indices)+1)) == sorted(indices) or \
                                   list(reversed(range(min(indices), max(indices)+1))) == sorted(indices, reverse=True):
                                    event_constitutive = [True]
                                else:
                                    event_constitutive = [False]

                            found.append({'event': 'cassette_exon',
                                          'C1': self.strand_case(n1, n3),
                                          'C2': self.strand_case(n3, n1),
                                          'A': n2,
                                          'Include1': self.strand_case(include1, include2),
                                          'Include2': self.strand_case(include2, include1),
                                          'Skip': skip,
                                          'event_id': 'CE_%s' % i,
                                          'constitutive': event_constitutive})
                            i += 1

            return found

        def tandem_cassette(self):
            """
            Strat: find all multi-exon skipping events, then, check if all of the middle
            exons only have two connections, and are connected to the exon before and after them
            """

            found = []
            event_index = 0

            if len(self.nodes) < 4:
                return []

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
                        # list of edges connecting n1 and n2
                        skips = n1.connects(n2)
                        if len(skips) >= 1:
                            conns = []
                            for k in range(j - i):
                                conns.append(self.nodes[i+k].connects(self.nodes[i+k+1]))

                            # checking all A exons connect together
                            if not all(len(x) > 0 for x in conns):
                                continue

                            # list of edges that connect C1 and A1
                            include1s = self.strand_case(n1.connects(self.nodes[i + 1]), self.nodes[j-1].connects(n2))
                            # list of edges that connect C2 and A_last
                            include2s = self.strand_case(self.nodes[j-1].connects(n2), n1.connects(self.nodes[i + 1]))

                            if len(include1s)>0 and len(include2s)>0:
                                for skip, include1, include2 in product(skips, include1s, include2s):
                                    # checking that first and last connections match skip coordinates
                                    if self.strand_case(include1.start == skip.start and include2.end == skip.end,
                                                        include1.end == skip.end and include2.start == skip.start):

                                        c1 = self.strand_case(n1, n2)
                                        c2 = self.strand_case(n2, n1)

                                        # update Module's seen junctions
                                        for edge in [skip]:
                                            if len(edge.lsvs) > 0:
                                                self.classified_lsvs.extend(edge.lsvs)
                                            self.classified_junctions.extend([include1, include2, edge])

                                        # build include list
                                        includes = []
                                        for exon_i in range(i+1, j-1):
                                            includes.append(self.nodes[exon_i].connects(self.nodes[exon_i+1])[0])

                                        found.append({'event': 'tandem_cassette', 'C1': c1,
                                                      'C2': c2,
                                                      #'As': self.nodes[i + 1:j], redundant?
                                                      'Skip': skip, 'Include1': include1,
                                                      'Include2': include2,
                                                      'Tandem_Exons': self.nodes[i + 1:j],
                                                      'Includes': includes,
                                                      'event_id':'TCE_%s' % event_index})
                                        event_index += 1
            return found

        def multi_exon_spanning(self):
            """
            Check if multi exon skipping occurs in this module.
            :return: boolean
            """
            found = []

            if len(self.nodes) < 4:
                return []

            #b = self.Filters.target_source_psi

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
                        skip = n1.connects(n2)
                        if skip:

                            for possible_i in range(i+1, j):
                                include1s = n1.connects(self.nodes[possible_i], ir=True)
                                if len(include1s) > 0:
                                    break
                            else:
                                include1s = (None,)

                            for possible_i in range(j-1, i, -1):
                                include2s = self.nodes[possible_i].connects(n2, ir=True)
                                if len(include2s) > 0:
                                    break
                            else:
                                include2s = (None,)
                            for sk in skip:
                                for include1 in include1s:
                                    for include2 in include2s:
                                        found.append({'event': 'multi_exon_spanning', 'C1': self.strand_case(n1, n2),
                                                      'C2': self.strand_case(n2, n1), 'As': self.nodes[i+1:j],
                                                      'Skip': sk, 'Include1': self.strand_case(include1, include2),
                                                      'Include2': self.strand_case(include2, include1),
                                                      })
                                        # update Module's seen junctions/lsvs iff we've appended to found
                                        if len(sk.lsvs) > 0:
                                            self.classified_lsvs.extend(sk.lsvs)
                                        self.classified_junctions.append(sk)
            return found

        def mutually_exclusive(self):
            """
            Check if mutually exclusive occurs in this module.
            :return: boolean
            """
            found = []
            #f = self.Filters.target_source_psi

            # for n1, n2, n3, n4 in combinations(self.nodes, 4):
            #     if n1.connects(n2, f) and n1.connects(n3, f) and n2.connects(n4, f) and n3.connects(n4, f):
            #         if not n2.connects(n3):
            #             return True

            for n1 in self.nodes[:-1]:
                for e1, e2, in combinations(n1.edges, 2):
                    # this check removes exitrons from consideration
                    if not (e1.start > e1.node.start and e1.end < e1.node.end) and \
                       not (e2.start > e2.node.start and e2.end < e2.node.end) and \
                       not (e1.ir or e2.ir):
                        if e1.node != e2.node and not (e1.node.connects(e2.node) or e2.node.connects(e1.node)):
                            for i1 in e1.node.edges:
                                for i2 in e2.node.edges:
                                    if i1.node == i2.node and not (i1.ir or i2.ir):
                                        skipA1 = self.strand_case(e2, i1)
                                        include1 = self.strand_case(e1, i2)
                                        skipA2 = self.strand_case(i1, e2)
                                        include2 = self.strand_case(i2, e1)
                                        # update Module's seen junctions
                                        for edge1, edge2 in zip([skipA1, skipA2], [include1, include2]):
                                            # the C1 source or C2 target LSVs only!
                                            shared_lsv = set(edge1.lsvs) & set(edge2.lsvs)
                                            if len(shared_lsv) == 1:
                                                self.classified_lsvs.append(shared_lsv.pop())
                                        self.classified_junctions.extend([skipA1, skipA2,include1, include2])
                                        found.append({'event': 'mutually_exclusive',
                                                      'C1': self.strand_case(n1, i1.node),
                                                      'C2': self.strand_case(i1.node, n1),
                                                      'A1': self.strand_case(e1.node, e2.node),
                                                      'A2': self.strand_case(e2.node, e1.node),
                                                      'Include1': include1,
                                                      'SkipA1': skipA1,
                                                      'Include2': include2,
                                                      'SkipA2': skipA2})

            return found

        def alternative_intron(self):
            """
            Check if intron retention occurs in this module.
            """
            found = []

            for n1, n2 in combinations(self.nodes, 2):
                fwd_connects = n1.connects(n2, only_ir=True)
                for edge in fwd_connects:
                    if edge.ir:
                        if len(n1.edges) > 1 or len(n2.back_edges) > 1:
                            #spliced = n1.connects(n2)
                            spliced = [j for j in n1.edges if not j.ir]
                            for junc in n2.back_edges:
                                if not any(j == junc for j in spliced) and not junc.ir:
                                    spliced.append(junc)
                            for junc in spliced:
                                if len(junc.lsvs)>0:
                                    self.classified_lsvs.extend(junc.lsvs)
                            self.classified_junctions.extend(spliced)
                            self.classified_junctions.append(edge)
                            for non_intron_junc in spliced:
                                found.append({'event': 'alternative_intron', 'C1': n1, 'C2': n2,
                                              'Intron': edge, 'Spliced': non_intron_junc})

            return found



        def alt5ss(self):

            found = []

            for n1, n2 in combinations(self.nodes, 2):
                connections = n1.connects(n2)
                #print(connections, n1)


                    
                    
                if len(connections) > 1:

                    for i in range(len(connections)):

                        e1 = connections[i]  # should be proximal for this case

                        for e2 in connections[:i] + connections[i+1:]:

                            # check that two of the ends meet at the same coordinate
                            if len(set((self.strand_case(x.end, x.start) for x in (e1, e2,)))) == 1:

                                # check that e1 is proximal and e2 is distal
                                if self.strand_case(e1.start, e2.end) > self.strand_case(e2.start, e1.end):
                                    proximal = e1
                                    distal = e2
                                    n_e1 = self.strand_case(n1, n2)
                                    n_e2 = self.strand_case(n2, n1)

                                    # update seen junctions in module
                                    # preferentially add target LSV, just like the TSV writer does
                                    associated_lsv = None
                                    for lid in proximal.lsvs:
                                        if associated_lsv is None:
                                            associated_lsv = lid
                                        if ":t:" in lid:
                                            associated_lsv = lid
                                    if associated_lsv is not None:
                                        self.classified_lsvs.append(associated_lsv)
                                    self.classified_junctions.append(proximal)
                                    self.classified_junctions.append(distal)
                                    found.append({'event': 'alt5ss', 'E1': n_e1, 'E2': n_e2,
                                                  'Proximal': proximal, 'Distal': distal})

            return found


        def alt3ss(self):

            found = []

            for n1, n2 in combinations(self.nodes, 2):
                connections = n1.connects(n2)
                if len(connections) > 1:

                    for i in range(len(connections)):

                        e1 = connections[i]  # should be proximal for this case
                        for e2 in connections[:i] + connections[i+1:]:

                            # check that two of the ends meet at the same coordinate
                            if len(set((self.strand_case(x.start, x.end) for x in (e1, e2,)))) == 1:

                                # check that e1 is proximal and e2 is distal
                                if self.strand_case(e2.end, e1.start) > self.strand_case(e1.end, e2.start):
                                    proximal = e1
                                    distal = e2
                                    n_e1 = self.strand_case(n1, n2)
                                    n_e2 = self.strand_case(n2, n1)

                                    # update seen junctions in module
                                    # preferentially add source LSV, just like the TSV writer does
                                    associated_lsv = None
                                    for lid in proximal.lsvs:
                                        if associated_lsv is None:
                                            associated_lsv = lid
                                        if ":s:" in lid:
                                            associated_lsv = lid
                                    if associated_lsv is not None:
                                        self.classified_lsvs.append(associated_lsv)
                                    self.classified_junctions.append(proximal)
                                    self.classified_junctions.append(distal)
                                    found.append({'event': 'alt3ss', 'E1': n_e1, 'E2': n_e2,
                                                  'Proximal': proximal, 'Distal': distal})
            return found


        def p_alt5ss(self):

            found = []

            for n1, n2, n3 in combinations(self.nodes, 3):

                skips = n1.connects(n3)

                if self.graph.strand == '+':
                    # iterate all nodes in between to check if they all connect with introns

                    include2s = n2.connects(n3)
                    include1s = n1.connects(n2, only_ir=True)

                    if not include1s:
                        start_idx = self.nodes.index(n1)
                        end_idx = self.nodes.index(n2)
                        for ni in range(start_idx, end_idx):
                            if not self.nodes[ni].connects(self.nodes[ni+1], only_ir=True):
                                break
                        else:
                            include1s = self.nodes[start_idx].connects(self.nodes[start_idx+1], only_ir=True)

                else:
                    include2s = n1.connects(n2)
                    include1s = n2.connects(n3, only_ir=True)
                    if not include1s:
                        start_idx = self.nodes.index(n2)
                        end_idx = self.nodes.index(n3)
                        for ni in range(start_idx, end_idx):
                            if not self.nodes[ni].connects(self.nodes[ni+1], only_ir=True):
                                break
                        else:
                            include1s = self.nodes[start_idx].connects(self.nodes[start_idx+1], only_ir=True)

                if include1s and include2s and skips:

                    for include1, include2, skip in product(include1s, include2s, skips):
                        if len(set((self.strand_case(x.end, x.start) for x in (skip, include2,)))) == 1:

                            # skip when there are any junctions connecting the two intron-connected exons
                            # which share coordinate with the skip junction
                            if self.graph.strand == '+':
                                if any(_j.start == skip.start for _j in n1.connects(n2)):
                                    continue
                            else:
                                if any(_j.end == skip.end for _j in n2.connects(n3)):
                                    continue

                            # update seen junctions in Module
                            if len(skip.lsvs) > 0:
                                self.classified_lsvs.extend(skip.lsvs)
                            self.classified_junctions.extend([skip, include1, include2])
                            found.append({'event': 'p_alt5ss', 'C1': n1, 'C2': n3,
                                      'A': n2, 'Include1': include1,
                                      'Include2': include2, 'Skip': skip})
                        # found.append({'event': 'alt5ss', 'E1': n1, 'E2': n2, 'Proximal': proximal, 'Distal': distal})

            return found

        def p_alt3ss(self):

            found = []

            #b = self.Filters.target_source_psi

            for n1, n2, n3 in combinations(self.nodes, 3):

                skips = n1.connects(n3)

                if self.graph.strand == '+':
                    include1s = n1.connects(n2)
                    include2s = n2.connects(n3, only_ir=True)

                    if not include2s:
                        start_idx = self.nodes.index(n2)
                        end_idx = self.nodes.index(n3)
                        for ni in range(start_idx, end_idx):
                            if not self.nodes[ni].connects(self.nodes[ni + 1], only_ir=True):
                                break
                        else:
                            include2s = self.nodes[start_idx].connects(self.nodes[start_idx + 1], only_ir=True)

                else:
                    include1s = n2.connects(n3)
                    include2s = n1.connects(n2, only_ir=True)

                    if not include2s:
                        start_idx = self.nodes.index(n1)
                        end_idx = self.nodes.index(n2)
                        for ni in range(start_idx, end_idx):
                            if not self.nodes[ni].connects(self.nodes[ni+1], only_ir=True):
                                break
                        else:
                            include2s = self.nodes[start_idx].connects(self.nodes[start_idx+1], only_ir=True)


                if include1s and include2s and skips:

                    for include1, include2, skip in product(include1s, include2s, skips):
                        if len(set((self.strand_case(x.start, x.end) for x in (include1, skip,)))) == 1:

                            # skip when there are any junctions connecting the two intron-connected exons
                            # which share coordinate with the skip junction
                            if self.graph.strand == '+':
                                if any(_j.end == skip.end for _j in n2.connects(n3)):
                                    continue
                            else:
                                if any(_j.start == skip.start for _j in n1.connects(n2)):
                                    continue

                            # update seen junctions in Module
                            if len(skip.lsvs) > 0:
                                self.classified_lsvs.extend(skip.lsvs)
                            self.classified_junctions.extend([skip, include1, include2])
                            found.append({'event': 'p_alt3ss', 'C1': n1, 'C2': n3,
                                          'A': n2, 'Include1': include1,
                                          'Include2': include2, 'Skip': skip})
                        # found.append({'event': 'alt5ss', 'E1': n1, 'E2': n2, 'Proximal': proximal, 'Distal': distal})

            return found

        def alt3and5ss(self):

            found = []

            for n1, n2 in combinations(self.nodes, 2):
                connections = n1.connects(n2)
                if connections and len(connections) > 1:
                    for e1, e2 in combinations(connections, 2):
                        if e1.start != e2.start and e1.end != e2.end:
                            # update seen junctions in Module
                            for edge in [e1, e2]:
                                if len(edge.lsvs) > 0:
                                    self.classified_lsvs.extend(edge.lsvs)
                                self.classified_junctions.append(edge)
                            found.append(
                                {'event': 'alt3and5ss', 'E1': self.strand_case(n1, n2),
                                 'E2': self.strand_case(n2, n1),
                                 'J1': self.strand_case(e1, e2), 'J2': self.strand_case(e2, e1)})
                    #
                    #
                    # closest_edge = connections[0]
                    # pop_i = 0
                    # for i, edge in enumerate(connections):
                    #     if edge.start > closest_edge.start:
                    #         closest_edge = edge
                    #         pop_i = i
                    # proximal = connections.pop(pop_i)
                    # for distal in connections:
                    #     found.append({'event': 'alt3ss', 'E1': n1, 'E2': n2, 'Proximal': proximal, 'Distal': distal})
            return found


        def p_alt_last_exon(self):
            found = []
            for node in self.nodes:

                if not node.is_half_exon:
                    continue

                all_half_exons = True
                for other_node in self.nodes:
                    connections = node.connects(other_node, ir=True) + other_node.connects(node, ir=True)

                    if connections:
                        if not other_node.is_half_exon:
                            all_half_exons = False
                        if self.graph.strand == '+':
                            if other_node.start > node.start:
                                break
                        else:
                            if other_node.start < node.start:
                                break
                else:
                    if all_half_exons:
                        continue

                    # first get the junction(s) connecting from the half_exon backwards
                    for junc_from_he in self.strand_case(node.back_edges, node.edges):

                        ref_exon = self.strand_case(self.graph.start_node(junc_from_he), self.graph.end_node(junc_from_he))

                        # then, look for junctions coming out of the found 'reference' exon pointed in the same
                        # direction as the half-exon
                        for junc_to_another_exon in self.strand_case(ref_exon.edges, ref_exon.back_edges):
                            other_exon = self.strand_case(self.graph.end_node(junc_to_another_exon),
                                                          self.graph.start_node(junc_to_another_exon))
                            if other_exon == node:
                                continue
                            shared_lsv = set(junc_from_he.lsvs) & set(junc_to_another_exon.lsvs)
                            # only classify if LSV exists
                            # this set is empty for example if a p_ale/p_afe 
                            # has the two required junctions but there were not 
                            # enough reads for the quantifier to see an LSV.
                            if len(shared_lsv) > 0:
                                this_lsv = shared_lsv.pop()
                                if not ":s:" in this_lsv:
                                    continue
                                if len(shared_lsv) == 1:
                                    self.classified_lsvs.append(this_lsv)


                            if len(junc_from_he) > len(junc_to_another_exon):
                                proximal = other_exon
                                distal = node
                                skipA1 = junc_from_he
                                skipA2 = junc_to_another_exon
                            else:
                                proximal = node
                                distal = other_exon
                                skipA1 = junc_to_another_exon
                                skipA2 = junc_from_he
                            found.append({'event': 'p_ale', 'Proximal': proximal,
                                          'Distal': distal, 'Reference': ref_exon,
                                          'SkipA2': skipA2,
                                          'SkipA1': skipA1})

                            self.classified_junctions.append(skipA1)
                            self.classified_junctions.append(skipA2)


            return found

        def p_alt_first_exon(self):
            found = []
            for node in self.nodes:

                if not node.is_half_exon:
                    continue

                all_half_exons = True
                for other_node in self.nodes:

                    connections = node.connects(other_node, ir=True) + other_node.connects(node, ir=True)
                    if connections:
                        if not other_node.is_half_exon:
                            all_half_exons = False
                        if self.graph.strand == '+':
                            if other_node.start < node.start:
                                break
                        else:
                            if other_node.start > node.start:
                                break
                else:
                    if all_half_exons:
                        continue

                    # first get the junction(s) connecting from the half_exon backwards
                    for junc_from_he in self.strand_case(node.edges, node.back_edges):

                        ref_exon = self.strand_case(self.graph.end_node(junc_from_he), self.graph.start_node(junc_from_he))

                        # then, look for junctions coming out of the found 'reference' exon pointed in the same
                        # direction as the half-exon
                        for junc_to_another_exon in self.strand_case(ref_exon.back_edges, ref_exon.edges):
                            other_exon = self.strand_case(self.graph.start_node(junc_to_another_exon),
                                                          self.graph.end_node(junc_to_another_exon))
                            if other_exon == node:
                                continue

                            shared_lsv = set(junc_from_he.lsvs) & set(junc_to_another_exon.lsvs)
                            # only classify if LSV exists
                            # this set is empty for example if a p_ale/p_afe 
                            # has the two required junctions but there were not 
                            # enough reads for the quantifier to see an LSV.
                            if len(shared_lsv) > 0:
                                this_lsv = shared_lsv.pop()
                                if not ":t:" in this_lsv:
                                    continue
                                if len(shared_lsv) == 1:
                                    self.classified_lsvs.append(this_lsv)

                            if len(junc_from_he) > len(junc_to_another_exon):
                                proximal = other_exon
                                distal = node
                                skipA1 = junc_to_another_exon
                                skipA2 = junc_from_he
                            else:
                                proximal = node
                                distal = other_exon
                                skipA1 = junc_from_he
                                skipA2 = junc_to_another_exon
                            
                            found.append({'event': 'p_afe', 'Proximal': proximal,
                                          'Distal': distal, 'Reference': ref_exon,
                                          'SkipA2': skipA2,
                                          'SkipA1': skipA1})
                            #print(found)
                            self.classified_junctions.append(skipA1)
                            self.classified_junctions.append(skipA2)


            return found


        def alt_first_exon(self):

            found = []
            for node in self.nodes:

                if node.is_half_exon:
                    continue

                for other_node in self.nodes:
                    connections = node.connects(other_node, ir=True) + other_node.connects(node, ir=True)
                    if connections:
                        if self.graph.strand == '+':
                            if node == self.nodes[0]:

                                break
                            if other_node.start < node.start:
                                break
                        else:
                            if node == self.nodes[-1]:

                                break
                            if other_node.start > node.start:
                                break
                else:

                    # first get the junction(s) connecting from the afe backwards
                    for junc_from_afe in self.strand_case(node.edges, node.back_edges):



                        ref_exon = self.strand_case(self.graph.end_node(junc_from_afe), self.graph.start_node(junc_from_afe))

                        # then, look for junctions coming out of the found 'reference' exon pointed in the same
                        # direction as the afe
                        for junc_to_another_exon in self.strand_case(ref_exon.back_edges, ref_exon.edges):
                            other_exon = self.strand_case(self.graph.start_node(junc_to_another_exon),
                                                          self.graph.end_node(junc_to_another_exon))
                            if other_exon == node:
                                continue

                            if len(junc_from_afe) > len(junc_to_another_exon):
                                proximal = other_exon
                                distal = node
                                skipA1 = junc_to_another_exon
                                skipA2 = junc_from_afe
                            else:
                                proximal = node
                                distal = other_exon
                                skipA1 = junc_from_afe
                                skipA2 = junc_to_another_exon

                            shared_lsv = set(junc_from_afe.lsvs) & set(junc_to_another_exon.lsvs)
                            if len(shared_lsv) == 1:
                                self.classified_lsvs.append(shared_lsv.pop())
                            found.append({'event': 'afe', 'Proximal': proximal,
                                  'Distal': distal, 'Reference': ref_exon,
                                  'SkipA2': skipA2,
                                  'SkipA1': skipA1})
                            self.classified_junctions.append(skipA2)
                            self.classified_junctions.append(skipA1)

            return found

        def alt_last_exon(self):

            found = []
            for node in self.nodes:

                if node.is_half_exon:
                    continue

                for other_node in self.nodes:
                    connections = node.connects(other_node, ir=True) + other_node.connects(node, ir=True)

                    if connections:
                        if self.graph.strand == '+':
                            if node == self.nodes[-1]:
                                break
                            if other_node.start > node.start:
                                break
                        else:
                            if node == self.nodes[0]:
                                break
                            if other_node.start < node.start:
                                break
                else:


                    # first get the junction(s) connecting from the half_exon backwards
                    for junc_from_ale in self.strand_case(node.back_edges, node.edges):

                        ref_exon = self.strand_case(self.graph.start_node(junc_from_ale), self.graph.end_node(junc_from_ale))

                        # then, look for junctions coming out of the found 'reference' exon pointed in the same
                        # direction as the half-exon
                        for junc_to_another_exon in self.strand_case(ref_exon.edges, ref_exon.back_edges):
                            other_exon = self.strand_case(self.graph.end_node(junc_to_another_exon),
                                                          self.graph.start_node(junc_to_another_exon))
                            if other_exon == node:
                                continue


                            if len(junc_from_ale) > len(junc_to_another_exon):
                                proximal = other_exon
                                distal = node
                                skipA1 = junc_from_ale
                                skipA2 = junc_to_another_exon
                            else:
                                proximal = node
                                distal = other_exon
                                skipA1 = junc_to_another_exon
                                skipA2 = junc_from_ale

                            shared_lsv = set(junc_from_ale.lsvs) & set(junc_to_another_exon.lsvs)
                            if len(shared_lsv) == 1:
                                self.classified_lsvs.append(shared_lsv.pop())
                            found.append({'event': 'ale', 'Proximal': proximal,
                                  'Distal': distal, 'Reference': ref_exon,
                                  'SkipA2': skipA2,
                                  'SkipA1': skipA1})
                            self.classified_junctions.append(skipA1)
                            self.classified_junctions.append(skipA2)

            return found

        def orphan_junction(self):
            found = []
            for n1, n2 in combinations(self.nodes, 2):

                if n1.is_half_exon and n2.is_half_exon:
                    conn = n1.connects(n2)
                    for _conn in conn:
                        # update seen lsvs in Module
                        if len(_conn.lsvs) > 0:
                            self.classified_lsvs.extend(_conn.lsvs)
                        self.classified_junctions.append(_conn)
                        found.append({'event': 'orphan_junction', 'A1': n1, 'A2': n2, 'Junc': _conn})
            return found


        def exitron(self):
            found = []
            # look for junctions which are in the same exon
            for node in self.nodes:
                for edge in node.edges:
                    if edge.start > node.start and edge.end < node.end:
                        found.append({'event': 'exitron', 'Exon': node, 'Junc': edge})
            return found

        def constitutive(self):
            found = []
            # why populate nodes to check and not just check all combos of 2?
            nodes_to_check = []
            if len(self.nodes[0].edges_no_exitrons) == 1:
                nodes_to_check.append(self.nodes[0])
            if len(self.nodes[-1].back_edges_no_exitrons) == 1:
                nodes_to_check.append(self.nodes[-1])
            for node in self.nodes[1:-1]:
                if len(node.edges_no_exitrons) == 1 and len(node.back_edges_no_exitrons) == 1:
                    nodes_to_check.append(node)

            for n1, n2 in combinations(nodes_to_check, 2):
                junc = n1.connects(n2, ir=True) # TBD; set ir=False?
                if junc and len(junc) == 1:

                    junc = junc[0]
                    config = ClassifyConfig()

                    if not junc.ir: # TBD; if set ir=False, no need to check here?
                        with SpliceGraph(config.splice_graph_file) as sg:
                            try:
                                junc_reads = next(sg.junction_reads_exp({'start': junc.start, 'end': junc.end,
                                                                    'gene_id':self.graph.gene_id},
                                                                     self.graph.experiment_names))['reads']
                            except StopIteration:
                                continue

                            if junc_reads >= config.keep_constitutive:
                                # update seen junctions in Module
                                if len(junc.lsvs) > 0:
                                    self.classified_lsvs.extend(junc.lsvs)
                                self.classified_junctions.append(junc)
                                found.append({'event': 'constitutive', 'Junc': junc, 'Spliced': junc,
                                              'C1': n2, 'C2': n1})

            return found

        def constitutive_intron(self):
            """
            Check if intron retention occurs in this module.
            """
            found = []

            for n1, n2 in combinations(self.nodes, 2):
                fwd_connects = n1.connects(n2, only_ir=True)
                for edge in fwd_connects:
                    if edge.ir:
                        if len(n1.edges) == 1 and len(n2.back_edges) == 1:
                            config = ClassifyConfig()
                            with SpliceGraph(config.splice_graph_file) as sg:
                                try:
                                    junc_reads = \
                                    next(sg.intron_retention_reads_exp({'start': edge.start + 1, 'end': edge.end - 1,
                                                                        'gene_id': self.graph.gene_id},
                                                                       self.graph.experiment_names))['reads']
                                except StopIteration:
                                    continue

                                if junc_reads >= config.keep_constitutive:
                                    # update seen junctions in Module
                                    if len(edge.lsvs) > 0:
                                        self.classified_lsvs.extend(edge.lsvs)
                                    self.classified_junctions.append(edge)
                                    found.append({'event': 'constitutive_intron', 'C1': n1, 'C2': n2,
                                              'Intron': edge, 'Spliced': edge})

            return found


        # Appropriate to be a static Class method?
        def module_is_constitutive(self):
            """
            Detects whether the module is a pair of constitutively (by junction or intron) joined exons ...
            :return: Boolean
            """
            if len(self.nodes) == 2:
                n1 = self.nodes[0]
                n2 = self.nodes[1]
                connection = n1.connects(n2)
                # if connected by single junction
                if connection and len(connection)==1:
                    # if constitutive...
                    if len(n1.edges) == 1 and len(n2.back_edges) == 1:
                        return True
            return False


        def mpe(self):
            """
            Identify splice graph splits and associated constitutive regions (upstream if source, downstream if target).

            Event dicts structured as follows:
            {
            'event': 'mpe_source' or 'mpe_target',
            'reference_exon': <reference exon node>
            'at_module_edge': Boolean : is the reference exon at the appropriate edge the module?
                (appropriate meaning left if source and right if target)
            'constitutive_regions': <[list of nodes (or introns) that are upstream constitutive]>
            'edge_of_transcript': True/False is ref exon first/last exon of a transcript, potentially?
            }
            """
            found = []
            # skip certain types of modules (i.e. constitutive)
            if self.module_is_constitutive():
                return found
            for ii in range(len(self.nodes)):
                thisnode = self.nodes[ii]
                # node has with 2+ (forward) edges, so single source
                if len(thisnode.edges) > 1:
                    constitutive_regions = []
                    backedges = thisnode.back_edges
                    #if len(backedges) == 1:
                        #print(" ## backedge of thisnode %s is %s " % (thisnode, backedges))
                    while len(backedges) <= 1:
                        if len(backedges) == 0:
                            break
                        # only one backedge...
                        backedge = backedges[0]
                        #print("backedge: %s" % backedge)
                        if backedge.ir:
                            constitutive_regions.append(backedge)
                        # "end_node" is actually upstream, b/c neg strand..
                        upstream_node = self.graph.start_node(backedge)
                        #print("upstream: %s, upstream_node.edges %s" % (upstream_node, upstream_node.edges))
                        # if upstream node only forward connects to thisnode
                        if len(upstream_node.edges) == 1:
                            constitutive_regions.append(upstream_node)
                        else:
                            break
                        backedges = upstream_node.back_edges
                    if self.graph.strand == "+":
                        event_type = "mpe_source"
                    else:
                        event_type = "mpe_target"
                    if ii == 0:
                        at_edge = True
                    else:
                        at_edge = False
                    if ii == (len(self.nodes) - 1):
                        at_opposite_edge = True
                    else:
                        at_opposite_edge = False
                    # is this potentially first or last exon in a transcript?
                    if len(thisnode.back_edges) == 0:
                        edge_of_transcript = True
                    else:
                        edge_of_transcript = False
                    if not at_opposite_edge:
                        found.append(
                            {
                            'event': event_type,
                            'reference_exon': thisnode,
                            'at_module_edge': at_edge,
                            'constitutive_regions': constitutive_regions,
                            'edge_of_transcript':edge_of_transcript
                            })
                # node has with 2+ (forward) edges, so single source
                if len(thisnode.back_edges) > 1:
                    constitutive_regions = []
                    backedges = thisnode.edges
                    # if len(backedges) == 1:
                    #     print(" ## edge of thisnode %s is %s " % (thisnode, backedges))
                    while len(backedges) <= 1:
                        if len(backedges) == 0:
                            break
                        # only one backedge...
                        backedge = backedges[0]
                        # print("backedge: %s" % backedge)
                        if backedge.ir:
                            constitutive_regions.append(backedge)
                        upstream_node = self.graph.end_node(backedge) # Fix?
                        # print("upstream: %s, upstream_node.back_edges %s" % (upstream_node, upstream_node.back_edges))
                        # if upstream node only forward connects to thisnode
                        if len(upstream_node.back_edges) == 1:
                            constitutive_regions.append(upstream_node)
                        else:
                            break
                        backedges = upstream_node.edges
                    if self.graph.strand == "+":
                        event_type = "mpe_target"
                    else:
                        event_type = "mpe_source"
                    if ii == (len(self.nodes) - 1):
                        at_edge = True
                    else:
                        at_edge = False
                    if ii == 0:
                        at_opposite_edge = True
                    else:
                        at_opposite_edge = False
                    # is this potentially first or last exon in a transcript?
                    if len(thisnode.edges) == 0:
                        edge_of_transcript = True
                    else:
                        edge_of_transcript = False
                    if not at_opposite_edge:
                        found.append(
                            {
                            'event': event_type,
                            'reference_exon': thisnode,
                            'at_module_edge': at_edge,
                            'constitutive_regions': constitutive_regions,
                            'edge_of_transcript':edge_of_transcript
                            })
            return found

        def mpe_regions(self):
            """
            Gets the mpe_source events and the mpe_target events
            :return: [mpe_single_sources] [mpe_single_targets]
            """
            return self.mpe()

        def as_types(self):
            """
            Helper function that returns a list of types found in this module.
            :return: list of AS types, flag is module is complex true or false
            """
            ret = []

            if ClassifyConfig().putative_multi_gene_regions:
                for region in self.p_multi_gene_regions:
                    region['event'] = 'p_multi_gene_region'
                    ret.append(region)
                return ret, False, len(ret)

            # print('---------------------------', self.idx, '--------------------------------')
            # print(self.nodes)
            # print(self.get_all_edges())
            # print([e.de_novo for e in self.get_all_edges()])
            as_type_dict = {
                # 'alt_downstream': self.alternate_downstream,
                # 'alt_upstream': self.alternate_upstream,
                'cassette_exon': self.cassette_exon,
                'mutually_exclusive': self.mutually_exclusive,
                'alternative_intron': self.alternative_intron,
                'alt3ss': self.alt3ss,
                'alt5ss': self.alt5ss,
                'alt3and5ss': self.alt3and5ss,
                'p_alt3ss': self.p_alt3ss,
                'p_alt5ss': self.p_alt5ss,
                'p_alt_last_exon': self.p_alt_last_exon,
                'p_alt_first_exon': self.p_alt_first_exon,
                'alt_last_exon': self.alt_last_exon,
                'alt_first_exon': self.alt_first_exon,
                'multi_exon_spanning': self.multi_exon_spanning,
                'tandem_cassette': self.tandem_cassette,
                'orphan_junction': self.orphan_junction,
            }

            if ClassifyConfig().keep_constitutive:
                as_type_dict['constitutive'] = self.constitutive
                as_type_dict['constitutive_intron'] = self.constitutive_intron
            event_counts = {
                'cassette_exon': 0,
                'mutually_exclusive': 0,
                'alternative_intron': 0,
                'alt3ss': 0,
                'alt5ss': 0,
                'alt3and5ss': 0,
                'p_alt3ss': 0,
                'p_alt5ss': 0,
                'p_ale': 0,
                'p_afe': 0,
                'afe': 0,
                'ale': 0,
                'multi_exon_spanning': 0,
                'tandem_cassette': 0,
                'orphan_junction': 0,
                'constitutive': 0,
                'constitutive_intron': 0,
            }

            complex = False

            for k, v in as_type_dict.items():
                # call function to check if each event exists
                res = v()
                ret += res

            for e in ret:
                event_counts[e['event']] += 1

            total_events = sum(event_counts.values())
            check_events = total_events - event_counts['constitutive'] - event_counts['constitutive_intron']

            if check_events == 2 and event_counts['multi_exon_spanning'] == 1 and \
                    event_counts['tandem_cassette'] == 1:
                complex = False
            elif check_events == 2 and event_counts['alternative_intron'] == 1 and\
                    (event_counts['p_alt3ss'] == 1 or event_counts['p_alt5ss'] == 1):
                complex = False
            elif check_events > 1:
                complex = True

            non_classified_lsvs = set(self.all_lsvs) - set(self.classified_lsvs)
            non_classified_junctions = set(self.get_all_edges()) - set(self.classified_junctions)

            # If there are LSVs that were not classified
            if len(non_classified_lsvs) > 0:
                # if there are junctions in module that were not classified
                if len(non_classified_junctions) > 0:
                    # then label this complex since the junction(s) isn't in a standard event definition
                    complex = True
                these_t_lsvs = self.target_lsv_ids
                these_s_lsvs = self.source_lsv_ids
                if len(these_t_lsvs) > 0 or len(these_s_lsvs) > 0:
                    # print(these_t_lsvs, these_s_lsvs)
                    # print(self.get_all_edges())
                    # create an 'other_event' for this module
                    ret += self.other_event(edges_of_interest=non_classified_junctions, lsvs=non_classified_lsvs)
                    total_events += 1
                    if complex:
                        ret[-1]['novel'] = True

            if HIDE_SUB_COMPLEX and complex:
                ret = []

            self.is_complex = complex

            return ret, complex, total_events



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
    #                 'alternative_intron', 'alt3ss', 'alt5ss', 'altTranscStart', 'altTranscEnd',
    #                 'multi_exon_spanning'}
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
            writer.alternative_intron()
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


