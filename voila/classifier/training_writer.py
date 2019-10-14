from voila.voila_log import voila_log
from voila.config import ClassifyConfig
import numpy as np
import h5py
from voila.classifier.tsv_writer import BaseTsvWriter
import os, csv
from itertools import combinations
import copy, glob

class TrainingWriter(BaseTsvWriter):

    def __init__(self, graph, gene_id):
        """

        :param output_path: The folder where all output TSV files will be written under
        :param graph: the Graph object of the gene
        """
        super().__init__(graph, gene_id)

        self.config = ClassifyConfig()
        self.avg_multival = True
        self.log = voila_log()

        self.graph = graph
        self.gene_id = gene_id

        if self.graph:
            self.modules = self.graph.modules()
            self.resolved_junctions = self._all_junctions()
            #self._split_exons()


    @staticmethod
    def tsv_names():
        config = ClassifyConfig()
        names = []
        if 'exons' in config.enabled_outputs:
            names.append('exons.tsv')
        if 'junctions' in config.enabled_outputs:
            names.append('junctions.tsv')
        if 'paths' in config.enabled_outputs:
            names.append('paths.tsv')
        return names

    @staticmethod
    def delete_hdf5s():
        config = ClassifyConfig()
        paths = glob.glob(config.directory + '/*.hdf5.*')
        for path in paths:
            os.remove(path)

    def start_all_headers(self):

        if 'exons' in self.config.enabled_outputs:
            headers = self.common_headers + ['Exon ID', 'Exon ID Start Coordinate', 'Exon ID End Coordinate']
            self.start_headers(headers, 'exons.tsv')
        if 'junctions' in self.config.enabled_outputs:
            headers = self.common_headers + ['Junc ID', 'Junc ID Start Coordinate', 'Junc ID End Coordinate',
                                             'Exon 1 ID', 'Exon 1 ID Start Coordinate', 'Exon 1 ID End Coordinate',
                                             'Exon 2 ID', 'Exon 2 ID Start Coordinate', 'Exon 2 ID End Coordinate',
                                             'Is Intron']
            self.start_headers(headers, 'junctions.tsv')
        if 'paths' in self.config.enabled_outputs:
            headers = self.common_headers + ['Path ID', 'Path', 'Edge Type', 'Exon Length mod 3']
            self.start_headers(headers, 'paths.tsv')

    def _get_node_id(self, start, end, idx=None):
        # if hasattr(node, 'num_times_split'):
        #     # should be the base node of a split, need to use extended ID
        #     node_id = "{gene_id}_{chr}_{strand}_{max_start}_{max_end}_{idx}".format(
        #         gene_id=self.gene_id,
        #         chr=self.graph.chromosome,
        #         strand=self.graph.strand,
        #         max_start=node.start,
        #         max_end=node.end,
        #         idx=1
        #     )
        # elif hasattr(node, 'split_index'):
        #     # a shorter node made from a split, need to use extended ID
        #     node_id = "{gene_id}_{chr}_{strand}_{max_start}_{max_end}_{idx}".format(
        #         gene_id=self.gene_id,
        #         chr=self.graph.chromosome,
        #         strand=self.graph.strand,
        #         max_start=node.maximal_start,
        #         max_end=node.maximal_end,
        #         idx=node.split_index
        #     )
        if idx is not None:
            node_id = "{gene_id}_{chr}_{strand}_{start}_{end}_{idx}".format(
                gene_id=self.gene_id,
                chr=self.graph.chromosome,
                strand=self.graph.strand,
                start=start,
                end=end,
                idx=idx
            )
        else:
            node_id = "{gene_id}_{chr}_{strand}_{start}_{end}".format(
                gene_id=self.gene_id,
                chr=self.graph.chromosome,
                strand=self.graph.strand,
                start=start,
                end=end
            )
        return node_id

    def paths_tsv(self):
        """
        Paths through a module

        format:
        Path ID:  "module_id"_"chr"_"strand"_"index"





        """
        with open(os.path.join(self.config.directory, 'paths.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            rows = []

            for module_idx, path_idx, paths_found, frameshift in self.paths_through_module():

                common_data = ["%s_%d" % (self.gene_id, module_idx), self.gene_id, self.graph.gene_name,
                               self.graph.chromosome, self.graph.strand]

                rows.append(common_data + ["%s_%d" % (self.gene_id, path_idx),
                                           '_'.join(junc.range_str() for junc in paths_found),
                                           '_'.join(junc.short_name for junc in paths_found), frameshift])



            # for module_idx, e1_start, e1_end, e2_start, e2_end, ir in self.resolved_junctions:
            #
            #     common_data = ["%s_%d" % (self.gene_id, module_idx), self.gene_id, self.graph.gene_name,
            #                    self.graph.chromosome, self.graph.strand]
            #
            #     if not (e1_start, e1_end,) in written_exons:
            #         rows.append(common_data + [self._get_node_id(e1_start, e1_end), e1_start, e1_end])
            #         written_exons.add((e1_start, e1_end,))
            #
            #     if not (e2_start, e2_end,) in written_exons:
            #         rows.append(common_data + [self._get_node_id(e2_start, e2_end), e2_start, e2_end])
            #         written_exons.add((e2_start, e2_end,))

                    # if module.graph.strand == '+':
                    #     rows = rows + tmp
                    # else:
                    #     rows = tmp + rows


            writer.writerows(rows)

    def exons_tsv(self):
        """
        exons.list format:
        Exon ID:  "gene_id"_"chr"_"strand"_"start"_"end".

        If there are two exon entries (an exon has two 5' or 3' splice sites) ,
         then the ID will have "gene_id"_"chr"_"strand"_"maximal-start"_"maximal-end"_1"
          and "gene_id"_"chr"_"strand"_"maximal-start"_"maximal-end"_2"

        The re-written job follows along from the all_junctions, as below.
        For each junction from the originally ran junctions, grab each left and right exon and
        remove duplicates. This should give the required information about the exons.

        """
        with open(os.path.join(self.config.directory, 'exons.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            rows = []
            written_exons = set()
            for module_idx, e1_start, e1_end, e2_start, e2_end, ir in self.resolved_junctions:

                common_data = ["%s_%d" % (self.gene_id, module_idx), self.gene_id, self.graph.gene_name,
                               self.graph.chromosome, self.graph.strand]

                if not (e1_start, e1_end,) in written_exons:
                    rows.append(common_data + [self._get_node_id(e1_start, e1_end), e1_start, e1_end])
                    written_exons.add((e1_start, e1_end,))

                if not (e2_start, e2_end,) in written_exons:
                    rows.append(common_data + [self._get_node_id(e2_start, e2_end), e2_start, e2_end])
                    written_exons.add((e2_start, e2_end,))

                    # if module.graph.strand == '+':
                    #     rows = rows + tmp
                    # else:
                    #     rows = tmp + rows


            writer.writerows(rows)

    def paths_through_module(self):
        """
        For each module, recursively find paths out of the last exon entered in that path and add then to the running
        list of junctions in the path until we reach the last exon in the module. Then output an ordered list of
        paths that were found in each iteration.

        This will end up skipping alt-first and "p-alt-first" because there are no junctions entering
        So we need to manually run the algorithm starting at these nodes as well as the first node.
        :return:
        """


        paths_found = []
        for module in self.modules:
            path_idx = 0
            #print('=====================')

            def add_junctions_out(node, prev_juncs, exon_length, is_module_length):
                nonlocal path_idx
                #print('iter', node, prev_juncs, node.edges, node == module.nodes[-1])
                if not node.edges or node == module.nodes[-1]:
                    if (not node.edges and node != module.nodes[-1]) or node.end == -1:
                        # takes care of checking for ALE cases and p_ALE cases (they will always be at the end)
                        is_module_length = False
                    # got to end of possible path
                    #print(node)
                    exon_length = exon_length + (node.end - prev_juncs[-1].end) + 1
                    if is_module_length:
                        frameshift = str(exon_length % 3)
                    else:
                        frameshift = 'N/A'
                    path_idx += 1
                    paths_found.append((module.idx, path_idx, prev_juncs, frameshift))
                else:
                    for junc in node.edges:
                        _prev_juncs = prev_juncs[:]
                        _prev_juncs.append(junc)
                        next_node = self.graph.end_node(junc)
                        add_junctions_out(next_node, _prev_juncs,
                                      exon_length + (junc.start - node.start) - (junc.end - next_node.start) + 1,
                                          is_module_length)


            # run over junctions from the start of the module
            add_junctions_out(module.nodes[0], [], 0, True)
            # find other possible starting nodes like afe, p_afe, and run over those
            for node in module.nodes[1:-1]:
                if not node.back_edges:
                    add_junctions_out(node, [], 0, False)

        return paths_found






    def junctions_tsv(self):
        with open(os.path.join(self.config.directory, 'junctions.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            rows = []
            for module_idx, e1_start, e1_end, e2_start, e2_end, ir in self.resolved_junctions:

                common_data = ["%s_%d" % (self.gene_id, module_idx), self.gene_id, self.graph.gene_name,
                               self.graph.chromosome, self.graph.strand]

                junc_start, junc_end = e1_end, e2_start

                junc_id = "{gene_id}_{chr}_{strand}_{start}_{end}".format(
                    gene_id=self.gene_id,
                    chr=self.graph.chromosome,
                    strand=self.graph.strand,
                    start=junc_start,
                    end=junc_end
                )

                # if common_data + [junc_id, edge.start, edge.end,
                #                      n1.idx, n1.start, n1.end,
                #                      n2.idx, n2.start, n2.end, edge.ir] in rows:
                #     print("DUP: %s " % str(common_data + [junc_id, edge.start, edge.end,
                #                      n1.idx, n1.start, n1.end,
                #                      n2.idx, n2.start, n2.end, edge.ir]))
                #    # assert False
                #print(edge, n1, n2, edge.end == n2.start, edge.start == n1.end)

                # assert edge.end == n2.start
                # assert edge.start == n1.end

                rows.append(common_data + [junc_id, junc_start, junc_end,
                                     "%d_%d" % (e1_start, e1_end), e1_start, e1_end,
                                     "%d_%d" % (e2_start, e2_end), e2_start, e2_end, ir])

            writer.writerows(rows)




    def adjacency_matrix(self):


        #as_types = {x.idx: x.as_types() for x in modules}



        # print([x for x in modules[0].nodes])
        # print([x for x in modules[0].get_all_edges()])
        files_to_create = ['identity']
        for _type in self.types2headers:
            if _type == 'psi':
                for filename in self.types2headers[_type]:
                    files_to_create.append(filename)
            elif _type == 'dpsi':
                for filename in self.types2headers[_type]:
                    files_to_create.append(filename)


        for filename in files_to_create:

            with h5py.File(os.path.join(self.config.directory, '%s.hdf5.%s' % (filename, self.pid)), "a") as hf:

                for module in self.modules:

                    num_nodes = len(module.nodes)

                    mat = np.empty(shape=(num_nodes, num_nodes))

                    for i, ny in enumerate(module.nodes):
                        for j, nx in enumerate(module.nodes):

                            quants = None
                            if i > j:
                                edges = nx.connects(ny)
                            else:
                                edges = ny.connects(nx)

                            # if len(edges) > 1:
                            #     self.log.warning("For Gene id %s ; Found multiple edges between nodes: %s" %
                            #                      (self.gene_id, str(edges)))


                            if filename == 'identity':
                                mat[i][j] = 1 if edges else 0
                            else:
                                if edges:
                                    #quants = self.quantifications(module, edge=edges[0])
                                    quant = self.edge_quant(module, edges[0], filename)
                                    mat[i][j] = quant
                                else:
                                    mat[i][j] = 0

                    dset = hf.create_dataset('%s_%s' % (self.gene_id, module.idx), data=mat)
                    dset.attrs['exons'] = " ".join((self._get_node_id(n.start, n.end) for n in module.nodes))


    def combine_hdf5s(self, dest_path, file_paths):
        with h5py.File(dest_path, "w") as out_hf:
            for hdf5_file in file_paths:
                with h5py.File(hdf5_file, "r") as in_hf:
                    for gene_id in in_hf['/']:
                        #print(gene_id)
                        gene_id = gene_id.encode('utf-8')
                        h5py.h5o.copy(in_hf.id, gene_id, out_hf.id, gene_id)
                os.remove(hdf5_file)




    def _all_junctions(self):
        """
        For each specific junction, find all possible exon combinations for each side
        (all path splits on that exon which have a more 'outer' coordinate (less on E1 for + strand,
        greater on E2 for - strand) Make a new row for each of these possibilities.
        :return:
        """
        junctions = set()
        for module in self.modules:
            for junc in module.get_all_edges(ir=True):
                e1 = self.graph.start_node(junc)
                e2 = self.graph.end_node(junc)

                # gather all possible points where a junction attaches ... do we care about fwd cases on e1 or just back?

                if self.graph.strand == '+':
                    e1_splice_sites = [j.end for j in e1.back_edges if j.end < junc.start]
                else:
                    e1_splice_sites = [j.end for j in e1.back_edges if j.end > junc.start]
                # same for e2, seems we may only care about fwd edges, not back, but verify
                if self.graph.strand == '+':
                    e2_splice_sites = [j.start for j in e2.edges if j.start > junc.end]
                else:
                    e2_splice_sites = [j.start for j in e2.edges if j.start < junc.end]

                # if the outermost edge of the exon is missing, add it
                if not e1.start in e1_splice_sites:
                    e1_splice_sites.append(e1.start)
                if not e2.end in e2_splice_sites:
                    e2_splice_sites.append(e2.end)


                # 'exponential' part, find all combinations of exons that can exist from combinations of e1 and e2
                for e1ss in e1_splice_sites:
                    for e2ss in e2_splice_sites:
                        # module_idx, e1_start, e1_end, e2_start, e2_end, ir
                        junc_row = (module.idx, e1ss, junc.start, junc.end, e2ss, junc.ir)
                        # if junc_row in junctions:
                        #     print("dup foind: " + str(junc_row))
                        junctions.add(junc_row)
                        # if junc.start == 228147699 and junc.end == 228148371:
                        #     print(junc_row)


        return junctions
                #     print(e1, e1_splice_sites)
                #     print(e2, e2_splice_sites)






    def _split_exons(self):
        """
        <DEFUNCT METHOD, PRESERVED IF WE MIGHT ANALYZE IT FOR RE-PURPOSE LATER>
        For training data, we need to make additional nodes for each time that there are two junctions
        in different positions in one exon, basically, any not at the ends after trimming. (alt3/5 ish)
        """
        for module in self.modules:

            # find non intronic edged that are not at the start or end of each exon (somewhere in the middle)
            # for i, node in enumerate(self.nodes):
            #     dupes_to_create = []
            #     for edge in node.edges:
            #         if not edge.ir and not edge.start == node.end:
            #             dupes_to_create.append(edge)
            #     for edge in node.back_edges:
            #         if not edge.ir and not edge.end == node.start:
            #             dupes_to_create.append(edge)
            #
            #     for new_exon in dupes_to_create:
            #         # for each of the found edges, we need to remove all dupe edges
            #         # clone the exon that many times, and add one of the dups junctions
            #         # if it is a back junction, we need to remove the current back junction


            for n1, n2 in combinations(module.nodes, 2):
                # look for cases with multiple connections
                fwd_connects = n1.connects(n2, ir=True)
                if len(fwd_connects) > 1:
                    # print(n1, n2)
                    # for all connections not the outermost, clone the node, remove all connections except that one,
                    # and trim the exon to it
                    for junc in fwd_connects:

                        if not junc.start == n1.end:
                            print('1 dup foind', junc, n1, n2)
                            #self._add_exon({'start': n1.start, 'end': junc.start})
                            dupe = self.graph.Node({'start': n1.start, 'end': junc.start})

                            if hasattr(n1, 'num_times_split'):
                                n1.num_times_split += 1
                            else:
                                n1.num_times_split = 1
                            dupe.split_index = n1.num_times_split + 1
                            dupe.maximal_start = n1.start
                            dupe.maximal_end = n1.end

                            junc.node = n2

                            junc.start_node = dupe

                            #print(junc)
                            #dupe.end = junc.start
                            dupe.edges = [junc]

                            for edge in n1.back_edges:
                                new_edge = self.graph.Edge({'start':edge.start, 'end':edge.end})
                                new_edge.node = dupe
                                new_edge.lsvs = edge.lsvs

                                new_edge.end_node = dupe

                                #index = self.graph.edges.index(edge)
                                #self.graph.edges.insert(index, edge)
                                dupe.back_edges.append(edge)
                                # need to get the node that the back edge connects to, and add the new edge to it
                                n0 = self.graph.start_node(edge)

                                n0.edges.append(new_edge)

                            #dupe.back_edges = n1.back_edges
                            dupe.idx = "%d_%d" % (dupe.start, dupe.end)
                            #dupes.append(dupe)
                            #print("added dupe fwd", n1.start, junc.start)
                            n1.edges.remove(junc)
                            #n2.back_edges.remove(junc)
                            # need to find the other end of all back edges, and make sure that they
                            # connect to the new junction as well

                            # need to append this node directly before / after the dupe node
                            index = module.nodes.index(n1)+1
                            module.nodes.insert(index, dupe)


                        if not junc.end == n2.start:
                            print('2 dup foind', junc, n1, n2)
                            #self._add_exon({'start': junc.end, 'end': n2.end})
                            #print("added dupe back", junc.end, n2.end)
                            dupe = self.graph.Node({'start': junc.end, 'end': n2.end})

                            if hasattr(n2, 'num_times_split'):
                                n2.num_times_split += 1
                            else:
                                n2.num_times_split = 1
                            dupe.split_index = n2.num_times_split + 1
                            dupe.maximal_start = n2.start
                            dupe.maximal_end = n2.end

                            junc.node = dupe
                            #print(junc, n2, dupe)
                            #junc.start_node = n2
                            junc.end_node = dupe



                            #print(junc, junc.node)
                            #dupe.end = junc.start
                            dupe.back_edges = [junc]
                            for edge in n2.edges:
                                new_edge = self.graph.Edge({'start':edge.start, 'end':edge.end})
                                new_edge.lsvs = edge.lsvs

                                new_edge.start_node = dupe

                                #index = self.graph.edges.index(edge)
                                #self.graph.edges.insert(index, edge)
                                dupe.edges.append(new_edge)
                                # need to get the node that the back edge connects to, and add the new edge to it
                                n0 = self.graph.end_node(edge)

                                n0.back_edges.append(new_edge)
                                new_edge.node = n0

                            #dupe.edges = n2.edges
                            dupe.idx = "%d_%d" % (dupe.start, dupe.end)
                            #dupes.append(dupe)
                            #n1.edges.remove(junc)
                            n2.back_edges.remove(junc)

                            index = module.nodes.index(n1)+1
                            module.nodes.insert(index, dupe)


