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
            self._split_exons()


    @staticmethod
    def tsv_names():
        names = ['exons.tsv', 'junctions.tsv']
        return names

    @staticmethod
    def delete_hdf5s():
        config = ClassifyConfig()
        paths = glob.glob(config.directory + '/*.hdf5.*')
        for path in paths:
            os.remove(path)

    def start_all_headers(self):

        if self.config.putative_multi_gene_regions:
            headers = ['Gene ID_Region', 'Gene ID', 'Gene Name', 'Chr', 'Strand', 'First Exon Start coord',
                       'First Exon End coord', 'Last Exon Start coord', "Last Exon End coord"]
            self.start_headers(headers, 'p_multi_gene_region.tsv')
            return

        headers = self.common_headers + ['Exon ID', 'Exon ID Start Coordinate', 'Exon ID End Coordinate']
        self.start_headers(headers, 'exons.tsv')
        headers = self.common_headers + ['Junc ID', 'Junc ID Start Coordinate', 'Junc ID End Coordinate',
                                         'Exon 1 ID', 'Exon 1 ID Start Coordinate', 'Exon 1 ID End Coordinate',
                                         'Exon 2 ID', 'Exon 2 ID Start Coordinate', 'Exon 2 ID End Coordinate']
        self.start_headers(headers, 'junctions.tsv')

    def exons_tsv(self):
        """
        exons.list format:
        Exon ID:  "gene_id"_"chr"_"strand"_"start"_"end".

        If there are two exon entries (an exon has two 5' or 3' splice sites) ,
         then the ID will have "gene_id"_"chr"_"strand"_"maximal-start"_"maximal-end"_1"
          and "gene_id"_"chr"_"strand"_"maximal-start"_"maximal-end"_2"
        """
        with open(os.path.join(self.config.directory, 'exons.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module in self.modules:

                common_data = ["%s_%d" % (self.gene_id, module.idx), self.gene_id, self.graph.gene_name,
                       self.graph.chromosome, self.graph.strand]

                rows = []
                for i, node in enumerate(module.nodes):



                    if hasattr(node, 'num_times_split'):
                        # should be the base node of a split, need to use extended ID
                        node_id = "{gene_id}_{chr}_{strand}_{max_start}_{max_end}_{idx}".format(
                            gene_id=self.gene_id,
                            chr=self.graph.chromosome,
                            strand=self.graph.strand,
                            max_start=node.start,
                            max_end=node.end,
                            idx=1
                        )
                    elif hasattr(node, 'split_index'):
                        # a shorter node made from a split, need to use extended ID
                        node_id = "{gene_id}_{chr}_{strand}_{max_start}_{max_end}_{idx}".format(
                            gene_id=self.gene_id,
                            chr=self.graph.chromosome,
                            strand=self.graph.strand,
                            max_start=node.maximal_start,
                            max_end=node.maximal_end,
                            idx=node.split_index
                        )
                    else:
                        node_id = "{gene_id}_{chr}_{strand}_{start}_{end}".format(
                            gene_id=self.gene_id,
                            chr=self.graph.chromosome,
                            strand=self.graph.strand,
                            start=node.start,
                            end=node.end
                        )
                    rows.append(common_data + [node_id, node.start, node.end])

                    # if module.graph.strand == '+':
                    #     rows = rows + tmp
                    # else:
                    #     rows = tmp + rows


                writer.writerows(rows)


    def junctions_tsv(self):
        with open(os.path.join(self.config.directory, 'junctions.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')
            for module in self.modules:

                common_data = ["%s_%d" % (self.gene_id, module.idx), self.gene_id, self.graph.gene_name,
                               self.graph.chromosome, self.graph.strand]

                rows = []

                for i, edge in enumerate(module.get_all_edges()):

                    # remove introns
                    if edge.ir:
                        continue

                    n1 = self.graph.start_node(edge)
                    n2 = self.graph.end_node(edge)
                    junc_id = "{gene_id}_{chr}_{strand}_{start}_{end}".format(
                        gene_id=self.gene_id,
                        chr=self.graph.chromosome,
                        strand=self.graph.strand,
                        start=edge.start,
                        end=edge.end
                    )

                    rows.append(common_data + [junc_id, edge.start, edge.end,
                                         n1.idx, n1.start, n1.end,
                                         n2.idx, n2.start, n2.end])

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
                    dset.attrs['exons'] = " ".join((n.idx for n in module.nodes))


    def combine_hdf5s(self, dest_path, file_paths):
        with h5py.File(dest_path, "w") as out_hf:
            for hdf5_file in file_paths:
                with h5py.File(hdf5_file, "r") as in_hf:
                    for gene_id in in_hf['/']:
                        #print(gene_id)
                        gene_id = gene_id.encode('utf-8')
                        h5py.h5o.copy(in_hf.id, gene_id, out_hf.id, gene_id)
                os.remove(hdf5_file)

    def _split_exons(self):
        """
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
                fwd_connects = n1.connects(n2, ir=False)
                if len(fwd_connects) > 1:
                    # print(n1, n2)
                    # for all connections not the outermost, clone the node, remove all connections except that one,
                    # and trim the exon to it
                    for junc in fwd_connects:
                        if not junc.start == n1.end:
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
                            #print(junc)
                            #dupe.end = junc.start
                            dupe.edges = [junc]
                            for edge in n1.back_edges:
                                new_edge = self.graph.Edge({'start':edge.start, 'end':edge.end})
                                new_edge.node = dupe
                                new_edge.lsvs = edge.lsvs
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
                            #print(junc, junc.node)
                            #dupe.end = junc.start
                            dupe.back_edges = [junc]
                            for edge in n2.edges:
                                new_edge = self.graph.Edge({'start':edge.start, 'end':edge.end})
                                new_edge.lsvs = edge.lsvs
                                #index = self.graph.edges.index(edge)
                                #self.graph.edges.insert(index, edge)
                                dupe.edges.append(new_edge)
                                # need to get the node that the back edge connects to, and add the new edge to it
                                n0 = self.graph.end_node(edge)

                                n0.back_edges.append(new_edge)
                                new_edge.node = n0

                            dupe.edges = n2.edges
                            dupe.idx = "%d_%d" % (dupe.start, dupe.end)
                            #dupes.append(dupe)
                            #n1.edges.remove(junc)
                            n2.back_edges.remove(junc)

                            index = module.nodes.index(n1)+1
                            module.nodes.insert(index, dupe)


