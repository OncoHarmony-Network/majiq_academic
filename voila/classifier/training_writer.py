from voila.config import ClassifyConfig
import numpy as np
import h5py
from voila.classifier.tsv_writer import BaseTsvWriter
import os, csv

class TrainingWriter(BaseTsvWriter):

    def __init__(self, graph, gene_id):
        """

        :param output_path: The folder where all output TSV files will be written under
        :param graph: the Graph object of the gene
        """
        super().__init__(graph, gene_id)

        self.config = ClassifyConfig()

        self.graph = graph
        self.gene_id = gene_id

        if self.graph:
            self.modules = self.graph.modules()

    @staticmethod
    def tsv_names():
        names = ['exons.tsv', 'junctions.tsv']
        return names

    @staticmethod
    def delete_hdf5s():
        config = ClassifyConfig()
        path = os.path.join(config.directory, 'adjacency_matrix.hdf5')
        if os.path.exists(path):
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
        with open(os.path.join(self.config.directory, 'exons.tsv.%s' % self.pid), 'a', newline='') as csvfile:
            writer = csv.writer(csvfile, dialect='excel-tab', delimiter='\t')

            for module in self.modules:

                common_data = ["%s_%d" % (self.gene_id, module.idx), self.gene_id, self.graph.gene_name,
                       self.graph.chromosome, self.graph.strand]

                rows = []
                for i, node in enumerate(module.nodes):

                    if i != len(module.nodes)-1:
                        fwd_connections = [e for e in node.edges if not e.ir]

                    else:
                        fwd_connections = []
                    if i != 0:
                        back_connections = [e for e in node.back_edges if not e.ir]

                    else:
                        back_connections = []

                    # remove introns
                    all_connections = []

                    if fwd_connections:
                        fwd_coords = set()
                        for edge in fwd_connections:
                            # we need to find all connections that do not share a start coord
                            if not edge.start in fwd_coords:
                                fwd_coords.add(edge.start)
                                all_connections.append(edge)
                    if back_connections:
                        back_coords = set()
                        for edge in back_connections:
                            # we need to find all connections that do not share a start coord
                            if not edge.end in back_coords:
                                back_coords.add(edge.end)
                                all_connections.append(edge)

                    #tmp = []
                    for j, edge in enumerate(all_connections, start=1):
                        rows.append(common_data + ['%s.%s' % (node.idx,
                                                              j), node.start, node.end])

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

                    rows.append(common_data + [edge.idx, edge.start, edge.end,
                                         n1.idx, n1.start, n1.end,
                                         n2.idx, n2.start, n2.end])

                writer.writerows(rows)

    def _mat_identity(self, nx, ny):
        if nx.connects(ny) or ny.connects(nx):
            return 1
        return 0

    def _mat_psi1(self, nx, ny):
        self.quantifications()
        if nx.connects(ny) or ny.connects(nx):
            return 1
        return 0

    def _mat_psi2(self, nx, ny):
        if nx.connects(ny) or ny.connects(nx):
            return 1
        return 0

    def _mat_dpsi(self, nx, ny):
        if nx.connects(ny) or ny.connects(nx):
            return 1
        return 0

    def adjacency_matrix(self):


        modules = self.graph.modules()
        #as_types = {x.idx: x.as_types() for x in modules}



        # print([x for x in modules[0].nodes])
        # print([x for x in modules[0].get_all_edges()])


        with h5py.File(os.path.join(self.config.directory, 'adjacency_matrix.hdf5.%s' % self.pid), "a") as hf:

            for mat_func in [(self._mat_identity, 'identity',)]:

                for module in modules:
                    num_nodes = len(module.nodes)

                    mat = np.empty(shape=(num_nodes, num_nodes))

                    for i, ny in enumerate(module.nodes):
                        for j, nx in enumerate(module.nodes):

                            mat[i][j] = mat_func[0](nx, ny)

                    hf.create_dataset('%s_%s/%s' % (self.gene_id, module.idx, mat_func[1]), data=mat)


    def combine_hdf5s(self, file_paths):
        with h5py.File(os.path.join(self.config.directory, 'adjacency_matrix.hdf5'), "w") as out_hf:
            for hdf5_file in file_paths:
                with h5py.File(hdf5_file, "r") as in_hf:

                    for gene_id in in_hf['/']:
                        #print(gene_id)
                        gene_id = gene_id.encode('utf-8')
                        h5py.h5o.copy(in_hf.id, gene_id, out_hf.id, gene_id)
                os.remove(hdf5_file)



