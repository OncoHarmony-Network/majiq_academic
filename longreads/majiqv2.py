from rna_voila.api import SpliceGraph
from rna_voila.api.view_splice_graph import ViewSpliceGraph
from graph import Graph, exon
from collections import namedtuple

fake_config = namedtuple('ClassifyConfig', "splice_graph_file")

class MajiqV2Reader:

    def __init__(self, path):
        self.path = path


    def parse_splicegraph(self, gene_id):
        with SpliceGraph(self.path) as sg:
            experiment_names = sg.experiment_names
            strand = sg.gene(gene_id)['strand']
        # with ViewSpliceGraph(path, splice_graph_file=path) as sg:
        #     self.view_gene = sg.view_gene(gene_id)
        #     print(self.view_gene)
            #gene_meta = sg.gene(gene_id)
            #print(gene_meta)

        self.meta = {'strand': strand, 'gene_id': gene_id, 'transcript_id': None}
        #config = fake_config(splice_graph_file=path)

        self.graph = Graph(gene_id, experiment_names, self.path)
        self.modules = self.graph.modules()



    def extent(self, gene_id):
        with ViewSpliceGraph(splice_graph_file=self.path) as sg:
            return sg.gene_start(gene_id), sg.gene_end(gene_id)

    def moduleExtent(self, module_idx):
        # must already have run parse_splicegraph()

        # we treat the module as
        start = float('inf')
        end = float('-inf')
        for edge in self.modules[module_idx].get_all_edges(ir=True):
            start = min(edge.start, start)
            end = max(edge.end, end)

        return start, end


    # def paths_through_gene(self):
    #     """
    #     For each module, recursively find paths out of the last exon entered in that path and add then to the running
    #     list of junctions in the path until we reach the last exon in the module. Then output an ordered list of
    #     paths that were found in each iteration.
    #
    #     This will end up skipping alt-first and "p-alt-first" because there are no junctions entering
    #     So we need to manually run the algorithm starting at these nodes as well as the first node.
    #     :return:
    #     """
    #
    #
    #     paths_found = []
    #     for module in self.modules:
    #         path_idx = 0
    #         #print('=====================')
    #
    #         def add_junctions_out(node, prev_juncs, exon_length, is_module_length):
    #             nonlocal path_idx
    #             #print('iter', node, prev_juncs, node.edges, node == module.nodes[-1])
    #             if not node.edges or node == module.nodes[-1]:
    #                 if (not node.edges and node != module.nodes[-1]) or node.end == -1:
    #                     # takes care of checking for ALE cases and p_ALE cases (they will always be at the end)
    #                     is_module_length = False
    #                 # got to end of possible path
    #                 #print(node)
    #                 exon_length = exon_length + (node.end - prev_juncs[-1].end) + 1
    #                 if is_module_length:
    #                     frameshift = str(exon_length % 3)
    #                 else:
    #                     frameshift = 'N/A'
    #                 path_idx += 1
    #                 paths_found.append((module.idx, path_idx, prev_juncs, frameshift))
    #             else:
    #                 for junc in node.edges:
    #                     _prev_juncs = prev_juncs[:]
    #                     _prev_juncs.append(junc)
    #                     next_node = self.graph.end_node(junc)
    #                     add_junctions_out(next_node, _prev_juncs,
    #                                       exon_length + (junc.start - node.start) - (junc.end - next_node.start) + 1,
    #                                       is_module_length)
    #
    #
    #         # run over junctions from the start of the module
    #         add_junctions_out(module.nodes[0], [], 0, True)
    #         # find other possible starting nodes like afe, p_afe, and run over those
    #         for node in module.nodes[1:-1]:
    #             if not node.back_edges:
    #                 add_junctions_out(node, [], 0, False)
    #
    #     return paths_found

    def getAllPathsUtil(self, u, d, visited, path, start=None, dont_append=False, is_denovo=False, has_reads=True):

        visited[u._idx]= True
        if start is None:
            start = u.start
        if not dont_append:
            path.append(exon(start, u.end))
        if u == d:

            yield path, is_denovo, has_reads
        else:

            for edge in u.edges:
                is_denovo = is_denovo or edge.is_de_novo()
                has_reads = has_reads and edge.has_reads
                _path = path[:]
                next_node = self.graph.end_node(edge)
                if edge.ir:
                    _path[-1] = exon(path[-1].start, next_node.end)
                    _start = next_node.end
                    yield from self.getAllPathsUtil(next_node, d, visited, _path, _start, True, is_denovo, has_reads)

                else:
                    _path[-1] = exon(path[-1].start, edge.start)
                    _start = edge.end
                    if visited[next_node._idx] == False:
                        yield from self.getAllPathsUtil(next_node, d, visited, _path, _start, False, is_denovo, has_reads)
                        #yield from self.getAllPathsUtil(next_node, d, visited, path)

        path.pop()
        visited[u._idx]= False


    # Prints all paths from 's' to 'd'
    def getAllPathsBetweenNodes(self, startNode, endNode):

        # Mark all the vertices as not visited
        visited =[False]*(len(self.graph.nodes))
        for i, node in enumerate(self.graph.nodes):
            node._idx = i

        # Create an array to store paths
        path = []

        # Call the recursive helper function to print all paths
        yield from self.getAllPathsUtil(startNode, endNode, visited, path)

    def getNumModules(self):
        return len(self.modules)

    def getAllPaths(self, module_idx=None):

        if module_idx is None:
            start_node_idx = 0
            end_node_idx = len(self.graph.nodes) - 1
        else:
            start_node_idx = self.modules[module_idx]._global_node_start_idx
            end_node_idx = self.modules[module_idx]._global_node_end_idx

        paths2search = [(self.graph.nodes[start_node_idx], self.graph.nodes[end_node_idx])]

        alt_starts = []
        alt_ends = []
        for node in self.graph.nodes[start_node_idx:end_node_idx]:
            if not node.edges and not node == self.graph.nodes[-1]:
                alt_ends.append(node)
            elif not node.back_edges and not node == self.graph.nodes[0]:
                alt_starts.append(node)

        # gather all possibilities from beginning to alt ends
        for alt_end in alt_ends:
            paths2search.append((self.graph.nodes[start_node_idx], alt_end))
        # # gather all possibilities from alt starts to end
        for alt_start in alt_starts:
            paths2search.append((alt_start, self.graph.nodes[end_node_idx]))
        # # finally, gather all possibilities from alt starts to alt ends
        for alt_start in alt_starts:
            for alt_end in alt_ends:
                paths2search.append((alt_start, alt_end))

        for start, end in paths2search:
            for path, is_denovo, has_reads in self.getAllPathsBetweenNodes(start, end):
                exons = tuple(exon(n.start, n.end) for n in path)
                yield exons, self.meta, is_denovo, has_reads


if __name__ == "__main__":
    sqlpath = '/tmp/sg_generated.sql'
    # parser = MajiqV2Reader(sqlpath)
    # parser.parse_splicegraph('generated')

    sqlpath = '/slowdata/lrdata/majiq/splicegraph.sql'
    parser = MajiqV2Reader(sqlpath)
    parser.parse_splicegraph("gene:ENSG00000109534")

    # for i in range(3):
    #     print(parser.modules[i].nodes)
    #     print(parser.modules[i]._global_node_start_idx, parser.modules[i]._global_node_end_idx)

    #print(parser.getNumModules())
    for path in parser.getAllPaths(module_idx=1):
        print(path)