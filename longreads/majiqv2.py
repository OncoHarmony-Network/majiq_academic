from rna_voila.api import SpliceGraph
from rna_voila.api.view_splice_graph import ViewSpliceGraph
from graph import Graph, exon
from graph import module as _module
from collections import namedtuple

fake_config = namedtuple('ClassifyConfig', "splice_graph_file")

class MajiqV2Reader:

    def __init__(self, path):
        self.path = path
        with SpliceGraph(self.path) as sg:
            self.experiment_names = sg.experiment_names
            #strand = sg.gene(gene_id)['strand']

    @property
    def gene_ids(self):
        with SpliceGraph(self.path) as sg:
            for entry in sg.genes():
                yield entry['id']

    def has_gene(self, gene_id):
        with SpliceGraph(self.path) as sg:
            return sg.gene(gene_id) is not None

    def allpaths_data(self, max_paths, modules=None, module_idx=None, majiq_module_extent=None):
        majiq_exons = set()
        majiq_denovo = {}
        majiq_has_reads = {}

        num_paths = 0
        for (ord_majiq_transcript, majiq_meta, denovo, has_reads) in self.getAllPaths(module_idx=module_idx):
            num_paths += 1
            if max_paths == 0 or num_paths > max_paths:
                raise RecursionError()
            if modules:
                set_key = tuple(exon(max(majiq_module_extent[0], e.start) if e.start > 0 else e.start, min(majiq_module_extent[1], e.end) if e.end > 0 else e.end) for e in ord_majiq_transcript)
            else:
                set_key = tuple(exon(e.start, e.end) for e in ord_majiq_transcript)
            majiq_exons.add(set_key)
            majiq_denovo[set_key] = denovo
            majiq_has_reads[set_key] = has_reads

        return majiq_exons, majiq_denovo, majiq_has_reads

    def parse_splicegraph(self, gene_id):

        # with ViewSpliceGraph(path, splice_graph_file=path) as sg:
        #     self.view_gene = sg.view_gene(gene_id)
        #     print(self.view_gene)
            #gene_meta = sg.gene(gene_id)
            #print(gene_meta)

        #self.meta = {'strand': strand, 'gene_id': gene_id, 'transcript_id': None}
        #config = fake_config(splice_graph_file=path)

        self.graph = Graph(gene_id, self.experiment_names, self.path)
        self.modules = self.graph.modules()

    def annotated_starts(self, gene_id):
        with SpliceGraph(self.path) as sg:
            return set(x['coordinate'] for x in sg.alt_starts(gene_id))

    def annotated_ends(self, gene_id):
        with SpliceGraph(self.path) as sg:
            return set(x['coordinate'] for x in sg.alt_ends(gene_id))

    def annotated_exons(self, gene_id):
        starts = []
        ends = []
        with SpliceGraph(self.path) as sg:
            for x in sg.exons(gene_id):
                if x['annotated']:
                    starts.append(x['start'])
                    ends.append(x['end'])
        return starts, ends

    def all_exons(self, gene_id):
        starts = []
        ends = []
        with SpliceGraph(self.path) as sg:
            for x in sg.exons(gene_id):
                starts.append(x['start'])
                ends.append(x['end'])
        return starts, ends


    def annotated_exons_order(self, gene_id):
        all_coords_order = []
        with SpliceGraph(self.path) as sg:
            for x in sg.exons(gene_id):
                if x['annotated']:
                    all_coords_order.append(exon(x['start'], x['end']))
        return all_coords_order
        
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

        return _module(start, end)

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

                _path = path[:]
                next_node = self.graph.end_node(edge)
                if edge.ir:
                    _path[-1] = exon(path[-1].start, next_node.end)
                    _start = next_node.end
                    _is_denovo = is_denovo or edge.is_de_novo()
                    _has_reads = has_reads and edge.has_reads
                    yield from self.getAllPathsUtil(next_node, d, visited, _path, _start, True, _is_denovo, _has_reads)

                else:
                    _path[-1] = exon(path[-1].start, edge.start)
                    _start = edge.end
                    if visited[next_node._idx] == False:
                        _is_denovo = is_denovo or edge.is_de_novo()
                        _has_reads = has_reads and edge.has_reads
                        yield from self.getAllPathsUtil(next_node, d, visited, _path, _start, False, _is_denovo, _has_reads)
                        #yield from self.getAllPathsUtil(next_node, d, visited, path)

        path.pop()
        visited[u._idx] = False


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
        if len(self.graph.nodes) < 2:
            return
        if module_idx is None:
            start_node_idx = 0
            end_node_idx = len(self.graph.nodes) - 1
        else:
            start_node_idx = self.modules[module_idx]._global_node_start_idx
            end_node_idx = self.modules[module_idx]._global_node_end_idx

        #print(start_node_idx, end_node_idx, self.graph.nodes)
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
                #print("PATH", exons)
                yield exons, {}, is_denovo, has_reads


if __name__ == "__main__":
    sqlpath = '/tmp/sg_generated.sql'
    # parser = MajiqV2Reader(sqlpath)
    # parser.parse_splicegraph('generated')

    sqlpath = '/tmp/sg_generated.sql'
    parser = MajiqV2Reader(sqlpath)
    #parser.parse_splicegraph("gene:ENSG00000109534")

    gene_id = 'ENSG00000046651.16'
    gene_id = 'ENSG00000000003.15'
    gene_id = 'mxe_long_alt_last_minus'
    #print(list(parser.annotated_starts(gene_id)))

    # for i in range(3):
    #     print(parser.modules[i].nodes)
    #     print(parser.modules[i]._global_node_start_idx, parser.modules[i]._global_node_end_idx)
    parser.parse_splicegraph(gene_id)

    #print(parser.getNumModules())
    #print(len(list(parser.getAllPaths())))
    for path in parser.getAllPaths():
        print(path)