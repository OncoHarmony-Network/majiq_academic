from bisect import bisect_left, bisect_right
from itertools import combinations
from pathlib import Path

from voila import constants
from voila.api import SpliceGraph, Matrix
from voila.api.matrix_utils import generate_means
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import networkx.drawing.nx_pydot as nxpydot

PSI_THRESHOLD = 0.01
DPSI_THRESHOLD = None



class _Region:
    def __init__(self, start, end):
        self.start = start
        self.end   = end

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end
    def __lt__(self, other):
        return self.start < other.start or (self.start == other.start and self.end < other.end)
    def __repr__(self):
        return "%s-%s" % (self.start, self.end)
    def __hash__(self):
        return hash(self.__repr__())

class Gene(_Region):

    def __init__(self, gene_id, start, end):
        self.exon_list = list()
        self.junction_list = list()
        self.gene_id = gene_id
        self.module_list = list()
        super(Gene, self).__init__(start, end)

    def add_exon(self, ex):
        self.exon_list.append(ex)

    def add_junction(self, j):
        self.junction_list.append(j)


    def reconstruct_graph(self):
        num_exons = len(self.exon_list)
        ex_idx = 0

        for j in self.junction_list:
            for xdx in range(ex_idx, num_exons):
                ex = self.exon_list[xdx]
                if ex.start <= j.start <= ex.end :
                    ex.add_out_junc(j)
                    j.add_donor(ex)

        for j in self.junction_list:
            for xdx in range(ex_idx, num_exons):
                ex = self.exon_list[xdx]
                if ex.start <= j.end <= ex.end:
                    ex.add_in_junc(j)
                    j.add_acceptor(ex)


    def get_modules(self):

        for m in self.module_list:
            yield m

    def generate_module(self):

        open_j = 0
        start_ex = None


        for ex in self.exon_list:
            open_j -= len(ex.i_js)
            if open_j == 0:
                if start_ex is not None:
                    m = Module(start_ex.start, ex.end, start_ex, ex)
                    self.module_list.append(m)
                    # print("new_module")
                elif len(ex.i_js) == 1:
                    start_ex = None
                    continue
                start_ex = ex
            elif open_j<0:
                print('ERROR')

            open_j += len(ex.o_js)

class LSV():

    def __init__(self, id, type):
        self.id = id
        self.type = type

    def add_Exon(self, ex):
        self.exon = ex


class Exon(_Region):

    def __init__(self, start, end):
        self.o_js = set()
        self.i_js = set()
        self.slsv = None
        self.tlsv = None
        super(Exon, self).__init__(start, end)

    def add_in_junc(self, j):
        self.i_js.add(j)

    def add_out_junc(self, j):
        self.o_js.add(j)

    def add_lsv(self, lsvObj, is_source):
        if is_source:
            self.slsv = lsvObj
        else:
            self.tlsv = lsvObj

    def is_first(self):
        pass

    def is_last(self):
        pass


class Junction(_Region):
    def __init__(self, start, end, intron=False):
        self.donor = None
        self.acceptor = None
        self.psi = 1
        self.is_intron = intron
        super(Junction, self).__init__(start, end)

    def add_donor(self, ex):
        self.donor = ex

    def add_acceptor(self, ex):
        self.acceptor = ex

    def is_intron(self):
        return self.is_intron


class Module(_Region):

    def  __init__(self, start, end, init_exon, end_exon):
        self.start_ex = init_exon
        self.end_ex   = end_exon
        super(Module, self).__init__(start, end)

    def __repr__(self):
        s = str(self.start_ex) + ' -> ' + str(self.end_ex)
        return s

    def classify_module(self):
        pass



def sg_to_netwrkx(nid, start_node, end_node):

    dd = {}

    def recr_traverse(G, node, end_node, dd):
        # print(node)
        if node == end_node:
            return
        for j in node.o_js:
            if j.acceptor is not None and str(j) not in dd:
                if str(start_node) not in dd:
                    G.add_node(j.acceptor)
                    dd[str(start_node)] = 1
                dd[str(j)] = 1
                G.add_edge(node, j.acceptor)
                recr_traverse(G, j.acceptor, end_node, dd)
        return
    plt.figure()
    G = nx.MultiDiGraph()
    dd[str(start_node)] = 1
    G.add_node(start_node)
    recr_traverse(G, start_node, end_node, dd)
    nxpydot.write_dot(G, '%s.dot'% nid)


def from_splicegraph_to_gene(gObj, sg_file, lsv_dict={}):

    with SpliceGraph(sg_file) as sg:
        for exon in sg.exons(gObj.gene_id):

            key = '%s-%s' % (exon['start'], exon['end'])
            jst = exon['start'] if exon['start']>0 else exon['end'] - 1
            jnd = exon['end'] if exon['end'] > 0 else exon['start'] + 1
            ex = Exon(jst, jnd)
            gObj.add_exon(ex)

            try:
                sk = 's:' + key
                ex.add_lsv(lsv_dict[sk], True)
            except KeyError:
                pass

            try:
                sk = 't:' + key
                ex.add_lsv(lsv_dict[sk], False)
            except KeyError:
                pass


        for junc in sg.junctions(gObj.gene_id):
            jst = junc['start'] if junc['start']>0 else junc['end'] - 1
            jnd = junc['end'] if junc['end'] > 0 else junc['start'] + 1
            j = Junction(jst, jnd)
            gObj.add_junction(j)

        for junc in sg.intron_retentions(gObj.gene_id):
            j = Junction(junc['start'], junc['end'], True)
            gObj.add_junction(j)

    gObj.exon_list.sort()
    gObj.junction_list.sort()
    gObj.reconstruct_graph()

    sg_to_netwrkx(gObj.gene_id, gObj.exon_list[0], gObj.exon_list[-1])



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('sg_file', help='Splicegraph file that comes from the build execution')
    parser.add_argument('voila_file', help='voila file')
    args = parser.parse_args()

    lsv_d = dict()
    with Matrix(Path(args.voila_file).expanduser()) as m:
        gene_ids = list(m.gene_ids)

    for gene_id in gene_ids:
        with Matrix(Path(args.voila_file).expanduser()) as m:
            for lsv_id in m.lsv_ids(gene_ids=[gene_id]):
                lsv = m.psi(lsv_id)
                # print (lsv_id, lsv.__dict__)
                ex_coord = ':'.join(lsv_id.split(':')[1:])
                lsv_d[ex_coord] = LSV(lsv_id, lsv._lsv_type)


        print(gene_id)
        gObj = Gene(gene_id, 0, 10000)
        from_splicegraph_to_gene(gObj, args.sg_file, lsv_d)
        gObj.generate_module()

        for mm in gObj.get_modules():
            pass
            # sg_to_netwrkx(str(mm), mm.start_ex, mm.end_ex)

