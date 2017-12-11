from voila.api import splice_graph_sql as sg, voila_sql as v


class SpliceGraph(sg.Genes, sg.Junctions, sg.Exons):
    pass


class Voila(v.Exons, v.Junctions, v.Lsvs):
    pass
