from voila.api.matrix_hdf5 import Psi, DeltaPsi, Heterogen
from voila.api.splice_graph_sqlite import Exons, Junctions, Genes, IntronRetentions


class SpliceGraph(Genes, Junctions, Exons, IntronRetentions):
    pass


class Matrix(DeltaPsi, Psi, Heterogen):
    pass
