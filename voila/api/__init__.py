from voila.api.matrix_hdf5 import Psi, DeltaPsi, Heterogen
from voila.api.splice_graph_sql import Exons, Junctions, Genes, IntronRetention


class SpliceGraph(Genes, Junctions, Exons, IntronRetention):
    pass


class Matrix(DeltaPsi, Psi, Heterogen):
    pass
