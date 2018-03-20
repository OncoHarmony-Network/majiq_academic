from voila.api.matrix_hdf5 import Psi, DeltaPsi, Heterogen
from voila.api.splice_graph_sql import Exons, Junctions, Genes


class SpliceGraph(Genes, Junctions, Exons):
    pass


class Matrix(DeltaPsi, Psi, Heterogen):
    pass
