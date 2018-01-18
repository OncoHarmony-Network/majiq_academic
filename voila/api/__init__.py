from voila.api.matrix_hdf5 import Psi, DeltaPsi
from voila.api.splice_graph_sql import Exons, Junctions, Genes
from voila.api.voila_hdf5 import VoilaHDF5


class SpliceGraph(Genes, Junctions, Exons):
    pass


class Voila(VoilaHDF5):
    pass


class Matrix(DeltaPsi, Psi):
    pass
