from rna_voila.api.matrix_hdf5 import Psi, DeltaPsi, Heterogen
from rna_voila.api.splice_graph import Exons, Junctions, Genes, IntronRetentions, AltStarts, AltEnds


class SpliceGraph(Genes, Junctions, Exons, IntronRetentions, AltStarts, AltEnds):
    pass


class Matrix(DeltaPsi, Psi, Heterogen):
    pass
