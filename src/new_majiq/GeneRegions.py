"""
ContigRegions.py

Parent class that wraps new_majiq internals for gene regions

Author: Joseph K Aicher
"""

import numpy as np

from new_majiq.Regions import Regions
from new_majiq.Genes import Genes


class GeneRegions(Regions):

    def __init__(self, gene_regions):
        super().__init__(gene_regions)
        return

    @property
    def _gene_regions(self):
        """ Underlying internals class for gene regions
        """
        return self._regions

    @property
    def genes(self) -> Genes:
        """ Genes for which these regions are defined
        """
        return Genes(self._parents)

    @property
    def gene_idx(self) -> np.ndarray:
        """ Index of gene on which region is defined
        """
        return self._gene_regions.gene_idx

    def slice_for_gene(self, gene_idx: int) -> slice:
        return self._slice_for_parent(gene_idx)
