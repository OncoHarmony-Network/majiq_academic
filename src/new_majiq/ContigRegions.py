"""
ContigRegions.py

Parent class that wraps new_majiq internals for contig regions

Author: Joseph K Aicher
"""

import numpy as np

from new_majiq.Regions import Regions
from new_majiq.Contigs import Contigs


class ContigRegions(Regions):

    def __init__(self, contig_regions):
        super().__init__(contig_regions)
        return

    @property
    def _contig_regions(self):
        """ Underlying internals class for contig regions
        """
        return self._regions

    @property
    def contigs(self) -> Contigs:
        """ Contigs for which these regions are defined
        """
        return Contigs(self._parents)

    @property
    def contig_idx(self) -> np.ndarray:
        """ Index of contig on which each region is defined
        """
        return self._contig_regions.contig_idx

    @property
    def strand(self) -> np.ndarray:
        """ Strand direction for each region
        """
        return self._contig_regions.strand

    def slice_for_contig(self, contig_idx: int) -> slice:
        return self._slice_for_parent(contig_idx)
