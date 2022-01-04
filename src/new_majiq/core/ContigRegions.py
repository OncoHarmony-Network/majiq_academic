"""
ContigRegions.py

Parent class that wraps new_majiq internals for contig regions

Author: Joseph K Aicher
"""

import numpy as np
import numpy.typing as npt

from .Contigs import Contigs
from .Regions import Regions


class ContigRegions(Regions):
    def __init__(self, contig_regions):
        super().__init__(contig_regions)
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}[{len(self.contigs)} contigs, {len(self)}]"

    @property
    def _contig_regions(self):
        """Underlying internals class for contig regions"""
        return self._regions

    @property
    def contigs(self) -> Contigs:
        """Contigs for which these regions are defined"""
        return Contigs(self._parents)

    @property
    def contig_idx(self) -> npt.NDArray[np.uint64]:
        """Index of contig on which each region is defined"""
        return self._contig_regions.contig_idx

    @property
    def strand(self) -> npt.NDArray[np.str_]:
        """Strand direction for each region"""
        return self._contig_regions.strand

    def slice_for_contig(self, contig_idx: int) -> slice:
        return self._slice_for_parent(contig_idx)
