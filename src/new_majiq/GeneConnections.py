"""
GeneConnections.py

Parent class that wraps new_majiq internals for gene connections (introns,
junctions)

Author: Joseph K Aicher
"""

from typing import Optional

import numpy as np

from new_majiq.Exons import Exons, _Exons
from new_majiq.GeneRegions import GeneRegions


class GeneConnections(GeneRegions):
    def __init__(self, gene_connections):
        super().__init__(gene_connections)
        return

    @property
    def _gene_connections(self):
        """Underlying internals class for gene connections"""
        return self._gene_regions

    def checksum(self):
        """Checksum including passed/simplified/connections status"""
        return self._gene_connections.checksum()

    def checksum_nodata(self):
        """Checksum only considering gene_idx, start, and end"""
        return self._gene_connections.checksum_nodata()

    def connect_exons(self, exons: Exons) -> None:
        """Connect regions to specified exons"""
        self._gene_connections.connect_exons(exons._exons)
        return

    @property
    def connected_exons(self) -> Optional[Exons]:
        """exons the connections are associated with (or None otherwise)"""
        raw: Optional[_Exons] = self._gene_connections.connected_exons
        return None if raw is None else Exons(raw)

    @property
    def denovo(self) -> np.ndarray:
        """Indicate if each connection is denovo or not"""
        return self._gene_connections.denovo

    @property
    def passed_build(self) -> np.ndarray:
        """Indicate if each connection passed build filters (reliable) or not"""
        return self._gene_connections.passed_build

    @property
    def simplified(self) -> np.ndarray:
        """Indicate if each connection is simplified or not"""
        return self._gene_connections.simplified

    @property
    def start_exon_idx(self) -> np.ndarray:
        """Indicate exon_idx associated with start coordinate"""
        return self._gene_connections.start_exon_idx

    @property
    def end_exon_idx(self) -> np.ndarray:
        """Indicate exon_idx associated with end coordinate"""
        return self._gene_connections.end_exon_idx

    def src_exon_idx(self, region_idx: Optional[np.ndarray] = None):
        if region_idx is None:
            region_idx = self._region_idx
        return self._gene_connections.src_exon_idx(region_idx)

    def dst_exon_idx(self, region_idx: Optional[np.ndarray] = None):
        if region_idx is None:
            region_idx = self._region_idx
        return self._gene_connections.dst_exon_idx(region_idx)

    def _pass_all(self) -> None:
        """Set all connections to have passed build"""
        self._gene_connections._pass_all()
        return

    def _simplify_all(self) -> None:
        """Set all connections to be simplified"""
        self._gene_connections._simplify_all()
        return

    def _unsimplify_all(self) -> None:
        """Set all connections to be unsimplified"""
        self._gene_connections._unsimplify_all()
        return
