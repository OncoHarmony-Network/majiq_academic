"""
GeneConnections.py

Parent class that wraps new_majiq internals for gene connections (introns,
junctions)

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

from new_majiq.GeneRegions import GeneRegions
from new_majiq.Exons import Exons, _Exons
from typing import (
    Optional,
    Union,
)
from pathlib import Path


class GeneConnections(GeneRegions):
    def __init__(self, gene_connections):
        super().__init__(gene_connections)
        return

    def checksum(self):
        """Checksum including passed/simplified/connections status"""
        return self._gene_introns.checksum()

    def checksum_nodata(self):
        """Checksum only considering gene_idx, start, and end"""
        return self._gene_introns.checksum_nodata()

    @property
    def _gene_connections(self):
        """Underlying internals class for gene connections"""
        return self._gene_regions

    def connect_exons(self, exons: Exons) -> None:
        """ Connect regions to specified exons
        """
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
