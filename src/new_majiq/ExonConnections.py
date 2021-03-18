"""
ExonConnections.py

ExonConnections for SpliceGraph

Author: Joseph K Aicher
"""

import numpy as np

from new_majiq.Exons import Exons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.Events import Events
from new_majiq.internals import ExonConnections as _ExonConnections

from typing import (
    Final,
    List,
)


class ExonConnections(object):
    """Tracks introns/junctions associated with each exon (stranded direction)"""

    def __init__(self, exon_connections: _ExonConnections):
        self._exon_connections: Final[_ExonConnections] = exon_connections
        return

    @property
    def exons(self) -> Exons:
        """underlying exons"""
        return Exons(self._exon_connections._exons)

    @property
    def introns(self) -> GeneIntrons:
        """underlying introns the exons are connected to"""
        return GeneIntrons(self._exon_connections._introns)

    @property
    def junctions(self) -> GeneJunctions:
        """underlying junctions the junctions are connected to"""
        return GeneJunctions(self._exon_connections._junctions)

    def lsvs(self) -> Events:
        """construct Events for all LSVs defined by these exon connections"""
        return Events(self._exon_connections.lsvs())

    def constitutive(self) -> Events:
        """construct Events for all constitutive events defined by ExonConnections"""
        return Events(self._exon_connections.constitutive())

    def event_id(self, ref_exon_idx: np.ndarray, event_type: np.ndarray) -> List[str]:
        """List of event identifiers for VOILA for specified events"""
        return self._exon_connections.event_id(ref_exon_idx, event_type)

    def event_description(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> List[str]:
        """List of event descriptions for VOILA for specified events"""
        return self._exon_connections.event_description(ref_exon_idx, event_type)
