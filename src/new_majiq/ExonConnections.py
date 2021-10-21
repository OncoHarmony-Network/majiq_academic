"""
ExonConnections.py

ExonConnections for SpliceGraph

Author: Joseph K Aicher
"""

from typing import Final, List

import numpy as np

import new_majiq.constants as constants
from new_majiq.Events import Events
from new_majiq.Exons import Exons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.internals import ExonConnections as _ExonConnections
from new_majiq.SimplifierGroup import SimplifierGroup, _SimplifierGroup


class ExonConnections(object):
    """Tracks introns/junctions associated with each exon (stranded direction)"""

    @classmethod
    def create_connecting(
        cls,
        exons: Exons,
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> "ExonConnections":
        """connect input exons, introns, and gene junctions"""
        return ExonConnections(
            _ExonConnections(
                exons._exons, introns._gene_introns, junctions._gene_junctions
            )
        )

    def __init__(self, exon_connections: _ExonConnections):
        self._exon_connections: Final[_ExonConnections] = exon_connections
        return

    def simplifier(self) -> SimplifierGroup:
        """Create simplifier group to unsimplify introns and junctions"""
        return SimplifierGroup(_SimplifierGroup(self._exon_connections))

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

    def lsvs(
        self, select_lsvs: constants.SelectLSVs = constants.DEFAULT_SELECT_LSVS
    ) -> Events:
        """construct Events for all LSVs defined by these exon connections"""
        if select_lsvs == constants.SelectLSVs.STRICT_LSVS:
            return Events(self._exon_connections.strict_lsvs())
        elif select_lsvs == constants.SelectLSVs.PERMISSIVE_LSVS:
            return Events(self._exon_connections.permissive_lsvs())
        elif select_lsvs == constants.SelectLSVs.SOURCE_LSVS:
            return Events(self._exon_connections.source_lsvs())
        elif select_lsvs == constants.SelectLSVs.TARGET_LSVS:
            return Events(self._exon_connections.target_lsvs())
        else:
            raise ValueError(
                f"Invalid {select_lsvs = }, must be from {list(constants.SelectLSVs)}"
            )

    def constitutive(self) -> Events:
        """construct Events for all constitutive events defined by ExonConnections"""
        return Events(self._exon_connections.constitutive())

    @staticmethod
    def _event_type_is_source(event_type: np.ndarray) -> np.ndarray:
        """convert array(dtype="S1") to array(dtype=bool) for vectorized internals"""
        event_type = np.array(event_type, copy=False)
        is_source = event_type == b"s"
        is_target = event_type == b"t"
        if not (is_source | is_target).all():
            raise ValueError("event_type has invalid values (must be b's' or b't')")
        return is_source

    def _events_for(self, ref_exon_idx: np.ndarray, event_type: np.ndarray) -> Events:
        """construct events for specified exons/event types"""
        return Events(
            self._exon_connections.events_for(
                ref_exon_idx, self._event_type_is_source(event_type)
            )
        )

    def has_intron(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> np.ndarray:
        """Indicate if selected events have a non-simplified intron"""
        return self._exon_connections.has_intron(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def event_size(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> np.ndarray:
        """Indicate number of connections in the event"""
        return self._exon_connections.event_size(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def passed(self, ref_exon_idx: np.ndarray, event_type: np.ndarray) -> np.ndarray:
        """Indicate if any of the connections in the event are passed"""
        return self._exon_connections.passed(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def redundant(self, ref_exon_idx: np.ndarray, event_type: np.ndarray) -> np.ndarray:
        """Indicate if the event is redundant (subset by a different event)"""
        return self._exon_connections.redundant(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def is_strict_LSV(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> np.ndarray:
        """Indicate if the event is a strict LSV

        (passed, event size > 1, nonredundant or mutually redundant source)
        """
        return self._exon_connections.is_strict_LSV(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def is_permissive_LSV(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> np.ndarray:
        """Indicate if the event is a permissive LSV

        (passed, event size > 1, not mutually redundant target)
        """
        return self._exon_connections.is_permissive_LSV(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def is_source_LSV(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> np.ndarray:
        """Indicate if the event is a source LSV

        (passed, event size > 1, event_type == 's')
        """
        return self._exon_connections.is_source_LSV(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def is_target_LSV(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> np.ndarray:
        """Indicate if the event is a target LSV

        (passed, event size > 1, event_type == 't')
        """
        return self._exon_connections.is_target_LSV(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def is_constitutive(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> np.ndarray:
        """Indicate if the event is constitutive (nonredundant with event_size == 1)"""
        return self._exon_connections.is_constitutive(
            ref_exon_idx, self._event_type_is_source(event_type)
        )

    def event_id(self, ref_exon_idx: np.ndarray, event_type: np.ndarray) -> List[str]:
        """List of event identifiers for VOILA for specified events"""
        return self._exon_connections.event_id(ref_exon_idx, event_type)

    def event_description(
        self, ref_exon_idx: np.ndarray, event_type: np.ndarray
    ) -> List[str]:
        """List of event descriptions for VOILA for specified events"""
        return self._exon_connections.event_description(ref_exon_idx, event_type)

    def src_introns_for(self, exon_idx: int) -> np.ndarray:
        """array of intron_idx that have exon_idx as src_exon

        Note
        ----
        This includes exitrons and simplified introns.
        """
        if exon_idx < 0 or exon_idx >= len(self.exons):
            raise ValueError("invalid exon_idx")
        idx_slice = slice(
            *self._exon_connections.src_intron_exon_offsets[exon_idx : (2 + exon_idx)]
        )
        return self._exon_connections.src_intron_idx[idx_slice]

    def dst_introns_for(self, exon_idx: int) -> np.ndarray:
        """array of intron_idx that have exon_idx as dst_exon

        Note
        ----
        This includes exitrons and simplified introns.
        """
        if exon_idx < 0 or exon_idx >= len(self.exons):
            raise ValueError("invalid exon_idx")
        idx_slice = slice(
            *self._exon_connections.dst_intron_exon_offsets[exon_idx : (2 + exon_idx)]
        )
        return self._exon_connections.dst_intron_idx[idx_slice]

    def src_junctions_for(self, exon_idx: int) -> np.ndarray:
        """array of junction_idx that have exon_idx as src_exon

        Note
        ----
        This includes exitrons and simplified junctions.
        """
        if exon_idx < 0 or exon_idx >= len(self.exons):
            raise ValueError("invalid exon_idx")
        idx_slice = slice(
            *self._exon_connections.src_junction_exon_offsets[exon_idx : (2 + exon_idx)]
        )
        return self._exon_connections.src_junction_idx[idx_slice]

    def dst_junctions_for(self, exon_idx: int) -> np.ndarray:
        """array of junction_idx that have exon_idx as dst_exon

        Note
        ----
        This includes exitrons and simplified junctions.
        """
        if exon_idx < 0 or exon_idx >= len(self.exons):
            raise ValueError("invalid exon_idx")
        idx_slice = slice(
            *self._exon_connections.dst_junction_exon_offsets[exon_idx : (2 + exon_idx)]
        )
        return self._exon_connections.dst_junction_idx[idx_slice]
