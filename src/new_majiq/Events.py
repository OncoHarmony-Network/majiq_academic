"""
Events.py

Wrapper around events

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.internals import Events as _Events

from pathlib import Path
from typing import (
    Final,
    Optional,
    Union,
)


class Events(object):
    """events and event connections"""

    def __init__(self, events: _Events):
        self._events: Final[_Events] = events
        return

    @property
    def introns(self) -> GeneIntrons:
        """Introns referenced by event connections"""
        return GeneIntrons(self._events.introns)

    @property
    def junctions(self) -> GeneJunctions:
        """Junctions referenced by event connections"""
        return GeneJunctions(self._events.junctions)

    @property
    def num_events(self) -> int:
        return self._events.num_events

    @property
    def e_idx(self) -> np.ndarray:
        return np.arange(self.num_events)

    @property
    def ref_exon_idx(self) -> np.ndarray:
        return self._events.ref_exon_idx

    @property
    def event_type(self) -> np.ndarray:
        return self._events.event_type

    @property
    def _offsets(self) -> np.ndarray:
        return self._events._offsets

    @property
    def ec_idx_start(self) -> np.ndarray:
        return self._events.connection_idx_start

    @property
    def ec_idx_end(self) -> np.ndarray:
        return self._events.connection_idx_end

    def connections_slice_for_event(self, event_idx: int) -> slice:
        return slice(
            self.connection_idx_start[event_idx],
            self.connection_idx_end[event_idx],
        )

    @property
    def num_connections(self) -> int:
        return self._events.num_connections

    @property
    def ec_idx(self) -> np.ndarray:
        return np.arange(self.num_connections)

    @property
    def is_intron(self) -> np.ndarray:
        return self._events.is_intron

    @property
    def connection_idx(self) -> np.ndarray:
        return self._events.idx

    def connection_gene_idx(self, ec_idx: Optional[np.ndarray] = None) -> np.ndarray:
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_gene_idx(ec_idx)

    def connection_start(self, ec_idx: Optional[np.ndarray] = None) -> np.ndarray:
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_start(ec_idx)

    def connection_end(self, ec_idx: Optional[np.ndarray] = None) -> np.ndarray:
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_end(ec_idx)

    @property
    def df(self) -> xr.Dataset:
        return xr.Dataset(
            {},
            {
                # events
                "e_idx": self.e_idx,
                "ref_exon_idx": ("e_idx", self.ref_exon_idx),
                "event_type": ("e_idx", self.event_type),
                # nice offsets
                "ec_idx_start": ("e_idx", self.ec_idx_start),
                "ec_idx_end": ("e_idx", self.ec_idx_end),
                # event connections
                "ec_idx": self.ec_idx,
                "is_intron": ("ec_idx", self.is_intron),
                "connection_idx": ("ec_idx", self.connection_idx),
            },
        )

    @property
    def df_events(self) -> xr.Dataset:
        return xr.Dataset(
            {},
            {
                "e_idx": self.e_idx,
                "ref_exon_idx": ("e_idx", self.ref_exon_idx),
                "event_type": ("e_idx", self.event_type),
                "ec_idx_start": ("e_idx", self.ec_idx_start),
                "ec_idx_end": ("e_idx", self.ec_idx_end),
            },
        )

    @property
    def df_event_connections(self) -> xr.Dataset:
        return xr.Dataset(
            {},
            {
                "ec_idx": self.ec_idx,
                "is_intron": ("ec_idx", self.is_intron),
                "connection_idx": ("ec_idx", self.connection_idx),
                "_offsets": ("e_offsets_idx", self._offsets),
            },
        )

    def to_netcdf(self, path: Union[str, Path], mode: str) -> None:
        """Save to specified netcdf file"""
        (
            self.df
            # drop indexes and nice offsets
            .drop_vars(["e_idx", "ec_idx_start", "ec_idx_end", "ec_idx"])
            # add raw offsets
            .assign_coords(_offsets=("e_offsets_idx", self._offsets))
            # add hash for introns/junctions
            .assign_attrs(
                intron_hash=self.introns.checksum(),
                junction_hash=self.junctions.checksum(),
            ).to_netcdf(path, mode, group=constants.NC_EVENTS)
        )
        return

    @classmethod
    def from_netcdf(
        cls,
        path: Union[str, Path],
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> "Events":
        df = xr.open_dataset(path, group=constants.NC_EVENTS)
        if df.intron_hash != introns.checksum():
            raise ValueError("Saved hash for introns does not match")
        if df.junction_hash != junctions.checksum():
            raise ValueError("Saved hash for junctions does not match")
        return Events(
            _Events(
                introns._gene_introns,
                junctions._gene_junctions,
                df.ref_exon_idx,
                df.event_type,
                df._offsets,
                df.is_intron,
                df.connection_idx,
            )
        )
