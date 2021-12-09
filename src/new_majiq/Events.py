"""
Events.py

Wrapper around events

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Final, NamedTuple, Optional, Union

import numpy as np
import pandas as pd
import xarray as xr

import new_majiq.constants as constants
from new_majiq.internals import Events as _Events
from new_majiq.internals import EventsAlign

from ._workarounds import _load_zerodim_variables
from .Contigs import Contigs
from .Exons import Exons
from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions
from .Genes import Genes


class UniqueEventsMasks(NamedTuple):
    """
    unique_events_mask: np.ndarray
        boolean mask into events that are unique (i.e. not found in other)
    shared_events_idx: np.ndarray
        index into events in other for shared events (corresponding to False
        values in unique_events_mask)
    """

    unique_events_mask: np.ndarray  # boolean mask into events that are unique
    shared_events_idx: np.ndarray  # index from nonunique to matching in other


class Events(object):
    """Collections of introns/junctions all starting or ending at the same exon

    Parameters
    ----------
    events: _Events
        Underlying object binding the internal C++ API

    See Also
    --------
    Events.from_zarr
    ExonConnections.lsvs
    ExonConnections.constitutive
    PsiCoverage.get_events
    PsiControlsSummary.get_events
    """

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
    def exons(self) -> Exons:
        """Exons used by these events"""
        exons = self.junctions.connected_exons
        assert isinstance(exons, Exons), "Events appears to not reference exons"
        return exons

    @property
    def genes(self) -> Genes:
        """Genes used by these events"""
        return self.junctions.genes

    @property
    def contigs(self) -> Contigs:
        """Contigs used by these events"""
        return self.genes.contigs

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
            self.ec_idx_start[event_idx],
            self.ec_idx_end[event_idx],
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

    def connection_contig_idx(self, ec_idx: Optional[np.ndarray] = None) -> np.ndarray:
        return self.genes.contig_idx[self.connection_gene_idx(ec_idx)]

    def connection_start(self, ec_idx: Optional[np.ndarray] = None) -> np.ndarray:
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_start(ec_idx)

    def connection_end(self, ec_idx: Optional[np.ndarray] = None) -> np.ndarray:
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_end(ec_idx)

    def connection_denovo(self, ec_idx: Optional[np.ndarray] = None) -> np.ndarray:
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_denovo(ec_idx)

    @property
    def connection_ref_exon_idx(self) -> np.ndarray:
        return np.repeat(self.ref_exon_idx, np.diff(self._offsets.astype(int)))

    def connection_other_exon_idx(
        self, ec_idx: Optional[np.ndarray] = None
    ) -> np.ndarray:
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_other_exon_idx(ec_idx)

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

    @property
    def save_df(self) -> xr.Dataset:
        """Dataset that is directly saved for Events"""
        return (
            self.df
            # drop indexes and nice offsets
            .drop_vars(["e_idx", "ec_idx_start", "ec_idx_end", "ec_idx"])
            # add raw offsets
            .assign_coords(_offsets=("e_offsets_idx", self._offsets))
            # add hash for introns/junctions
            .assign_attrs(
                intron_hash=self.introns.checksum(),
                junction_hash=self.junctions.checksum(),
            )
        )

    def to_zarr(
        self, path: Union[str, Path], mode: str, consolidated: bool = True
    ) -> None:
        """Save to specified zarr file"""
        self.save_df.pipe(lambda x: x.chunk(x.sizes)).pipe(
            _load_zerodim_variables
        ).to_zarr(
            path,
            mode=mode,
            group=constants.NC_EVENTS,
            consolidated=consolidated,
        )
        return

    @classmethod
    def from_zarr(
        cls,
        path: Union[str, Path],
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> "Events":
        with xr.open_zarr(path, group=constants.NC_EVENTS) as df:
            if df.intron_hash != introns.checksum():
                raise ValueError("Saved hash for introns does not match")
            if df.junction_hash != junctions.checksum():
                raise ValueError("Saved hash for junctions does not match")
            return Events(
                _Events(
                    introns._gene_introns,
                    junctions._gene_junctions,
                    df.ref_exon_idx.values,
                    df.event_type.values,
                    df._offsets.values,
                    df.is_intron.values,
                    df.connection_idx.values,
                )
            )

    def __getitem__(self, event_mask) -> "Events":
        """Subset of events corresponding to boolean event mask"""
        event_mask = np.array(event_mask, copy=False, dtype=bool)  # make sure array
        if event_mask.ndim != 1:
            raise ValueError("event_mask must be 1-dimensional")
        elif len(event_mask) != self.num_events:
            raise ValueError("event_mask must match events")
        event_size = (self.ec_idx_end - self.ec_idx_start).astype(int)
        subset_size = event_size[event_mask]
        subset_offsets = np.empty(1 + len(subset_size), dtype=np.uint64)
        subset_offsets[0] = 0
        subset_offsets[1:] = np.cumsum(subset_size)
        ec_mask = np.repeat(event_mask, event_size)
        return Events(
            _Events(
                self.introns._gene_introns,
                self.junctions._gene_junctions,
                self.ref_exon_idx[event_mask],
                self.event_type[event_mask],
                subset_offsets,
                self.is_intron[ec_mask],
                self.connection_idx[ec_mask],
            )
        )

    def unique_events_mask(self, other: "Events") -> UniqueEventsMasks:
        """Get events unique to self vs other, and indexes to shared events"""
        aligned = EventsAlign(self._events, other._events)
        unique_mask = np.ones(self.num_events, dtype=bool)
        unique_mask[aligned.left_event_idx] = False  # if aligned, not unique
        return UniqueEventsMasks(
            unique_events_mask=unique_mask,
            shared_events_idx=np.array(aligned.right_event_idx, copy=True),
        )

    @property
    def ec_dataframe(self) -> pd.DataFrame:
        """Annotated event connections as dataframe"""
        gene_idx = self.connection_gene_idx()
        other_exon_idx = self.connection_other_exon_idx()
        return pd.DataFrame(
            dict(
                seqid=np.array(self.contigs.seqid)[self.connection_contig_idx()],
                strand=self.genes.strand[gene_idx],
                gene_name=np.array(self.genes.gene_name)[gene_idx],
                gene_id=np.array(self.genes.gene_id)[gene_idx],
                event_type=np.repeat(
                    self.event_type, np.diff(self._offsets.astype(int))
                ),
                ref_exon_start=self.exons.start[self.connection_ref_exon_idx],
                ref_exon_end=self.exons.end[self.connection_ref_exon_idx],
                start=self.connection_start(),
                end=self.connection_end(),
                is_denovo=self.connection_denovo(),
                is_intron=self.is_intron,
                other_exon_start=self.exons.start[other_exon_idx],
                other_exon_end=self.exons.end[other_exon_idx],
            ),
            index=pd.Index(np.arange(self.num_connections), name="ec_idx"),
        ).assign(
            event_type=lambda df: df.event_type.str.decode("utf-8"),
            strand=lambda df: df.strand.str.decode("utf-8"),
        )
