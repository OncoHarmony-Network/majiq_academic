"""
Events.py

Wrapper around events

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import Final, NamedTuple, Optional, Union

import numpy as np
import numpy.typing as npt
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
    """Masks betwen two :py:class:`Events` objects over e_idx for unique/shared events

    unique_events_mask: npt.NDArray[bool]
        boolean mask into events that are unique (i.e. not found in other)
    shared_events_idx: npt.NDArray[int]
        index into events in other for shared events (corresponding to False
        values in unique_events_mask)

    See Also
    --------
    Events.unique_events_mask
    """

    # boolean mask into events that are unique
    unique_events_mask: npt.NDArray[np.bool_]
    # index from nonunique to matching in other
    shared_events_idx: npt.NDArray[np.uint64]


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

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}["
            f"{self.num_events} events, {self.num_connections} connections]"
        )

    @property
    def introns(self) -> GeneIntrons:
        """:py:class:`Introns` over which events defined"""
        return GeneIntrons(self._events.introns)

    @property
    def junctions(self) -> GeneJunctions:
        """:py:class:`Junctions` over which events defined"""
        return GeneJunctions(self._events.junctions)

    @property
    def exons(self) -> Exons:
        """:py:class:`Exons` over which events defined"""
        return Exons(self._events.exons)

    @property
    def genes(self) -> Genes:
        """:py:class:`Genes` over which events defined"""
        return self.junctions.genes

    @property
    def contigs(self) -> Contigs:
        """:py:class:`Contigs` over which events defined"""
        return self.genes.contigs

    @property
    def num_events(self) -> int:
        """Number of events"""
        return self._events.num_events

    @property
    def e_idx(self) -> npt.NDArray[np.int64]:
        """Index over unique events"""
        return np.arange(self.num_events)

    @property
    def ref_exon_idx(self) -> npt.NDArray[np.uint64]:
        """Index into self.exons for reference exon of each unique event"""
        return self._events.ref_exon_idx

    @property
    def event_type(self) -> npt.NDArray[np.str_]:
        """Indicator if source ('s') or target ('b') for each unique event"""
        return self._events.event_type

    @property
    def _offsets(self) -> npt.NDArray[np.uint64]:
        """Offsets array for events into event connections"""
        return self._events._offsets

    @property
    def ec_idx_start(self) -> npt.NDArray[np.uint64]:
        """First index into event connections (ec_idx) for each unique event"""
        return self._events.connection_idx_start

    @property
    def ec_idx_end(self) -> npt.NDArray[np.uint64]:
        """One-past-end index into event connections (ec_idx) for each unique event"""
        return self._events.connection_idx_end

    def ec_idx_slice_for_event(self, e_idx: int) -> slice:
        """Get slice into event connections (ec_idx) for specified event"""
        return slice(self.ec_idx_start[e_idx], self.ec_idx_end[e_idx])

    @property
    def _gene_offsets(self) -> npt.NDArray[np.uint64]:
        """Offsets array for genes into events"""
        return self._events._gene_offsets

    @property
    def e_idx_start(self) -> npt.NDArray[np.uint64]:
        """First index into events (e_idx) for each gene"""
        return self._events.event_idx_start

    @property
    def e_idx_end(self) -> npt.NDArray[np.uint64]:
        """One-past-end index into events (e_idx) for each gene"""
        return self._events.event_idx_end

    def e_idx_slice_for_gene(self, gene_idx: int) -> slice:
        """Get slice into events (e_idx) for specified gene"""
        return slice(self.e_idx_start[gene_idx], self.e_idx_end[gene_idx])

    def slice_for_gene(self, gene_idx: int) -> slice:
        """Get slice into events (e_idx) for specified gene"""
        return self.e_idx_slice_for_gene(gene_idx)

    @cached_property
    def event_size(self) -> npt.NDArray[np.uint64]:
        """Number of event connections for each unique event"""
        return self.ec_idx_end - self.ec_idx_start

    def broadcast_eidx_to_ecidx(self, x: npt.ArrayLike, axis: int = 0) -> npt.NDArray:
        """Broadcast `x` over events to event connections

        Parameters
        ----------
        x: array_like
            Array of length self.num_events that will have values per-event
            repeated for each connection in the event
        axis: int
            Axis to broadcast from events to event connections (must have shape
            equal to `num_events`)

        Returns
        -------
        array
            with values of `x` repeated for each event connection
        """
        x = np.array(x, copy=False)
        try:
            if x.shape[axis] != self.num_events:
                raise ValueError("x must have length equal to the number of events")
        except IndexError:
            raise ValueError(f"x must have {axis = } to broadcast over")
        return np.take(x, self.connection_e_idx.view(np.int64), axis=axis)

    def connections_slice_for_event(self, event_idx: int) -> slice:
        """Get slice into event connections for event with specified index

        Parameters
        ----------
        event_idx: int
            Index of single event to get slice into event connections for

        Returns
        -------
        slice
        """
        return slice(
            self.ec_idx_start[event_idx],
            self.ec_idx_end[event_idx],
        )

    @property
    def num_connections(self) -> int:
        """Total number of connections over all events

        Total number of connections over all events (double counting if in
        source and target events)
        """
        return self._events.num_connections

    @property
    def ec_idx(self) -> npt.NDArray[np.int64]:
        """Index over event connections"""
        return np.arange(self.num_connections)

    @property
    def connection_e_idx(self) -> npt.NDArray[np.uint64]:
        """Index into events for each event connection"""
        return self._events.connection_event_idx

    @property
    def is_intron(self) -> npt.NDArray[np.bool_]:
        """Indicator if an intron or junction for each event connection"""
        return self._events.is_intron

    @property
    def connection_idx(self) -> npt.NDArray[np.uint64]:
        """Index into self.introns or self.junctions for each event connection"""
        return self._events.idx

    def connection_gene_idx(
        self, ec_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.uint64]:
        """Index into self.genes for selected event connections

        Parameters
        ----------
        ec_idx: Optional[npt._ArrayLikeInt_co]
            Indexes of selected event connections. If None, select all event
            connections in order

        Returns
        -------
        npt.NDArray[np.uint64]
        """
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_gene_idx(ec_idx)

    def connection_contig_idx(
        self, ec_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.uint64]:
        """Index into self.contigs for selected event connections

        Parameters
        ----------
        ec_idx: Optional[npt._ArrayLikeInt_co]
            Indexes of selected event connections. If None, select all event
            connections in order

        Returns
        -------
        npt.NDArray[np.uint64]
        """
        return self.genes.contig_idx[self.connection_gene_idx(ec_idx)]

    def connection_start(
        self, ec_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.int64]:
        """Start coordinate for each selected event connection

        Parameters
        ----------
        ec_idx: Optional[npt._ArrayLikeInt_co]
            Indexes of selected event connections. If None, select all event
            connections in order

        Returns
        -------
        npt.NDArray[np.int64]
        """
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_start(ec_idx)

    def connection_end(
        self, ec_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.int64]:
        """End coordinate for each selected event connection

        Parameters
        ----------
        ec_idx: Optional[npt._ArrayLikeInt_co]
            Indexes of selected event connections. If None, select all event
            connections in order

        Returns
        -------
        npt.NDArray[np.int64]
        """
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_end(ec_idx)

    def connection_denovo(
        self, ec_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.bool_]:
        """Indicator if connection was denovo for each selected event connection

        Parameters
        ----------
        ec_idx: Optional[npt._ArrayLikeInt_co]
            Indexes of selected event connections. If None, select all event
            connections in order

        Returns
        -------
        npt.NDArray[np.bool_]
        """
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_denovo(ec_idx)

    @property
    def connection_ref_exon_idx(self) -> npt.NDArray[np.uint64]:
        """Index into self.exons for reference exon for each event connection"""
        return self.broadcast_eidx_to_ecidx(self.ref_exon_idx)

    def connection_other_exon_idx(
        self, ec_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.uint64]:
        """Index into self.exons for nonreference exon for each event connection"""
        if ec_idx is None:
            ec_idx = self.ec_idx
        return self._events.connection_other_exon_idx(ec_idx)

    @property
    def df(self) -> xr.Dataset:
        """:py:class:`xr.Dataset` with event and event connections information"""
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
        """:py:class:`xr.Dataset` with event information"""
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
        """:py:class:`xr.Dataset` with event connections information"""
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
        """:py:class:`xr.Dataset` that is directly saved for Events"""
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
        """Save :py:class:`Events` to specified path"""
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
    def from_arrays(
        cls,
        introns: GeneIntrons,
        junctions: GeneJunctions,
        ref_exon_idx: npt._ArrayLikeInt_co,
        event_type: npt._ArrayLikeStr_co,
        offsets: npt._ArrayLikeInt_co,
        is_intron: npt._ArrayLikeBool_co,
        connection_idx: npt._ArrayLikeInt_co,
    ) -> "Events":
        """Create :class:`Events` from connections and input arrays

        Create :class:`Events` from connections (:class:`GeneIntrons` and
        :class:`GeneJunctions`) and input arrays
        """
        return Events(
            _Events(
                introns._gene_introns,
                junctions._gene_junctions,
                ref_exon_idx,
                event_type,
                offsets,
                is_intron,
                connection_idx,
            )
        )

    @classmethod
    def from_zarr(
        cls,
        path: Union[str, Path],
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> "Events":
        """Load :py:class:`Events` from specified path"""
        with xr.open_zarr(path, group=constants.NC_EVENTS) as df:
            df.load()
            if df.intron_hash != introns.checksum():
                raise ValueError("Saved hash for introns does not match")
            if df.junction_hash != junctions.checksum():
                raise ValueError("Saved hash for junctions does not match")
            return Events.from_arrays(
                introns,
                junctions,
                df.ref_exon_idx.values,
                df.event_type.values,
                df._offsets.values,
                df.is_intron.values,
                df.connection_idx.values,
            )

    def __getitem__(self, event_mask) -> "Events":
        """Subset :py:class:`Events` corresponding to boolean event mask

        Returns
        -------
        Events
        """
        event_mask = np.array(event_mask, copy=False, dtype=bool)  # make sure array
        if event_mask.ndim != 1:
            raise ValueError("event_mask must be 1-dimensional")
        elif len(event_mask) != self.num_events:
            raise ValueError("event_mask must match events")
        subset_size = self.event_size[event_mask]
        subset_offsets = np.empty(1 + len(subset_size), dtype=np.uint64)
        subset_offsets[0] = 0
        np.cumsum(subset_size, out=subset_offsets[1:])
        ec_mask = self.broadcast_eidx_to_ecidx(event_mask)
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
        """Get :py:class:`UniqueEventsMasks` with shared events and events unique to self

        Parameters
        ----------
        other: Events

        Returns
        -------
        UniqueEventsMasks
        """
        aligned = EventsAlign(self._events, other._events)
        unique_mask = np.ones(self.num_events, dtype=bool)
        unique_mask[aligned.left_event_idx] = False  # if aligned, not unique
        return UniqueEventsMasks(
            unique_events_mask=unique_mask,
            shared_events_idx=np.array(aligned.right_event_idx, copy=True),
        )

    @property
    def ec_dataframe(self) -> pd.DataFrame:
        """:py:class:`pd.DataFrame` over event connections detailing genomic information"""
        gene_idx = self.connection_gene_idx()
        other_exon_idx = self.connection_other_exon_idx()
        return pd.DataFrame(
            dict(
                seqid=np.array(self.contigs.seqid)[self.connection_contig_idx()],
                strand=self.genes.strand[gene_idx],
                gene_name=np.array(self.genes.gene_name)[gene_idx],
                gene_id=np.array(self.genes.gene_id)[gene_idx],
                event_type=self.broadcast_eidx_to_ecidx(self.event_type),
                ref_exon_start=self.exons.start[self.connection_ref_exon_idx],
                ref_exon_end=self.exons.end[self.connection_ref_exon_idx],
                start=self.connection_start(),
                end=self.connection_end(),
                is_denovo=self.connection_denovo(),
                is_intron=self.is_intron,
                other_exon_start=self.exons.start[other_exon_idx],
                other_exon_end=self.exons.end[other_exon_idx],
            ),
            index=pd.Index(self.ec_idx, name="ec_idx"),
        ).assign(
            event_type=lambda df: df.event_type.str.decode("utf-8"),
            strand=lambda df: df.strand.str.decode("utf-8"),
        )
