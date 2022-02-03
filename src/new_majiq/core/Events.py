"""
Events.py

Wrapper around events

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import Final, List, NamedTuple, Optional, Sequence, Union

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

    def index(
        self, ref_exon_idx: npt._ArrayLikeInt_co, event_type: npt._ArrayLikeStr_co
    ) -> npt.NDArray[np.int64]:
        """Get index for specified event in self (-1 if not present)

        Raises ValueError if ref_exon_idx is outside of valid range
        (i.e. >= len(exons))
        """
        from .ExonConnections import ExonConnections

        is_source = ExonConnections._event_type_is_source(event_type)
        return self._events.index(ref_exon_idx, is_source)

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

    def event_id(self, e_idx: Optional[npt._ArrayLikeInt_co] = None) -> List[str]:
        """List of event identifiers for VOILA for specified events"""
        if e_idx is None:
            e_idx = self.e_idx
        return self._events.event_id(e_idx)

    def has_intron(
        self, e_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.bool_]:
        """Indicate if selected events have an intron"""
        if e_idx is None:
            e_idx = self.e_idx
        return self.is_intron[self.ec_idx_end.view(np.int64)[e_idx] - 1]

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

    def select_eidx_to_select_ecidx(
        self, select_eidx: npt._ArrayLikeInt_co
    ) -> npt.NDArray[np.int64]:
        """Given index array for e_idx, get index array for associated ec_idx

        Parameters
        ----------
        select_eidx: array[int]
            1D array of indexes into select events. Values must be between 0
            and `self.num_events`.

        Returns
        -------
        array[int]
            1D array of indexes into event connections associated with
            `select_eidx`

        Notes
        -----
        If you have a mask over events, use broadcast_eidx_to_ecidx instead
        (i.e. don't convert to select_eidx using np.where).

        See Also
        --------
        Events.connections_slice_for_event : for 0D/scalar case
        """
        # make sure select_eidx is an array
        select_eidx = np.array(select_eidx, copy=False)
        # make sure it is 1D
        if select_eidx.ndim != 1:
            raise ValueError("select_eidx must be 1D")
        # convert to mask over events
        select_eidx_mask = np.zeros(self.num_events, dtype=bool)
        select_eidx_mask[select_eidx] = True
        # broadcast to connections
        select_ecidx_mask = self.broadcast_eidx_to_ecidx(select_eidx_mask)
        # extract indexes from mask
        return np.where(select_ecidx_mask)[0]

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

    def ec_dataframe(
        self,
        annotated_exons: Optional[Exons] = None,
        annotated_introns: Optional[GeneIntrons] = None,
    ) -> pd.DataFrame:
        """:py:class:`pd.DataFrame` over event connections detailing genomic information

        Parameters
        ----------
        annotated_exons: Optional[Exons]
            If specified, use to reduce the number of exons called as denovo
        annotated_introns: Optional[GeneIntrons]
            If specified, use to reduce the number of introns called as denovo

        Returns
        -------
        pd.DataFrame
        """
        gene_idx = self.connection_gene_idx()
        other_exon_idx = self.connection_other_exon_idx()
        # get connection denovo status
        connection_denovo = self.connection_denovo()
        connection_denovo[self.is_intron] = self.introns.is_denovo(
            self.connection_idx[self.is_intron], annotated_introns=annotated_introns
        )
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
                ref_exon_denovo=self.exons.is_denovo(
                    self.connection_ref_exon_idx, annotated_exons=annotated_exons
                ),
                is_denovo=connection_denovo,
                is_intron=self.is_intron,
                other_exon_start=self.exons.start[other_exon_idx],
                other_exon_end=self.exons.end[other_exon_idx],
                other_exon_denovo=self.exons.is_denovo(
                    other_exon_idx, annotated_exons=annotated_exons
                ),
            ),
            index=pd.Index(self.ec_idx, name="ec_idx"),
        ).assign(
            event_type=lambda df: df.event_type.str.decode("utf-8"),
            strand=lambda df: df.strand.str.decode("utf-8"),
        )

    def merge_dataframes(
        self,
        df_seq: Sequence[pd.DataFrame],
        events_seq: Sequence["Events"],
        annotated_exons: Optional[Exons] = None,
        annotated_introns: Optional[GeneIntrons] = None,
    ) -> pd.DataFrame:
        """Merge df_seq (index ec_idx, matching events_seq) onto self events

        Merge sequence of tables corresponding to sequence of events onto
        events shared with self. If there are overlapping events, values from
        first tables will be used preferentially over those from later tables.

        Notes
        -----

        - Resulting table will have union of events and columns from the input
          tables.
        - Resulting table will be annotated with `self.ec_dataframe` (this
          ensures that any outdated information from those columns are updated
          to current events).

        Parameters
        ----------
        df_seq: Sequence[pd.DataFrame]
            Sequence of dataframes with index variable "ec_idx"
        events_seq: Sequence[Events]
            Sequence of events with same length as `df_seq`. Records in
            `df_seq[i]` refer to connections in `events_seq[i]` for each `i`.
            Each element must share genes with `self`.
        annotated_exons: Optional[Exons]
            If specified, use to reduce the number of exons called as denovo
        annotated_introns: Optional[GeneIntrons]
            If specified, use to reduce the number of introns called as denovo

        Returns
        -------
        pd.DataFrame
            Union of records in `df_seq` that belong to events shared between
            `self` and `events_seq`, preferring values found in later tables.
            Event information from `self.ec_dataframe` will be
            added/overwritten to table.
            An additional column `events_pass` will be added indicating which
            element of df_seq/events_seq each record came from.
        """
        if len(df_seq) != len(events_seq):
            raise ValueError(
                "(df_seq, events_seq) have mismatched lengths"
                f" ({len(df_seq) = }, {len(events_seq) = })"
            )
        result = pd.DataFrame()
        # iterate over sequences in reverse order
        for events_pass, (df, events) in reversed(
            list(enumerate(zip(df_seq, events_seq), 1))
        ):
            # get overlap between self and events
            overlap = self.unique_events_mask(events)
            # how to map from events to self, which pass
            events_to_self = pd.DataFrame(
                dict(
                    self_ec_idx=np.where(
                        self.broadcast_eidx_to_ecidx(~overlap.unique_events_mask)
                    )[0],
                ),
                index=pd.Index(
                    events.select_eidx_to_select_ecidx(overlap.shared_events_idx),
                    name="ec_idx",
                ),
            ).assign(events_pass=events_pass)
            # update result with table from this event
            result = result.combine_first(
                events_to_self.join(df, on="ec_idx", how="inner").set_index(
                    "self_ec_idx"
                )
            )
        # annotate with ec_dataframe
        ec_dataframe = self.ec_dataframe(
            annotated_exons=annotated_exons, annotated_introns=annotated_introns
        )
        result = ec_dataframe.join(
            result.drop(columns=ec_dataframe.columns, errors="ignore"),
            how="inner",
            on="ec_idx",
        )
        # add e_idx as column as well
        result["e_idx"] = self.connection_e_idx[result.index]
        # reorder columns
        REORDER_COLUMNS = ["e_idx", "events_pass"]
        result = result[
            [
                *REORDER_COLUMNS,
                *(x for x in result.columns if x not in REORDER_COLUMNS),
            ]
        ]
        return result
