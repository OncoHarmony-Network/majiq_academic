"""
EventsCoverage.py

Coverage over events from SJ bins

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Final, List, Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import new_majiq.constants as constants
from new_majiq.experiments import bam_experiment_name
from new_majiq.internals import EventsCoverage as _EventsCoverage

from ._workarounds import _load_zerodim_variables
from .Events import Events
from .ExonConnections import ExonConnections
from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions
from .SJExperiment import SJExperiment


class EventsCoverage(object):
    """coverage over events for an experiment"""

    def __init__(
        self,
        events_coverage: _EventsCoverage,
        bam_path: str,
        bam_version: str,
    ):
        self._events_coverage: Final[_EventsCoverage] = events_coverage
        self._bam_path: Final[str] = bam_path
        self._bam_version: Final[str] = bam_version
        return

    @property
    def bam_path(self) -> str:
        return self._bam_path

    @property
    def bam_version(self) -> str:
        return self._bam_version

    @property
    def events(self) -> Events:
        return Events(self._events_coverage._events)

    @property
    def numreads(self) -> npt.NDArray[np.float32]:
        """Scaled readrate for each event connection"""
        return self._events_coverage.numreads

    @property
    def numbins(self) -> npt.NDArray[np.float32]:
        """Number of nonzero bins for each event connection"""
        return self._events_coverage.numbins

    @property
    def bootstraps(self) -> npt.NDArray[np.float32]:
        """Bootstrapped read coverage replicates for each event connection"""
        return self._events_coverage.bootstraps

    @property
    def _df(self) -> xr.Dataset:
        return xr.Dataset(
            {
                "numreads": ("ec_idx", self.numreads),
                "numbins": ("ec_idx", self.numbins),
                "bootstraps": (("ec_idx", "bootstrap_replicate"), self.bootstraps),
            },
            {
                "ec_idx": self.events.ec_idx,
            },
            {
                "bam_path": self.bam_path,
                "bam_version": self.bam_version,
            },
        )

    @property
    def df(self) -> xr.Dataset:
        return xr.merge(
            (self._df, self.events.df), join="exact", combine_attrs="no_conflicts"
        )

    def to_zarr(self, path: Union[str, Path], consolidated: bool = True) -> None:
        """Serialize to zarr format"""
        if Path(path).exists():
            raise ValueError(
                f"Will not save EventsCoverage to existing file {path}."
                " Please delete and try again if desired or pick a different"
                " output path."
            )
        # save events, events coverage
        self.events.to_zarr(path, "w", consolidated=False)
        self._df.drop_vars("ec_idx").pipe(_load_zerodim_variables).to_zarr(
            path,
            mode="a",
            group=constants.NC_EVENTSCOVERAGE,
            consolidated=consolidated,
        )
        return

    def to_legacy_majiq(
        self, path: Union[str, Path], exon_connections: Optional[ExonConnections] = None
    ) -> None:
        """Save EventsCoverage to *.majiq file compatible with MAJIQ v2

        Parameters
        ----------
        path: Union[str, Path]
            Where to save output file
        exon_connections: Optional[ExonConnections]
            If specified, use exon connections to obtain full lsv descriptions
        """
        event_id: List[str]
        event_description: List[str]
        if exon_connections is None:
            gene_id = np.array(self.events.genes.gene_id)[
                self.events.exons.gene_idx[self.events.ref_exon_idx]
            ]
            ref_start = self.events.exons.start[self.events.ref_exon_idx]
            ref_end = self.events.exons.end[self.events.ref_exon_idx]
            event_id = [
                f"{g}:{t.decode('utf-8')}:"
                f"{s if s >= 0 else 'na'}-{e if e >= 0 else 'na'}"
                for g, t, s, e in zip(
                    gene_id, self.events.event_type, ref_start, ref_end
                )
            ]
            event_description = [
                f"{t.decode('utf-8')}|na{'|i' if has_intron else ''}"
                for t, has_intron in zip(
                    self.events.event_type,
                    self.events.is_intron[self.events.ec_idx_end - 1],
                )
            ]
        else:
            event_id = exon_connections.event_id(
                self.events.ref_exon_idx, self.events.event_type
            )
            event_description = exon_connections.event_description(
                self.events.ref_exon_idx, self.events.event_type
            )
        # arrays in output files
        lsv_types = np.array(
            [
                (eid.encode(), edesc.encode())
                for eid, edesc in zip(event_id, event_description)
            ],
            dtype=np.dtype("|S250, |S250"),
        )
        meta = np.array(
            [(bam_experiment_name(self.bam_path).encode(), self.bam_version.encode())],
            dtype=np.dtype("|S250, |S25"),
        )
        coverage = self.bootstraps
        junc_info = np.array(
            [
                (eid, start, end, reads, bins)
                for eid, start, end, reads, bins in zip(
                    self.events.broadcast_eidx_to_ecidx(event_id),
                    self.events.connection_start(),
                    self.events.connection_end(),
                    self.numreads,
                    self.numbins,
                )
            ],
            dtype=np.dtype("|S250, u4, u4, f4, f4"),
        )
        with open(path, "w+b") as handle:
            np.savez(
                handle,
                lsv_types=lsv_types,
                junc_info=junc_info,
                coverage=coverage,
                meta=meta,
            )
        return

    @classmethod
    def from_zarr(
        cls,
        path: Union[str, Path],
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> "EventsCoverage":
        events = Events.from_zarr(path, introns, junctions)
        with xr.open_zarr(path, group=constants.NC_EVENTSCOVERAGE) as df:
            return EventsCoverage(
                _EventsCoverage(
                    events._events,
                    df.numreads.values,
                    df.numbins.values,
                    df.bootstraps.values,
                ),
                df.bam_path,
                df.bam_version,
            )

    @classmethod
    def from_events_and_sj(
        cls,
        events: Events,
        sj: SJExperiment,
        num_bootstraps: int = constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        pvalue_threshold: float = constants.DEFAULT_COVERAGE_STACK_PVALUE,
    ) -> "EventsCoverage":
        return EventsCoverage(
            _EventsCoverage.from_sj(
                events._events,
                sj.junctions._sj_junctionsbins,
                sj.introns._sj_intronsbins,
                num_bootstraps,
                pvalue_threshold,
            ),
            sj.original_path,
            sj.original_version,
        )
