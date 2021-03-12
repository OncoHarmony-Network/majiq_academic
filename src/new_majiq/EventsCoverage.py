"""
EventsCoverage.py

Coverage over events from SJ bins

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.Events import Events
from new_majiq.SJJunctionsBins import SJJunctionsBins
from new_majiq.SJIntronsBins import SJIntronsBins
from new_majiq.internals import EventsCoverage as _EventsCoverage

from pathlib import Path
from typing import (
    Final,
    Union,
)


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
    def numreads(self) -> np.ndarray:
        """Scaled readrate for each event connection"""
        return self._events_coverage.numreads

    @property
    def numbins(self) -> np.ndarray:
        """Number of nonzero bins for each event connection"""
        return self._events_coverage.numbins

    @property
    def bootstraps(self) -> np.ndarray:
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

    def to_netcdf(self, path: Union[str, Path]) -> None:
        """Serialize to netcdf format"""
        if Path(path).exists():
            raise ValueError(
                f"Will not save EventsCoverage to existing file {path}."
                " Please delete and try again if desired or pick a different"
                " output path."
            )
        # save events, events coverage
        self.events.to_netcdf(path, "w")
        self._df.drop_vars("ec_idx").to_netcdf(
            path, "a", group=constants.NC_EVENTSCOVERAGE
        )
        return

    @classmethod
    def from_netcdf(
        cls,
        path: Union[str, Path],
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> "EventsCoverage":
        events = Events.from_netcdf(path, introns, junctions)
        df = xr.open_dataset(path, group=constants.NC_EVENTSCOVERAGE)
        return EventsCoverage(
            _EventsCoverage(
                events._events,
                df.numreads,
                df.numbins,
                df.bootstraps,
            ),
            df.bam_path,
            df.bam_version,
        )

    @classmethod
    def from_events_and_sj(
        cls,
        events: Events,
        sj_junctions: SJJunctionsBins,
        sj_introns: SJIntronsBins,
        num_bootstraps: int = constants.DEFAULT_BUILD_NUM_BOOTSTRAPS,
        pvalue_threshold: float = constants.DEFAULT_BUILD_STACK_PVALUE,
    ) -> "EventsCoverage":
        # only accept if bam_path/version are same in sj_junctions/introns
        if sj_junctions.original_path != sj_introns.original_path:
            raise ValueError(
                "sj_junctions and sj_introns do not share original bam path"
            )
        if sj_junctions.original_version != sj_introns.original_version:
            raise ValueError("sj_junctions and sj_introns not from same majiq version")
        return EventsCoverage(
            _EventsCoverage.from_sj(
                events._events,
                sj_junctions._sj_junctionsbins,
                sj_introns._sj_intronsbins,
                num_bootstraps,
                pvalue_threshold,
            ),
            sj_junctions.original_path,
            sj_junctions.original_version,
        )
