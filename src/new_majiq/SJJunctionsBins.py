"""
SJJunctionsBins.py

Bins with raw read coverage from input experiments

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.SJBinsReads import SJBinsReads
from new_majiq.SJJunctions import SJJunctions
from new_majiq.internals import SJJunctionsBins as _SJJunctionsBins
from new_majiq.internals import ExperimentStrandness
from new_majiq.version import version

from typing import (
    Final,
    Union,
)
from pathlib import Path


class SJJunctionsBins(SJBinsReads):
    def __init__(
        self,
        sj_junctionsbins: _SJJunctionsBins,
        original_path: str,
        original_version: str,
        original_time: str,
    ):
        super().__init__(sj_junctionsbins)
        self._original_path: Final[str] = original_path
        self._original_version: Final[str] = original_version
        self._original_time: Final[str] = original_time
        return

    @property
    def _sj_junctionsbins(self) -> _SJJunctionsBins:
        return self._sj_binsreads

    @property
    def original_path(self) -> str:
        return self._original_path

    @property
    def original_version(self) -> str:
        return self._original_version

    @property
    def original_time(self) -> str:
        return self._original_time

    @property
    def sjb_idx(self):
        return self._sjbin_idx

    @property
    def regions(self) -> SJJunctions:
        return SJJunctions(self._regions)

    @property
    def sjb_idx_start(self) -> np.ndarray:
        return self._region_idx_start

    @property
    def sjb_idx_end(self) -> np.ndarray:
        return self._region_idx_end

    @property
    def _df(self) -> xr.Dataset:
        """xr.Dataset view of SJJunctionsBins read counts"""
        return xr.Dataset(
            {
                "bin_reads": ("sjb_idx", self.bin_reads),
            },
            {
                "sjb_idx": self.sjb_idx,
                "bin_idx": ("sjb_idx", self.bin_idx),
                "_offsets": ("sj_offsets_idx", self._offsets),
            },
            {
                "total_bins": self.total_bins,
                "original_path": self.original_path,
                "original_version": self.original_version,
                "original_time": self.original_time,
            },
        )

    @property
    def df(self) -> xr.Dataset:
        """View of SJJunctionsBins (combined with underlying regions)"""
        return xr.merge(
            (
                self._df.drop_vars("_offsets"),
                self.regions.df.assign_coords(
                    sjb_idx_start=("sj_idx", self.sjb_idx_start),
                    sjb_idx_end=("sj_idx", self.sjb_idx_end),
                ),
            ),
            join="exact",
            combine_attrs="no_conflicts",
        )

    @classmethod
    def from_bam(
        cls,
        path: Union[str, Path],
        strandness: ExperimentStrandness = constants.DEFAULT_BAM_STRANDNESS,
        nthreads: int = constants.DEFAULT_BAM_NTHREADS,
    ) -> "SJJunctionsBins":
        """Load SJJunctionsBins from BAM file

        Parameters
        ----------
        path: Union[str, Path]
            Path to input BAM file
        strandness: ExperimentStrandness
            The strand-specificity of the experiment
        nthreads: int
            Number of threads to use to read in BAM file
        """
        # TODO: save/serialize strandness
        path = str(Path(path).resolve())
        original_version = version()
        original_time = str(np.datetime64("now"))
        return SJJunctionsBins(
            _SJJunctionsBins.from_bam(path, strandness, nthreads),
            path,
            original_version,
            original_time,
        )

    def to_zarr(self, path: Union[str, Path]) -> None:
        """Serialize to zarr format"""
        if Path(path).exists():
            raise ValueError(
                f"Will not save SJJunctionsBins to existing file {path}."
                " Please delete and try again if desired, or please pick a"
                " different output path."
            )
        self.regions.contigs.to_zarr(path, "w")
        self.regions.to_zarr(path, "a")
        self._df.drop_vars("sjb_idx").to_zarr(
            path, mode="a", group=constants.NC_SJJUNCTIONSBINS
        )
        return

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "SJJunctionsBins":
        """Load SJJunctionsBins from zarr format"""
        regions = SJJunctions.from_zarr(path)
        with xr.open_zarr(path, group=constants.NC_SJJUNCTIONSBINS) as df:
            return SJJunctionsBins(
                _SJJunctionsBins(
                    regions._sj_junctions,
                    df.bin_reads.values,
                    df.bin_idx.values,
                    df._offsets.values,
                    df.total_bins,
                ),
                df.original_path,
                df.original_version,
                df.original_time,
            )
