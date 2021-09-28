"""
SJIntronsBins.py

Bins with raw read coverage from input experiments

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import xarray as xr

import new_majiq.constants as constants
from new_majiq._version import version as nm_version
from new_majiq._workarounds import _load_zerodim_variables
from new_majiq.Exons import Exons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.internals import ExperimentStrandness
from new_majiq.internals import SJIntronsBins as _SJIntronsBins
from new_majiq.logger import get_logger
from new_majiq.SJBinsReads import SJBinsReads
from new_majiq.SJIntrons import SJIntrons


class SJIntronsBins(SJBinsReads):
    def __init__(
        self,
        sj_intronsbins: _SJIntronsBins,
        strandness: ExperimentStrandness,
        original_path: str,
        original_version: str,
        original_time: str,
    ):
        super().__init__(
            sj_intronsbins, strandness, original_path, original_version, original_time
        )
        return

    @property
    def _sj_intronsbins(self) -> _SJIntronsBins:
        return self._sj_binsreads

    @property
    def sib_idx(self):
        return self._sjbin_idx

    @property
    def regions(self) -> SJIntrons:
        return SJIntrons(self._regions)

    @property
    def sib_idx_start(self) -> np.ndarray:
        return self._region_idx_start

    @property
    def sib_idx_end(self) -> np.ndarray:
        return self._region_idx_end

    @property
    def _df(self) -> xr.Dataset:
        """xr.Dataset view of SJIntronsBins read counts"""
        return xr.Dataset(
            {
                "bin_reads": ("sib_idx", self.bin_reads),
            },
            {
                "sib_idx": self.sib_idx,
                "bin_idx": ("sib_idx", self.bin_idx),
                "_offsets": ("si_offsets_idx", self._offsets),
            },
            {
                "total_bins": self.total_bins,
                "strandness": self.strandness.name,
                "original_path": self.original_path,
                "original_version": self.original_version,
                "original_time": self.original_time,
            },
        )

    @property
    def df(self) -> xr.Dataset:
        """View of SJIntronsBins (combined with underlying regions)"""
        return xr.merge(
            (
                self._df.drop_vars("_offsets"),
                self.regions.df.assign_coords(
                    sib_idx_start=("si_idx", self.sib_idx_start),
                    sib_idx_end=("si_idx", self.sib_idx_end),
                ),
            ),
            join="exact",
            combine_attrs="no_conflicts",
        )

    @classmethod
    def from_bam(
        cls,
        path: Union[str, Path],
        total_bins: int,
        exons: Exons,
        gene_introns: GeneIntrons,
        strandness: ExperimentStrandness = constants.DEFAULT_BAM_STRANDNESS,
        nthreads: int = constants.DEFAULT_BAM_NTHREADS,
    ) -> "SJIntronsBins":
        """Load SJIntronsBins from BAM file

        Parameters
        ----------
        path: Union[str, Path]
            Path to input BAM file
        total_bins: int
            Maximum possible bin for intron coverage (usually obtained from
            SJJunctionsBins)
        exons: Exons
            exons used to define sj_introns on which to infer coverage
        gene_introns: GeneIntrons
            introns used to identify sj_introns that can only count for
            annotated introns
        strandness: ExperimentStrandness
            The strand-specificity of the experiment
        nthreads: int
            Number of threads to use to read in BAM file
        """
        # TODO: save information about splicegraph used
        path = str(Path(path).resolve())
        original_version = nm_version
        original_time = str(np.datetime64("now"))
        return SJIntronsBins(
            _SJIntronsBins.from_bam(
                path,
                total_bins,
                exons._exons,
                gene_introns._gene_introns,
                strandness,
                nthreads,
            ),
            strandness,
            path,
            original_version,
            original_time,
        )

    def to_zarr(self, path: Union[str, Path], consolidated: bool = True) -> None:
        """Serialize to zarr format

        Will only write to existing file if has SJJunctionsBins with same
        original path and version
        """
        if Path(path).exists():
            try:
                with xr.open_zarr(path, group=constants.NC_SJJUNCTIONSBINS) as x:
                    if x.original_path != self.original_path:
                        raise ValueError(
                            f"Will not save SJIntronsBins to existing file {path}"
                            " which appears to have junctions reads from"
                            f" path {x.original_path}"
                            f" vs introns from path {self.original_path}"
                        )
                    if x.original_version != self.original_version:
                        raise ValueError(
                            f"Will not save SJIntronsBins to existing file {path}"
                            " which appears to have junctions reads from"
                            f" version {x.original_version}"
                            f" vs introns from version {self.original_version}"
                        )
            except (OSError, AttributeError):
                raise ValueError(
                    f"Will not save SJIntronsBins to existing file {path}"
                    " since it does not appear to have matching junctions reads"
                )
            # otherwise, appears to have compatible SJJunctionsBins which
            # implies that they have matching contigs
        # otherwise
        self.regions.contigs.to_zarr(path, "w", consolidated=False)
        self.regions.to_zarr(path, "a", consolidated=False)
        self._df.drop_vars("sib_idx").pipe(lambda x: x.chunk(x.sizes)).pipe(
            _load_zerodim_variables
        ).to_zarr(
            path,
            mode="a",
            group=constants.NC_SJINTRONSBINS,
            consolidated=consolidated,
        )
        return

    @classmethod
    def from_regions(
        cls,
        introns: SJIntrons,
        total_bins: int,
        strandness: ExperimentStrandness = ExperimentStrandness.NONE,
        original_path: str = "<none>",
        original_version: str = nm_version,
        original_time: Optional[str] = None,
    ) -> "SJIntronsBins":
        """Empty SJIntronsBins matched to input introns"""
        if original_time is None:
            original_time = str(np.datetime64("now"))
        return SJIntronsBins(
            _SJIntronsBins(
                introns._sj_introns,
                [],
                [],
                np.zeros(1 + len(introns), dtype=np.uint64),
                total_bins,
            ),
            strandness,
            original_path,
            original_version,
            original_time,
        )

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "SJIntronsBins":
        """Load SJIntronsBins from zarr format"""
        regions = SJIntrons.from_zarr(path)
        with xr.open_zarr(path, group=constants.NC_SJINTRONSBINS) as df:
            try:
                strandness = ExperimentStrandness(ord(df.strandness[0]))
            except AttributeError:
                get_logger().warning(
                    f"SJJunctionsBins in {path} did not save strandness"
                    " -> defaulting to NONE"
                )
                strandness = ExperimentStrandness.NONE
            return SJIntronsBins(
                _SJIntronsBins(
                    regions._sj_introns,
                    df.bin_reads.values,
                    df.bin_idx.values,
                    df._offsets.values,
                    df.total_bins,
                ),
                strandness,
                df.original_path,
                df.original_version,
                df.original_time,
            )
