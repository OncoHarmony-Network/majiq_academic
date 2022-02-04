"""
SJJunctionsBins.py

Bins with raw read coverage from input experiments

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import new_majiq.constants as constants
from new_majiq._version import version as nm_version
from new_majiq.internals import ExperimentStrandness
from new_majiq.internals import SJJunctionsBins as _SJJunctionsBins
from new_majiq.logger import get_logger

from .SJBinsReads import SJBinsReads
from .SJJunctions import SJJunctions


class SJJunctionsBins(SJBinsReads):
    """Per-bin read coverage over junctions

    Parameters
    ----------
    sj_intronsbins: _SJIntronsBins
        Underlying object binding the internal C++ API
    strandness: ExperimentStrandness
        strandness that was used for the experiment
    original_path: str
        Original path to BAM file from which coverage was derived
    original_version, original_time: str
        Version of MAJIQ and time BAM file was processed

    See Also
    --------
    SJExperiment.from_bam
    SJJunctionsBins.from_bam
    SJJunctionsBins.from_zarr
    """

    def __init__(
        self,
        sj_junctionsbins: _SJJunctionsBins,
        strandness: ExperimentStrandness,
        original_path: str,
        original_version: str,
        original_time: str,
    ):
        super().__init__(
            sj_junctionsbins, strandness, original_path, original_version, original_time
        )
        return

    def __repr__(self) -> str:
        """String representation of SJJunctionsBins"""
        return (
            f"SJJunctionsBins(strandness={self.strandness.name},"
            f" original_path=.../{Path(self.original_path).name})"
        )

    @property
    def _sj_junctionsbins(self) -> _SJJunctionsBins:
        return self._sj_binsreads

    def project_reads(
        self,
        new_junctions: SJJunctions,
        new_strandness: ExperimentStrandness,
        flip_strand: bool = False,
    ) -> "SJJunctionsBins":
        """Project reads from self onto matching junctions in new_junctions

        Parameters
        ----------
        new_junctions: SJJunctions
            Junctions bins should be defined against
        new_strandness: ExperimentStrandness
            What strandness the resulting SJJunctionsBins should say they are. Needed since we might be converting from stranded to unstranded
        flip_strand: bool
            map junctions from self to opposite strand rather than matching strands
        """
        return SJJunctionsBins(
            self._sj_junctionsbins.project_reads(
                new_junctions._sj_junctions, flip_strand
            ),
            new_strandness,
            self.original_path,
            self.original_version,
            self.original_time,
        )

    def to_unstranded(self) -> "SJJunctionsBins":
        """Convert stranded junction reads to unstranded junctions"""
        return self.project_reads(
            self.regions.to_unstranded(), ExperimentStrandness.NONE
        )

    def flip_strand(self) -> "SJJunctionsBins":
        """Flip +/- strand reads to -/+ strand reads"""
        if self.strandness == ExperimentStrandness.FORWARD:
            new_strandness = ExperimentStrandness.REVERSE
        elif self.strandness == ExperimentStrandness.REVERSE:
            new_strandness = ExperimentStrandness.FORWARD
        else:
            new_strandness = ExperimentStrandness.NONE
        return self.project_reads(
            self.regions.flip_strand(), new_strandness, flip_strand=True
        )

    @property
    def sjb_idx(self):
        return self._sjbin_idx

    @property
    def regions(self) -> SJJunctions:
        return SJJunctions(self._regions)

    @property
    def sjb_idx_start(self) -> npt.NDArray[np.uint64]:
        return self._region_idx_start

    @property
    def sjb_idx_end(self) -> npt.NDArray[np.uint64]:
        return self._region_idx_end

    @property
    def _df(self) -> xr.Dataset:
        """xr.Dataset view of SJJunctionsBins read counts"""
        return xr.Dataset(
            {
                "bin_reads": ("sjb_idx", self.bin_reads),
            },
            {
                "bin_idx": ("sjb_idx", self.bin_idx),
                "_offsets": ("sj_offsets_idx", self._offsets),
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
        """View of SJJunctionsBins (combined with underlying regions)"""
        return xr.merge(
            (
                # arrays defining sparse representation
                self._df.drop_vars(["bin_reads", "bin_idx", "_offsets"]),
                # sparse representation
                xr.Dataset(
                    {"bin_reads": (("sj_idx", "bin_idx"), self.bin_reads_2d)},
                    {"bin_idx": np.arange(self.total_bins)},
                ),
                # information about junctions
                self.regions.df,
            ),
            join="override",
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
        path = str(Path(path).resolve())
        original_version = nm_version
        original_time = str(np.datetime64("now"))
        return SJJunctionsBins(
            _SJJunctionsBins.from_bam(path, strandness, nthreads),
            strandness,
            path,
            original_version,
            original_time,
        )

    def to_zarr(self, path: Union[str, Path], consolidated: bool = True) -> None:
        """Serialize to zarr format"""
        if Path(path).exists():
            raise ValueError(
                f"Will not save SJJunctionsBins to existing file {path}."
                " Please delete and try again if desired, or please pick a"
                " different output path."
            )
        self.regions.contigs.to_zarr(
            path, "w", group=constants.NC_SJJUNCTIONSCONTIGS, consolidated=False
        )
        self.regions.to_zarr(path, "a", consolidated=False)
        self._df.pipe(lambda x: x.chunk(x.sizes)).to_zarr(
            path,
            mode="a",
            group=constants.NC_SJJUNCTIONSBINS,
            consolidated=consolidated,
        )
        return

    @classmethod
    def original_path_from_zarr(cls, path: Union[str, Path]) -> str:
        """Extract bam path from zarr with SJJunctionsBins"""
        with xr.open_zarr(path, group=constants.NC_SJJUNCTIONSBINS) as df:
            return df.original_path

    @classmethod
    def from_arrays(
        cls,
        regions: SJJunctions,
        bin_reads: npt._ArrayLikeInt_co,
        bin_idx: npt._ArrayLikeInt_co,
        offsets: npt._ArrayLikeInt_co,
        total_bins: int,
        strandness: ExperimentStrandness,
        original_path: str = "<none>",
        original_version: str = nm_version,
        original_time: Optional[str] = None,
    ) -> "SJJunctionsBins":
        """Create :class:`SJJunctionsBins` from :class:`SJJunctions`, arrays, metadata"""
        if original_time is None:
            original_time = str(np.datetime64("now"))
        return SJJunctionsBins(
            _SJJunctionsBins(
                regions._sj_junctions, bin_reads, bin_idx, offsets, total_bins
            ),
            strandness,
            original_path,
            original_version,
            original_time,
        )

    @classmethod
    def with_summaries(
        cls,
        regions: SJJunctions,
        numreads: npt._ArrayLikeInt_co,
        numbins: npt._ArrayLikeInt_co,
        total_bins: int,
        strandness: Optional[ExperimentStrandness] = None,
    ) -> "SJJunctionsBins":
        """Create :class:`SJJunctionsBins` with uniform coverage over nonzero bins"""
        # cast to integer if necessary
        numreads = np.array(numreads, copy=False)
        numbins = np.array(numbins, copy=False)
        # check valid sizes
        if numreads.ndim != 1 or numbins.ndim != 1:
            raise ValueError("numreads, numbins must be 1D")
        if len(regions) != numreads.shape[0] or len(regions) != numbins.shape[0]:
            raise ValueError("numreads, numbins must match regions in length")
        # check valid values
        if np.any(numreads < numbins):
            raise ValueError("numreads must be at least numbins")
        if np.any(numbins < 0) or np.any(numbins > total_bins):
            raise ValueError(f"numbins must be in [0, {total_bins = }]")
        # offsets are determined by numbins
        offsets = np.empty(1 + len(regions), dtype=np.int64)
        offsets[0] = 0
        np.cumsum(numbins, out=offsets[1:])
        # how many reads per bin (lower value, deal with remainder after)
        reads_per_bin = numreads // numbins
        # how many bins will get an extra read?
        remainder_bins = numreads % numbins
        # construct bin_reads so that remainders go on end each slice per region
        nonzero_bin_idx = np.repeat(offsets[1:], numbins) + np.arange(
            -1, -1 - offsets[-1], -1
        )
        bin_reads = np.repeat(reads_per_bin, numbins) + np.where(
            nonzero_bin_idx < np.repeat(remainder_bins, numbins), 1, 0
        )
        # make bin_idx arbitrarily.
        # To make sure there are no repeats, we tile the possible indexes.
        tile_ct = -(len(bin_reads) // -total_bins)  # ceil(bin_reads / total_bins)
        bin_idx = np.tile(np.arange(total_bins), tile_ct)[: len(bin_reads)]
        # guess strandness if not specified (REVERSE vs NONE)
        if strandness is None:
            strandness = (
                ExperimentStrandness.REVERSE
                if np.any(regions.strand != b".")
                else ExperimentStrandness.NONE
            )
        # use these arrays to construct simple SJJunctionsBins
        return SJJunctionsBins.from_arrays(
            regions,
            bin_reads.astype(np.int32, copy=False),
            bin_idx.astype(np.int32, copy=False),
            offsets.astype(np.uint64, copy=False),
            total_bins,
            strandness,
        )

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "SJJunctionsBins":
        """Load SJJunctionsBins from zarr format"""
        regions = SJJunctions.from_zarr(path)
        with xr.open_zarr(path, group=constants.NC_SJJUNCTIONSBINS) as df:
            df.load()
            try:
                strandness = ExperimentStrandness(ord(df.strandness[0]))
            except AttributeError:
                get_logger().warning(
                    f"SJJunctionsBins in {path} did not save strandness"
                    " -> defaulting to NONE"
                )
                strandness = ExperimentStrandness.NONE
            return SJJunctionsBins.from_arrays(
                regions,
                df.bin_reads.values,
                df.bin_idx.values,
                df._offsets.values,
                df.total_bins,
                strandness,
                original_path=df.original_path,
                original_version=df.original_version,
                original_time=df.original_time,
            )
