"""
SJIntronsBins.py

Bins with raw read coverage from input experiments

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Final, Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import new_majiq.constants as constants
from new_majiq._version import version as nm_version
from new_majiq.internals import ExperimentStrandness
from new_majiq.internals import SJIntronsBins as _SJIntronsBins
from new_majiq.logger import get_logger

from ._workarounds import _load_zerodim_variables
from .Exons import Exons
from .GeneIntrons import GeneIntrons
from .SJBinsReads import SJBinsReads
from .SJIntrons import SJIntrons


class SJIntronsBins(SJBinsReads):
    """Per-bin read coverage over introns

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
    gene_introns_checksum, exons_checksum: Optional[int]
        Checksum of gene introns, exons used to define regions for coverage

    Notes
    -----
    Intron coverage is dependent on the gene introns and exons used to define
    regions on which to assess coverage.

    See Also
    --------
    SJExperiment.from_bam
    SJIntronsBins.from_bam
    SJIntronsBins.from_zarr
    """

    def __init__(
        self,
        sj_intronsbins: _SJIntronsBins,
        strandness: ExperimentStrandness,
        original_path: str,
        original_version: str,
        original_time: str,
        gene_introns_checksum: Optional[int] = None,
        exons_checksum: Optional[int] = None,
    ):
        super().__init__(
            sj_intronsbins, strandness, original_path, original_version, original_time
        )
        # if provided, checksums indicating which introns/exons used in construction
        self.gene_introns_checksum: Final[Optional[int]] = gene_introns_checksum
        self.exons_checksum: Final[Optional[int]] = exons_checksum
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
    def sib_idx_start(self) -> npt.NDArray[np.uint64]:
        return self._region_idx_start

    @property
    def sib_idx_end(self) -> npt.NDArray[np.uint64]:
        return self._region_idx_end

    @property
    def _df(self) -> xr.Dataset:
        """xr.Dataset view of SJIntronsBins read counts"""
        checksums = {}
        if self.gene_introns_checksum is not None:
            checksums["gene_introns_checksum"] = self.gene_introns_checksum
        if self.exons_checksum is not None:
            checksums["exons_checksum"] = self.exons_checksum
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
                **checksums,  # type: ignore[arg-type]
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
            gene_introns_checksum=gene_introns.checksum_nodata(),
            exons_checksum=exons.checksum(),
        )

    def to_zarr(
        self,
        path: Union[str, Path],
        consolidated: bool = True,
        check_experiment_if_exists: bool = True,
    ) -> None:
        """Serialize to zarr format

        Will only write to existing file if has SJJunctionsBins with same
        original path and version

        Parameters
        ----------
        path: Union[str, Path]
            Path for output zarr file
        consolidated: bool
            Should the zarr file be consolidated (e.g. we don't expect array
            metadata to change) at the end of function?
        check_experiment_if_exists: bool
            If zarr file already exists at path, check that it has
            SJJunctionsBins and that the original path and version match
        """
        if check_experiment_if_exists and Path(path).exists():
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
    def from_arrays(
        cls,
        regions: SJIntrons,
        bin_reads: npt._ArrayLikeFloat_co,
        bin_idx: npt._ArrayLikeInt_co,
        offsets: npt._ArrayLikeInt_co,
        total_bins: int,
        strandness: ExperimentStrandness,
        original_path: str = "<none>",
        original_version: str = nm_version,
        original_time: Optional[str] = None,
        gene_introns_checksum: Optional[int] = None,
        exons_checksum: Optional[int] = None,
    ) -> "SJIntronsBins":
        """Create :class:`SJIntronsBins` from :class:`SJIntrons`, arrays, metadata"""
        if original_time is None:
            original_time = str(np.datetime64("now"))
        return SJIntronsBins(
            _SJIntronsBins(
                regions._sj_introns, bin_reads, bin_idx, offsets, total_bins
            ),
            strandness,
            original_path,
            original_version,
            original_time,
            gene_introns_checksum=gene_introns_checksum,
            exons_checksum=exons_checksum,
        )

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
        return SJIntronsBins.from_arrays(
            introns,
            np.array([], dtype=np.float32),
            np.array([], dtype=np.int32),
            np.zeros(1 + len(introns), dtype=np.uint64),
            total_bins,
            strandness,
            original_path=original_path,
            original_version=original_version,
            original_time=original_time,
        )

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "SJIntronsBins":
        """Load SJIntronsBins from zarr format"""
        regions = SJIntrons.from_zarr(path)
        with xr.open_zarr(path, group=constants.NC_SJINTRONSBINS) as df:
            df.load()
            try:
                strandness = ExperimentStrandness(ord(df.strandness[0]))
            except AttributeError:
                get_logger().warning(
                    f"SJJunctionsBins in {path} did not save strandness"
                    " -> defaulting to NONE"
                )
                strandness = ExperimentStrandness.NONE
            return SJIntronsBins.from_arrays(
                regions,
                df.bin_reads.values,
                df.bin_idx.values,
                df._offsets.values,
                df.total_bins,
                strandness,
                original_path=df.original_path,
                original_version=df.original_version,
                original_time=df.original_time,
                gene_introns_checksum=df.attrs.get("gene_introns_checksum"),
                exons_checksum=df.attrs.get("exons_checksum"),
            )
