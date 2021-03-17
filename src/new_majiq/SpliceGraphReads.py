"""
SpliceGraphReads.py

Number of reads from SJ bins

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr
import dask.array as da

import new_majiq.constants as constants

from new_majiq.experiments import bam_experiment_name
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.SJIntronsBins import SJIntronsBins
from new_majiq.SJJunctionsBins import SJJunctionsBins
from new_majiq.internals import SpliceGraphReads as _SpliceGraphReads

from functools import cached_property
from pathlib import Path
from typing import (
    Final,
    List,
    Union,
)


class MultiSpliceGraphReads(object):
    """Open up multiple SpliceGraphReads objects all at once"""

    def __init__(self, sg_reads_paths: Union[str, List[Union[str, Path]]]):
        """ wraps xr.open_mfdataset to open up multiple SpliceGraphReads

        Parameters
        ----------
        sg_reads_paths: Union[str, List[Union[str, Path]]]
            string glob of paths to open or explicit list of files to open
        """
        self._df: Final[xr.Dataset] = xr.open_mfdataset(
            sg_reads_paths,
            concat_dim="experiment",
            compat="equals",
            engine="zarr",
            join="exact",
            group=constants.NC_SGREADS,
        )
        return

    @property
    def experiments(self) -> np.ndarray:
        """strings array of experiments"""
        return self._df.experiment.values

    @property
    def num_introns(self) -> int:
        return self._df.dims["gi_idx"]

    @property
    def num_junctions(self) -> int:
        return self._df.dims["gj_idx"]

    def select_experiment(
        self,
        experiment: str,
        gi_idx: Union[int, slice] = slice(None),
        gj_idx: Union[int, slice] = slice(None),
    ) -> xr.Dataset:
        """Load subset of specified experiment into memory"""
        return self._df.sel(experiment=experiment).isel(gi_idx=gi_idx, gj_idx=gj_idx).compute()

    def summarize_experiments(
        self,
        reduction: str = "median",
        gi_idx: Union[int, slice] = slice(None),
        gj_idx: Union[int, slice] = slice(None),
    ) -> xr.Dataset:
        """Perform xarray reduction over experiments for specified subset"""
        return getattr(self._df.isel(gi_idx=gi_idx, gj_idx=gj_idx), reduction)("experiment").compute()

    @cached_property
    def intron_checksum(self) -> int:
        """summarize intron hash: shared unique value or -1"""
        self._df.intron_hash.load()
        hashes = set(self._df.intron_hash.values)
        return hashes.pop() if len(hashes) == 1 else -1

    @cached_property
    def junction_checksum(self) -> int:
        """summarize junction hash: shared unique value or -1"""
        self._df.junction_hash.load()
        hashes = set(self._df.junction_hash.values)
        return hashes.pop() if len(hashes) == 1 else -1


class SpliceGraphReads(object):
    """reads over introns/junctions in splicegraph"""

    def __init__(
        self,
        sg_reads: _SpliceGraphReads,
        bam_path: str,
        bam_version: str,
    ):
        self._sg_reads: Final[_SpliceGraphReads] = sg_reads
        self._bam_path: Final[str] = bam_path
        self._bam_version: Final[str] = bam_version
        return

    @property
    def bam_path(self) -> str:
        return self._bam_path

    @property
    def experiment_name(self) -> str:
        return bam_experiment_name(self.bam_path)

    @property
    def bam_version(self) -> str:
        return self._bam_version

    @property
    def introns(self) -> GeneIntrons:
        return GeneIntrons(self._sg_reads._introns)

    @property
    def junctions(self) -> GeneJunctions:
        return GeneJunctions(self._sg_reads._junctions)

    @property
    def introns_reads(self) -> np.ndarray:
        return self._sg_reads.introns_reads

    @property
    def junctions_reads(self) -> np.ndarray:
        return self._sg_reads.junctions_reads

    @property
    def _df(self) -> xr.Dataset:
        return xr.Dataset(
            {
                "introns_reads": (
                    ("experiment", "gi_idx"),
                    self.introns_reads[np.newaxis],
                ),
                "junctions_reads": (
                    ("experiment", "gj_idx"),
                    self.junctions_reads[np.newaxis],
                ),
            },
            {
                "experiment": [self.experiment_name],
            },
            {
                "bam_path": self.bam_path,
                "bam_version": self.bam_version,
            },
        )

    @property
    def df_introns(self) -> xr.Dataset:
        return self.introns.df.assign(
            introns_reads=("gi_idx", self.introns_reads)
        ).assign_attrs(
            bam_path=self.bam_path,
            bam_version=self.bam_version,
        )

    @property
    def df_junctions(self) -> xr.Dataset:
        return self.junctions.df.assign(
            junctions_reads=("gj_idx", self.junctions_reads)
        ).assign_attrs(
            bam_path=self.bam_path,
            bam_version=self.bam_version,
        )

    def to_zarr(self, path: Union[str, Path], mode: str) -> None:
        """Save to specified zarr file"""
        (
            self._df
            .chunk(constants.NC_SGREADS_CHUNKS)  # type: ignore
            # add hash for introns/junctions
            .assign_coords(
                intron_hash=("experiment", [self.introns.checksum_nodata()]),
                junction_hash=("experiment", [self.junctions.checksum_nodata()]),
            ).to_zarr(path, mode=mode, group=constants.NC_SGREADS)
        )
        return

    @classmethod
    def from_zarr(
        cls,
        path: Union[str, Path],
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> "SpliceGraphReads":
        """Load from zarr file"""
        df = xr.open_zarr(path, group=constants.NC_SGREADS).squeeze("experiment")
        if df.intron_hash != introns.checksum_nodata():
            raise ValueError("Saved hash for introns does not match")
        if df.junction_hash != junctions.checksum_nodata():
            raise ValueError("Saved hash for junctions does not match")
        return SpliceGraphReads(
            _SpliceGraphReads(
                introns._gene_introns,
                junctions._gene_junctions,
                df.introns_reads.values,
                df.junctions_reads.values,
            ),
            df.bam_path,
            df.bam_version,
        )

    @classmethod
    def from_connections_and_sj(
        cls,
        introns: GeneIntrons,
        junctions: GeneJunctions,
        sj_introns: SJIntronsBins,
        sj_junctions: SJJunctionsBins,
    ) -> "SpliceGraphReads":
        # only accept if bam_path/version are same in sj_junctions/introns
        if sj_junctions.original_path != sj_introns.original_path:
            raise ValueError(
                "sj_junctions and sj_introns do not share original bam path"
            )
        if sj_junctions.original_version != sj_introns.original_version:
            raise ValueError("sj_junctions and sj_introns not from same majiq version")
        return SpliceGraphReads(
            _SpliceGraphReads.from_sj(
                introns._gene_introns,
                junctions._gene_junctions,
                sj_introns._sj_intronsbins,
                sj_junctions._sj_junctionsbins,
            ),
            sj_junctions.original_path,
            sj_junctions.original_version,
        )
