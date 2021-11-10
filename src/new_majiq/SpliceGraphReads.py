"""
SpliceGraphReads.py

Number of reads from SJ bins

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Final,
    Hashable,
    List,
    Optional,
    Union,
    cast,
)

import dask.array as da
import numpy as np
import xarray as xr
from dask.delayed import Delayed
from dask.distributed import progress

import new_majiq.constants as constants
from new_majiq._workarounds import _load_zerodim_variables
from new_majiq.experiments import bam_experiment_name
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.internals import SpliceGraphReads as _SpliceGraphReads

if TYPE_CHECKING:
    from new_majiq.SJExperiment import SJExperiment


class SpliceGraphReads(object):
    """Coverage over introns and junctions for arbitrary number of samples

    Parameters
    ----------
    df: xr.Dataset
        Data variables:
            introns_reads[prefix, gi_idx]
            junctions_reads[prefix, gj_idx]
            intron_hash[prefix]
            junction_hash[prefix]
        Coordinates:
            prefix[prefix]
    """

    def __init__(self, df: xr.Dataset):
        for required in (
            "introns_reads",
            "junctions_reads",
            "intron_hash",
            "junction_hash",
        ):
            if required not in df.data_vars.keys():
                raise ValueError(
                    f"SpliceGraphReads.df does not have required variable {required}"
                )
        if "prefix" not in df.indexes.keys():
            raise ValueError(
                "SpliceGraphReads.df does not have required index for prefix dimension"
            )
        if set(df.introns_reads.dims) != {"prefix", "gi_idx"}:
            raise ValueError(
                "SpliceGraphReads.df.introns_reads must have dimensions (prefix, gi_idx)"
            )
        if set(df.junctions_reads.dims) != {"prefix", "gj_idx"}:
            raise ValueError(
                "SpliceGraphReads.df.junctions_reads must have dimensions (prefix, gj_idx)"
            )
        self.df: Final[xr.Dataset] = df
        return

    def __getitem__(self, prefixes) -> "SpliceGraphReads":
        """Get subset of SpliceGraphReads corresponding to selected prefixes"""
        if isinstance(prefixes, str):
            prefixes = [prefixes]
        return SpliceGraphReads(self.df.sel(prefix=prefixes))

    def summarize(
        self, new_prefix: str, reduction: str = "median"
    ) -> "SpliceGraphReads":
        """Perform xarray reduction over experiments"""
        df = (
            getattr(self.df.drop_vars(["intron_hash", "junction_hash"]), reduction)(
                "prefix"
            )
            .assign(
                intron_hash=self.intron_checksum,
                junction_hash=self.junction_checksum,
            )
            .expand_dims(prefix=[new_prefix])
            .assign_attrs(summary_function=reduction, prefixes_summarized=self.prefixes)
        )
        return SpliceGraphReads(df)

    @property
    def num_introns(self) -> int:
        return self.df.sizes["gi_idx"]

    @property
    def num_junctions(self) -> int:
        return self.df.sizes["gj_idx"]

    @property
    def num_prefixes(self) -> int:
        return self.df.sizes["prefix"]

    @property
    def prefixes(self) -> List[str]:
        return self.df["prefix"].values.tolist()

    @cached_property
    def intron_checksum(self) -> int:
        """summarize intron hash: shared unique value or -1"""
        hashes = set(self.df.intron_hash.load().values)
        return hashes.pop() if len(hashes) == 1 else -1

    @cached_property
    def junction_checksum(self) -> int:
        """summarize junction hash: shared unique value or -1"""
        hashes = set(self.df.junction_hash.load().values)
        return hashes.pop() if len(hashes) == 1 else -1

    @property
    def introns_reads(self) -> xr.DataArray:
        return self.df.introns_reads

    @property
    def junctions_reads(self) -> xr.DataArray:
        return self.df.junctions_reads

    def to_zarr(
        self,
        path: Union[str, Path],
        mode: str = "w",
        chunksize: int = constants.NC_SGREADS_CHUNKS,
        show_progress: bool = False,
    ) -> None:
        """Save to specified zarr file"""
        save_df_future = cast(
            Delayed,
            self.df.chunk(chunksize)
            .pipe(_load_zerodim_variables)
            .to_zarr(
                path,
                mode=mode,
                group=constants.NC_SGREADS,
                consolidated=True,
                compute=False,
            ),
        )
        if show_progress:
            save_df_future = save_df_future.persist()
            progress(save_df_future)
        else:
            save_df_future.compute()
        return

    def to_zarr_slice(
        self,
        path: Union[str, Path],
        prefix_slice: slice,
        chunksize: int = constants.NC_SGREADS_CHUNKS,
    ) -> None:
        """Save SpliceGraphReads to specified slice of existing zarr file

        Save SpliceGraphReads to specified slice of existing zarr file.
        Typically run after SpliceGraphReads.to_zarr_slice_init()

        Parameters
        ----------
        path: Union[str, Path]
            Path for existing output Zarr for SpliceGraphReads output
        prefix_slice:
            slice (of prefix dimension) of existing store to save to
        chunksize: int
            Chunksize over introns, junctions dimensions
        """
        self.df.drop_vars("prefix").to_zarr(
            path, group=constants.NC_SGREADS, region=dict(prefix=prefix_slice)
        )
        return

    @classmethod
    def to_zarr_slice_init(
        cls,
        path: Union[str, Path],
        prefixes: List[str],
        num_introns: int,
        num_junctions: int,
        chunksize: int = constants.NC_SGREADS_CHUNKS,
        reads_dtype: type = np.float32,
        attrs: Dict[Hashable, Any] = dict(),
    ) -> None:
        """Init zarr for SpliceGraphReads over many prefixes for multithreaded write

        Initialize zarr for SpliceGraphReads over many prefixes for
        multithreaded write. For use with to_zarr_slice()

        Parameters
        ----------
        path: Union[str, Path]
            Path for output Zarr for SpliceGraphReads output
        prefixes: List[str]
            Values for the prefix dimension coordinate
        num_introns, num_junctions: int
            Sizes for the intron, junction dimensions
        chunksize: int
            Chunksize over introns, junctions dimensions
        reads_dtype: type
            What type to use for junctions_reads, introns_reads
        attrs: Dict[Hashable, Any]
            attributes to store with the zarr file
        """
        hashes_arr = da.empty(len(prefixes), dtype=int, chunks=1)
        introns_arr = da.empty(
            (len(prefixes), num_introns), dtype=reads_dtype, chunks=(1, chunksize)
        )
        junctions_arr = da.empty(
            (len(prefixes), num_junctions), dtype=reads_dtype, chunks=(1, chunksize)
        )
        # save basic metadata
        xr.Dataset(
            dict(
                introns_reads=(("prefix", "gi_idx"), introns_arr),
                junctions_reads=(("prefix", "gj_idx"), junctions_arr),
                junction_hash=("prefix", hashes_arr),
                intron_hash=("prefix", hashes_arr),
            ),
        ).to_zarr(
            path,
            mode="w",
            compute=False,
            group=constants.NC_SGREADS,
            consolidated=False,
        )
        # save prefixes, attributes
        xr.Dataset({}, {"prefix": prefixes}, attrs).to_zarr(
            path, mode="a", group=constants.NC_SGREADS, consolidated=True
        )
        return

    @classmethod
    def from_zarr(
        cls, path: Union[str, Path, List[Union[str, Path]]]
    ) -> "SpliceGraphReads":
        """Load one or more SpliceGraphReads files together at once

        Load one or more SpliceGraphReads files together at once. If they have
        overlapping prefixes, data will be loaded from the first file with the
        given prefix.
        """
        if not isinstance(path, list):
            path = [path]
        df = xr.open_mfdataset(
            path,
            engine="zarr",
            group=constants.NC_SGREADS,
            combine="nested",
            concat_dim="prefix",
            join="override",
            compat="override",
            coords="minimal",
            data_vars="minimal",
        )
        if len(path) > 1:
            # attributes are defined by path[0]. We'd rather just have none
            df.attrs.clear()
        return SpliceGraphReads(df)

    def _to_internals(
        self,
        introns: GeneIntrons,
        junctions: GeneJunctions,
        prefix: Optional[str] = None,
    ) -> _SpliceGraphReads:
        df: xr.Dataset
        if prefix is None:
            if self.num_prefixes > 1:
                raise ValueError(
                    "no prefix was specified for SpliceGraphReads with > 1 prefix"
                )
            df = self.df.squeeze("prefix")
        else:
            df = self.df.sel(prefix=prefix)
        if df.intron_hash.values[()] != introns.checksum_nodata():
            raise ValueError("Intron hash does not match")
        if df.junction_hash.values[()] != junctions.checksum_nodata():
            raise ValueError("Junction hash does not match")
        return _SpliceGraphReads(
            introns._gene_introns,
            junctions._gene_junctions,
            df.introns_reads.values,
            df.junctions_reads.values,
        )

    @classmethod
    def _from_internals(
        cls, sgreads: _SpliceGraphReads, bam_path: str, bam_version: str
    ) -> "SpliceGraphReads":
        df = xr.Dataset(
            {
                "junctions_reads": (
                    ("prefix", "gj_idx"),
                    np.array([sgreads.junctions_reads]),
                ),
                "introns_reads": (
                    ("prefix", "gi_idx"),
                    np.array([sgreads.introns_reads]),
                ),
                "junction_hash": ("prefix", [sgreads._junctions.checksum_nodata()]),
                "intron_hash": ("prefix", [sgreads._introns.checksum_nodata()]),
            },
            {
                "prefix": [bam_experiment_name(bam_path)],
            },
            {
                "bam_path": bam_path,
                "bam_version": bam_version,
            },
        )
        return SpliceGraphReads(df)

    @staticmethod
    def _internals_from_connections_and_sj(
        introns: GeneIntrons,
        junctions: GeneJunctions,
        sj: "SJExperiment",
    ) -> _SpliceGraphReads:
        return _SpliceGraphReads.from_sj(
            introns._gene_introns,
            junctions._gene_junctions,
            sj.introns._sj_intronsbins,
            sj.junctions._sj_junctionsbins,
        )

    @classmethod
    def from_connections_and_sj(
        cls,
        introns: GeneIntrons,
        junctions: GeneJunctions,
        sj: "SJExperiment",
    ) -> "SpliceGraphReads":
        sgreads: _SpliceGraphReads = cls._internals_from_connections_and_sj(
            introns, junctions, sj
        )
        return cls._from_internals(sgreads, sj.original_path, sj.original_version)
