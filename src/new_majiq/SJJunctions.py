"""
SJJunctions.py

contig regions not intersected by gene exons (stranded or unstranded), noting
if intersecting annotated junction or not (part of annotated exon)

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.Contigs import Contigs
from new_majiq.ContigRegions import ContigRegions
from new_majiq.internals import SJJunctions as _SJJunctions
from typing import (
    Optional,
    Union,
)
from pathlib import Path


class SJJunctions(ContigRegions):
    def __init__(self, sj_junctions: _SJJunctions):
        super().__init__(sj_junctions)
        return

    @property
    def _sj_junctions(self) -> _SJJunctions:
        """expose underlying internals representation of SJJunctions"""
        return self._contig_regions

    @property
    def sj_idx(self) -> np.ndarray:
        return self._region_idx

    @property
    def df(self) -> xr.Dataset:
        """view on underlying SJ junctions as xarray Dataset"""
        return xr.Dataset(
            {},
            {
                "sj_idx": self.sj_idx,
                "contig_idx": ("sj_idx", self.contig_idx),
                "start": ("sj_idx", self.start),
                "end": ("sj_idx", self.end),
                "strand": ("sj_idx", self.strand),
            },
        )

    def to_netcdf(self, path: Union[str, Path], mode: str) -> None:
        """Serialize to netcdf format. Note contigs need to be saved separately"""
        self.df.to_netcdf(path, mode, group=constants.NC_SJJUNCTIONS)
        return

    @classmethod
    def from_netcdf(
        cls,
        path: Union[str, Path],
        contigs: Optional[Contigs] = None,
    ) -> "SJJunctions":
        """Read SJJunctions from netcdf file

        Parameters
        ----------
        path: Union[str, Path]
            path to netcdf file
        contigs: Optional[Contigs]
            contigs on which the junctions are defined. If None, try loading
            from netcdf file.
        """
        df = xr.open_dataset(path, group=constants.NC_SJJUNCTIONS)
        if contigs is None:
            contigs = Contigs.from_netcdf(path)
        return SJJunctions(
            _SJJunctions(
                contigs._contigs,
                df.contig_idx,
                df.start,
                df.end,
                df.strand,
            )
        )
