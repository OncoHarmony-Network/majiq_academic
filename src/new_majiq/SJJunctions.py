"""
SJJunctions.py

contig regions not intersected by gene exons (stranded or unstranded), noting
if intersecting annotated junction or not (part of annotated exon)

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Optional, Union

import numpy as np
import xarray as xr

import new_majiq.constants as constants
from new_majiq._workarounds import _load_zerodim_variables
from new_majiq.ContigRegions import ContigRegions
from new_majiq.Contigs import Contigs
from new_majiq.internals import SJJunctions as _SJJunctions


class SJJunctions(ContigRegions):
    def __init__(self, sj_junctions: _SJJunctions):
        super().__init__(sj_junctions)
        return

    @property
    def _sj_junctions(self) -> _SJJunctions:
        """expose underlying internals representation of SJJunctions"""
        return self._contig_regions

    def to_unstranded(self) -> "SJJunctions":
        """Get unique junctions ignoring strand (mark all unstranded)"""
        return SJJunctions(self._sj_junctions.to_unstranded())

    def flip_strand(self) -> "SJJunctions":
        """Get junctions with strands flipped in sorted order"""
        return SJJunctions(self._sj_junctions.flip_strand())

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

    def to_zarr(
        self, path: Union[str, Path], mode: str, consolidated: bool = True
    ) -> None:
        """Serialize to zarr format. Note contigs need to be saved separately"""
        self.df.drop_vars("sj_idx").pipe(lambda x: x.chunk(x.sizes)).pipe(
            _load_zerodim_variables
        ).to_zarr(
            path,
            mode=mode,
            group=constants.NC_SJJUNCTIONS,
            consolidated=consolidated,
        )
        return

    @classmethod
    def from_zarr(
        cls,
        path: Union[str, Path],
        contigs: Optional[Contigs] = None,
    ) -> "SJJunctions":
        """Read SJJunctions from zarr file

        Parameters
        ----------
        path: Union[str, Path]
            path to zarr file
        contigs: Optional[Contigs]
            contigs on which the junctions are defined. If None, try loading
            from zarr file.
        """
        if contigs is None:
            contigs = Contigs.from_zarr(path, group=constants.NC_SJJUNCTIONSCONTIGS)
        with xr.open_zarr(path, group=constants.NC_SJJUNCTIONS) as df:
            return SJJunctions(
                _SJJunctions(
                    contigs._contigs,
                    df.contig_idx.values,
                    df.start.values,
                    df.end.values,
                    df.strand.values,
                )
            )
