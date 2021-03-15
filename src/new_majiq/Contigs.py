"""
Contigs.py

Wrap new_majiq.internals.Contigs with Python for serialization/convenience
functions

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from typing import (
    Final,
    List,
    Optional,
    Union,
)
from new_majiq.internals import Contigs as _Contigs
from pathlib import Path


class Contigs(object):
    """Contigs/chromosomes used by a splicegraph"""

    def __init__(self, contigs: _Contigs):
        self._contigs: Final[_Contigs] = contigs
        return

    def checksum(self):
        return self._contigs.checksum()

    def __len__(self) -> int:
        """Number of contigs"""
        return len(self._contigs)

    @property
    def contig_idx(self) -> np.ndarray:
        return np.arange(len(self))

    @property
    def seqid(self) -> List[str]:
        return self._contigs.seqid

    def annotate_contig_idx(self, df: xr.Dataset) -> xr.Dataset:
        """For now, just add seqid to df using df.contig_idx"""
        return df.assign_coords(
            seqid=(df.contig_idx.dims, np.array(self.seqid)[df.contig_idx]),
        )

    @property
    def df(self) -> xr.Dataset:
        return xr.Dataset(
            {},
            {
                "contig_idx": self.contig_idx,
                "seqid": ("contig_idx", self.seqid),
            },
        )

    def __getitem__(self, seqid: str) -> int:
        """Get contig_idx for a given seqid"""
        try:
            return self._contigs[seqid]
        except IndexError:
            raise KeyError(f"{seqid = } not found in Contigs instance")

    def to_zarr(self, path: Union[str, Path], mode: str) -> None:
        self.df.drop_vars("contig_idx").to_zarr(
            path, mode=mode, group=constants.NC_CONTIGS
        )
        return

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "Contigs":
        with xr.open_zarr(path, group=constants.NC_CONTIGS) as df:
            return Contigs(_Contigs(df.seqid.values.tolist()))
