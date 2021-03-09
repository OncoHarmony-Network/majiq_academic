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


class Contigs(_Contigs):
    """ Contigs/chromosomes used by a splicegraph
    """

    def __init__(self, contigs: _Contigs):
        self._contigs: Final[_Contigs] = contigs
        return

    @classmethod
    def FromSeqids(self, seqids: List[str]) -> "Contigs":
        return Contigs(_Contigs(seqids))

    @property
    def seqid(self) -> List[str]:
        return self._contigs.seqid

    def __len__(self) -> int:
        """ Number of contigs
        """
        return len(self._contigs)

    @property
    def contig_idx(self):
        return np.arange(len(self))

    @property
    def df(self) -> xr.Dataset:
        return xr.Dataset(
            {},
            {
                "contig_idx": self.contig_idx,
                "seqid": ("contig_idx", self.seqid),
            }
        )

    def __getitem__(self, seqid: str) -> Optional[int]:
        """ Get contig_idx for a given seqid. None if no such seqid.
        """
        raise NotImplementedError("Need to expose safe_idx in internals")

    def to_netcdf(self, path: Union[str, Path], mode: str) -> None:
        self.df.to_netcdf(path, mode, constants.NC_CONTIGS)
        return

    @classmethod
    def from_netcdf(cls, path: Union[str, Path]) -> "Contigs":
        df = xr.open_dataset(path, group=constants.NC_CONTIGS)
        return Contigs(_Contigs(df.seqid.values.tolist()))
