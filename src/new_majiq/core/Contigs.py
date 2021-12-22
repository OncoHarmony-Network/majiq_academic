"""
Contigs.py

Wrap new_majiq.internals.Contigs with Python for serialization/convenience
functions

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Final, List, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import new_majiq.constants as constants
from new_majiq.internals import Contigs as _Contigs

from ._workarounds import _load_zerodim_variables


class Contigs(object):
    """Collection of contigs/chromosomes on which genes can be defined

    Parameters
    ----------
    contigs: _Contigs
        Underlying object binding the internal C++ API

    See Also
    --------
    Contigs.from_list
    Contigs.from_zarr
    """

    def __init__(self, contigs: _Contigs):
        self._contigs: Final[_Contigs] = contigs
        return

    def checksum(self):
        return self._contigs.checksum()

    def __len__(self) -> int:
        """Number of contigs"""
        return len(self._contigs)

    @property
    def contig_idx(self) -> npt.NDArray[np.int64]:
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

    def to_zarr(
        self,
        path: Union[str, Path],
        mode: str,
        group: str = constants.NC_CONTIGS,
        consolidated: bool = True,
    ) -> None:
        self.df.drop_vars("contig_idx").pipe(lambda x: x.chunk(x.sizes)).pipe(
            _load_zerodim_variables
        ).to_zarr(
            path,
            mode=mode,
            group=group,
            consolidated=consolidated,
        )
        return

    @classmethod
    def from_list(cls, seqid: List[str]) -> "Contigs":
        """Create :class:`Contigs` for the list of seqids

        Parameters
        ----------
        seqid: List[str]
            Unique names of independent contigs
        """
        return Contigs(_Contigs(seqid))

    @classmethod
    def from_zarr(
        cls, path: Union[str, Path], group: str = constants.NC_CONTIGS
    ) -> "Contigs":
        """Load :class:`Contigs` from specified path

        Parameters
        ----------
        path: Union[str, Path]
            Path where contigs are stored in zarr format
        group: str
            group in zarr store where contigs are stored. Default is where
            contigs are stored for splicegraphs, but some stores (e.g. for
            :class:`SJExperiment`) need to hold two independent sets of
            contigs, which requires having two separate group names, which this
            allows
        """
        with xr.open_zarr(path, group=group) as df:
            return Contigs.from_list(df.seqid.values.tolist())
