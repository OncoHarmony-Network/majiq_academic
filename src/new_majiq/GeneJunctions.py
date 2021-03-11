"""
GeneJunctions.py

Container of gene junctions that are open intervals over genes

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from typing import (
    Optional,
    Union,
)
from new_majiq.internals import GeneJunctions as _GeneJunctions
from new_majiq.GeneConnections import GeneConnections
from new_majiq.Genes import Genes
from pathlib import Path


class GeneJunctions(GeneConnections):
    def __init__(self, gene_junctions: _GeneJunctions):
        super().__init__(gene_junctions)
        return

    def checksum(self):
        return self._gene_junctions.checksum()

    @property
    def _gene_junctions(self) -> _GeneJunctions:
        """Underlying internals representation"""
        return self._gene_connections

    @property
    def gj_idx(self) -> np.ndarray:
        return self._region_idx

    @property
    def df(self) -> xr.Dataset:
        """xr.Dataset view of gene junctions data"""
        return xr.Dataset(
            {},
            {
                "gj_idx": self.gj_idx,
                "gene_idx": ("gj_idx", self.gene_idx),
                "start": ("gj_idx", self.start),
                "end": ("gj_idx", self.end),
                "denovo": ("gj_idx", self.denovo),
                "passed_build": ("gj_idx", self.passed_build),
                "simplified": ("gj_idx", self.simplified),
                "start_exon_idx": ("gj_idx", self.start_exon_idx),
                "end_exon_idx": ("gj_idx", self.end_exon_idx),
            },
        )

    def to_netcdf(self, path: Union[str, Path], mode: str) -> None:
        """Serialize to netcdf format. Note genes need to be saved separately"""
        self.df.drop_vars("gj_idx").to_netcdf(path, mode, group=constants.NC_GENEJUNCTIONS)
        return

    @classmethod
    def from_netcdf(
        self,
        path: Union[str, Path],
        genes: Optional[Genes] = None,
    ) -> "GeneJunctions":
        """Read exons from netcdf file

        Parameters
        ----------
        path: Union[str, Path]
            path to netcdf file
        genes: Optional[Genes]
            genes on which the junctions are defined. If None, try loading from
            netcdf file. Note that new_majiq checks if objects refer to the
            same genes (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        df = xr.open_dataset(path, group=constants.NC_GENEJUNCTIONS)
        if genes is None:
            genes = Genes.from_netcdf(path)
        return GeneJunctions(
            _GeneJunctions(
                genes._genes,
                df.gene_idx,
                df.start,
                df.end,
                df.denovo,
                df.passed_build,
                df.simplified,
            )
        )
