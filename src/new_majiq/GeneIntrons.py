"""
GeneIntrons.py

Container of gene introns that are open intervals over genes

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from typing import (
    Optional,
    Union,
)
from new_majiq.internals import GeneIntrons as _GeneIntrons
from new_majiq.GeneConnections import GeneConnections
from new_majiq.Genes import Genes
from pathlib import Path


class GeneIntrons(GeneConnections):
    def __init__(self, gene_introns: _GeneIntrons):
        super().__init__(gene_introns)
        return

    @property
    def _gene_introns(self) -> _GeneIntrons:
        """Underlying internals representation"""
        return self._gene_connections

    @property
    def gi_idx(self) -> np.ndarray:
        return self._region_idx

    @property
    def df(self) -> xr.Dataset:
        """xr.Dataset view of gene introns data"""
        return xr.Dataset(
            {},
            {
                "gi_idx": self.gi_idx,
                "gene_idx": ("gi_idx", self.gene_idx),
                "start": ("gi_idx", self.start),
                "end": ("gi_idx", self.end),
                "denovo": ("gi_idx", self.denovo),
                "passed_build": ("gi_idx", self.passed_build),
                "simplified": ("gi_idx", self.simplified),
                "start_gi_idx": ("gi_idx", self.start_exon_idx),
                "end_gi_idx": ("gi_idx", self.end_exon_idx),
            },
        )

    def to_netcdf(self, path: Union[str, Path], mode: str) -> None:
        """Serialize to netcdf format. Note genes need to be saved separately"""
        self.df.to_netcdf(path, mode, group=constants.NC_GENEINTRONS)
        return

    @classmethod
    def from_netcdf(
        self,
        path: Union[str, Path],
        genes: Optional[Genes] = None,
    ) -> "GeneIntrons":
        """Read exons from netcdf file

        Parameters
        ----------
        path: Union[str, Path]
            path to netcdf file
        genes: Optional[Genes]
            genes on which the introns are defined. If None, try loading from
            netcdf file. Note that new_majiq checks if objects refer to the
            same genes (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        df = xr.open_dataset(path, group=constants.NC_GENEINTRONS)
        if genes is None:
            genes = Genes.from_netcdf(path)
        return GeneIntrons(
            _GeneIntrons(
                genes._genes,
                df.gene_idx,
                df.start,
                df.end,
                df.denovo,
                df.passed_build,
                df.simplified,
            )
        )
