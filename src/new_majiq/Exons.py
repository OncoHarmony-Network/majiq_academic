"""
Exons.py

Exons for genes

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.Genes import Genes
from new_majiq.GeneRegions import GeneRegions
from new_majiq.internals import Exons as _Exons
from typing import (
    Optional,
    TYPE_CHECKING,
    Union,
)
from pathlib import Path

if TYPE_CHECKING:
    from new_majiq.GeneJunctions import GeneJunctions


class Exons(GeneRegions):
    def __init__(self, exons: _Exons):
        super().__init__(exons)
        return

    def infer_with_junctions(self, junctions: "GeneJunctions") -> "Exons":
        """Infer denovo exons/extended exon boundaries given denovo junctions"""
        from new_majiq.internals import SpliceGraph as _SpliceGraph
        return Exons(
            _SpliceGraph.infer_exons(self._exons, junctions._gene_junctions)
        )

    def checksum(self):
        return self._exons.checksum()

    @property
    def _exons(self) -> _Exons:
        """expose underlying internals representation of Exons"""
        return self._gene_regions

    @property
    def exon_idx(self) -> np.ndarray:
        return self._region_idx

    @property
    def annotated_start(self) -> np.ndarray:
        """Annotated coordinates start. If denovo, -1

        Note that denovo behavior is different than previous versions of MAJIQ
        """
        return self._exons.annotated_start

    @property
    def annotated_end(self) -> np.ndarray:
        """Annotated coordinates end. If denovo, -1

        Note that denovo behavior is different than previous versions of MAJIQ
        """
        return self._exons.annotated_end

    def is_denovo(self, exon_idx: Optional[np.ndarray] = None) -> np.ndarray:
        """Return denovo status of exon"""
        if exon_idx is None:
            exon_idx = self.exon_idx
        return self._exons.is_denovo(exon_idx)

    @property
    def df(self) -> xr.Dataset:
        """view on underlying exons as xarray Dataset"""
        return xr.Dataset(
            {},
            {
                "exon_idx": self.exon_idx,
                "gene_idx": ("exon_idx", self.gene_idx),
                "start": ("exon_idx", self.start),
                "end": ("exon_idx", self.end),
                "annotated_start": ("exon_idx", self.annotated_start),
                "annotated_end": ("exon_idx", self.annotated_end),
            },
        )

    def to_netcdf(self, path: Union[str, Path], mode: str) -> None:
        """Serialize to netcdf format. Note genes need to be saved separately"""
        self.df.drop_vars("exon_idx").to_netcdf(path, mode, group=constants.NC_EXONS)
        return

    @classmethod
    def from_netcdf(
        self,
        path: Union[str, Path],
        genes: Optional[Genes] = None,
    ) -> "Exons":
        """Read exons from netcdf file

        Parameters
        ----------
        path: Union[str, Path]
            path to netcdf file
        genes: Optional[Genes]
            genes on which the exons are defined. If None, try loading from
            netcdf file. Note that new_majiq checks if objects refer to the
            same genes (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        df = xr.open_dataset(path, group=constants.NC_EXONS)
        if genes is None:
            genes = Genes.from_netcdf(path)
        return Exons(
            _Exons(
                genes._genes,
                df.gene_idx,
                df.start,
                df.end,
                df.annotated_start,
                df.annotated_end,
            )
        )
