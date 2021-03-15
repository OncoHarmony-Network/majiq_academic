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
    TYPE_CHECKING,
    Union,
)
from new_majiq.internals import GeneIntrons as _GeneIntrons
from new_majiq.GeneConnections import GeneConnections
from new_majiq.Genes import Genes
from new_majiq.Exons import Exons
from pathlib import Path

if TYPE_CHECKING:
    from new_majiq.GroupIntronsGenerator import GroupIntronsGenerator


class GeneIntrons(GeneConnections):
    def __init__(self, gene_introns: _GeneIntrons):
        super().__init__(gene_introns)
        return

    def build_group(self) -> "GroupIntronsGenerator":
        """Create build group to update these introns to be passed or not"""
        from new_majiq.GroupIntronsGenerator import GroupIntronsGenerator

        return GroupIntronsGenerator(self)

    def filter_passed(
        self,
        keep_annotated: bool = constants.DEFAULT_BUILD_KEEP_ANNOTATED_IR,
        discard_denovo: bool = not constants.DEFAULT_BUILD_DENOVO_IR,
    ) -> "GeneIntrons":
        """Get subset of introns that passed build filters

        Parameters
        ----------
        keep_annotated: bool
            Keep all annotated introns regardless of whether they passed
        discard_denovo: bool
            Discard all denovo introns regardless of whether they passed
        """
        return GeneIntrons(
            self._gene_introns.filter_passed(keep_annotated, discard_denovo)
        )

    def potential_introns(self, exons: Exons) -> "GeneIntrons":
        """Get potential gene introns from exons, noting annotated introns"""
        return GeneIntrons(self._gene_introns.potential_introns(exons._exons))

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
                "start_exon_idx": ("gi_idx", self.start_exon_idx),
                "end_exon_idx": ("gi_idx", self.end_exon_idx),
            },
        )

    def to_netcdf(self, path: Union[str, Path], mode: str) -> None:
        """Serialize to netcdf format. Note genes need to be saved separately"""
        self.df.drop_vars("gi_idx").to_netcdf(
            path, mode, group=constants.NC_GENEINTRONS
        )
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
        if genes is None:
            genes = Genes.from_netcdf(path)
        with xr.open_dataset(path, group=constants.NC_GENEINTRONS) as df:
            return GeneIntrons(
                _GeneIntrons(
                    genes._genes,
                    df.gene_idx.values,
                    df.start.values,
                    df.end.values,
                    df.denovo.values,
                    df.passed_build.values,
                    df.simplified.values,
                )
            )
