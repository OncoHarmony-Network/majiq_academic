"""
Exons.py

Exons for genes

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import numpy as np
import xarray as xr

import new_majiq.constants as constants
from new_majiq._workarounds import _load_zerodim_variables
from new_majiq.GeneRegions import GeneRegions
from new_majiq.Genes import Genes
from new_majiq.internals import Exons as _Exons

if TYPE_CHECKING:
    from new_majiq.GeneIntrons import GeneIntrons
    from new_majiq.GeneJunctions import GeneJunctions


class Exons(GeneRegions):
    def __init__(self, exons: _Exons):
        super().__init__(exons)
        return

    def infer_with_junctions(self, junctions: "GeneJunctions") -> "Exons":
        """Infer denovo exons/extended exon boundaries given denovo junctions"""
        from new_majiq.internals import SpliceGraph as _SpliceGraph

        return Exons(_SpliceGraph.infer_exons(self._exons, junctions._gene_junctions))

    def checksum(self):
        return self._exons.checksum()

    @property
    def _exons(self) -> _Exons:
        """expose underlying internals representation of Exons"""
        return self._gene_regions

    @property
    def exon_idx(self) -> np.ndarray:
        return self._region_idx

    def potential_introns(
        self, make_simplified: bool = constants.DEFAULT_BUILD_DENOVO_SIMPLIFIED
    ) -> "GeneIntrons":
        """denovo, nonpassed introns corresponding to these exons

        Parameters
        ----------
        make_simplified: bool
            Simplified status of resulting introns
        """
        from new_majiq.GeneIntrons import GeneIntrons

        return GeneIntrons(self._exons.potential_introns(make_simplified))

    def empty_introns(self) -> "GeneIntrons":
        """Get empty introns to go with exons

        Note
        ----
        Currently doing unncessary work by inferring all introns.
        But, it gets the job done since potential_introns starts off denovo and
        not passed
        """
        return self.potential_introns().filter_passed()

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

    def is_exon_extension(self, exon_idx: Optional[np.ndarray] = None) -> np.ndarray:
        """Return if exon(s) have exon extension"""
        if exon_idx is None:
            exon_idx = self.exon_idx
        return self._exons.is_exon_extension(exon_idx)

    def is_full_exon(self, exon_idx: Optional[np.ndarray] = None) -> np.ndarray:
        """Return if exon(s) are full exons"""
        if exon_idx is None:
            exon_idx = self.exon_idx
        return self._exons.is_full_exon(exon_idx)

    def is_half_exon(self, exon_idx: Optional[np.ndarray] = None) -> np.ndarray:
        """Return if exon(s) are half exons"""
        if exon_idx is None:
            exon_idx = self.exon_idx
        return self._exons.is_half_exon(exon_idx)

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

    def to_zarr(
        self, path: Union[str, Path], mode: str, consolidated: bool = True
    ) -> None:
        """Serialize to zarr format. Note genes need to be saved separately"""
        self.df.drop_vars("exon_idx").pipe(lambda x: x.chunk(x.sizes)).pipe(
            _load_zerodim_variables
        ).to_zarr(path, mode=mode, group=constants.NC_EXONS, consolidated=consolidated)
        return

    @classmethod
    def from_zarr(
        self,
        path: Union[str, Path],
        genes: Optional[Genes] = None,
    ) -> "Exons":
        """Read exons from zarr file

        Parameters
        ----------
        path: Union[str, Path]
            path to zarr file
        genes: Optional[Genes]
            genes on which the exons are defined. If None, try loading from
            zarr file. Note that new_majiq checks if objects refer to the
            same genes (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        if genes is None:
            genes = Genes.from_zarr(path)
        with xr.open_zarr(path, group=constants.NC_EXONS) as df:
            return Exons(
                _Exons(
                    genes._genes,
                    df.gene_idx.values,
                    df.start.values,
                    df.end.values,
                    df.annotated_start.values,
                    df.annotated_end.values,
                )
            )
