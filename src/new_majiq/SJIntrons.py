"""
SJIntrons.py

contig regions not intersected by gene exons (stranded or unstranded), noting
if intersecting annotated intron or not (part of annotated exon)

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.Contigs import Contigs
from new_majiq.Exons import Exons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.ContigRegions import ContigRegions
from new_majiq.internals import SJIntrons as _SJIntrons
from typing import (
    Optional,
    Union,
)
from pathlib import Path


class SJIntrons(ContigRegions):
    def __init__(self, sj_introns: _SJIntrons):
        super().__init__(sj_introns)
        return

    @property
    def _sj_introns(self) -> _SJIntrons:
        """expose underlying internals representation of SJIntrons"""
        return self._contig_regions

    @property
    def si_idx(self) -> np.ndarray:
        return self._region_idx

    @property
    def annotated(self) -> np.ndarray:
        """Indicates if region associated with annotated intron"""
        return self._sj_introns.annotated

    @property
    def df(self) -> xr.Dataset:
        """view on underlying SJ introns as xarray Dataset"""
        return xr.Dataset(
            {},
            {
                "si_idx": self.si_idx,
                "contig_idx": ("si_idx", self.contig_idx),
                "start": ("si_idx", self.start),
                "end": ("si_idx", self.end),
                "strand": ("si_idx", self.strand),
                "annotated": ("si_idx", self.annotated),
            },
        )

    def to_zarr(self, path: Union[str, Path], mode: str) -> None:
        """Serialize to zarr format. Note contigs need to be saved separately"""
        self.df.drop_vars("si_idx").to_zarr(
            path, mode=mode, group=constants.NC_SJINTRONS
        )
        return

    @classmethod
    def from_exons_and_introns(
        cls, exons: Exons, gene_introns: GeneIntrons, stranded: bool
    ) -> "SJIntrons":
        """Construct SJIntrons from exons and gene introns"""
        return SJIntrons(
            _SJIntrons.from_exons_and_introns(
                exons._exons, gene_introns._gene_introns, stranded
            )
        )

    @classmethod
    def from_contigs(cls, contigs: Contigs) -> "SJIntrons":
        """Empty SJIntrons linked to specified contigs"""
        return SJIntrons(_SJIntrons(contigs._contigs, [], [], [], [], []))

    @classmethod
    def from_zarr(
        cls,
        path: Union[str, Path],
        contigs: Optional[Contigs] = None,
    ) -> "SJIntrons":
        """Read SJIntrons from zarr file

        Parameters
        ----------
        path: Union[str, Path]
            path to zarr file
        contigs: Optional[Contigs]
            contigs on which the introns are defined. If None, try loading from
            zarr file.
        """
        if contigs is None:
            contigs = Contigs.from_zarr(path)
        with xr.open_zarr(path, group=constants.NC_SJINTRONS) as df:
            return SJIntrons(
                _SJIntrons(
                    contigs._contigs,
                    df.contig_idx.values,
                    df.start.values,
                    df.end.values,
                    df.strand.values,
                    df.annotated.values,
                )
            )
