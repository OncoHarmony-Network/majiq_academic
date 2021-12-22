"""
GeneJunctions.py

Container of gene junctions that are open intervals over genes

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import new_majiq.constants as constants
from new_majiq.internals import GeneJunctions as _GeneJunctions

from ._workarounds import _load_zerodim_variables
from .Exons import Exons
from .GeneConnections import GeneConnections
from .Genes import Genes

if TYPE_CHECKING:
    from .PassedJunctions import GroupJunctionsGenerator, PassedJunctionsGenerator


class GeneJunctions(GeneConnections):
    """Collection of junctions per gene and their coordinates, flags, and exons"""

    def __init__(self, gene_junctions: _GeneJunctions):
        super().__init__(gene_junctions)
        return

    def build_group(self, exons: Exons) -> "GroupJunctionsGenerator":
        """Create :py:class:`GroupJunctionsGenerator` starting from these junctions and exons

        Parameters
        ----------
        exons: Exons
            Exons over the same genes as the junctions, which enable
            identification of most likely genes to which novel junctions belong
        """
        from .PassedJunctions import GroupJunctionsGenerator

        return GroupJunctionsGenerator(self, exons)

    def builder(self) -> "PassedJunctionsGenerator":
        """Create :py:class:`PassedJunctionsGenerator` starting from these junctions"""
        from .PassedJunctions import PassedJunctionsGenerator

        return PassedJunctionsGenerator(self)

    @property
    def _gene_junctions(self) -> _GeneJunctions:
        """Underlying internals representation"""
        return self._gene_connections

    @property
    def gj_idx(self) -> npt.NDArray[np.int64]:
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

    def to_zarr(
        self, path: Union[str, Path], mode: str, consolidated: bool = True
    ) -> None:
        """Serialize to zarr format. Note genes need to be saved separately"""
        self.df.drop_vars(["gj_idx", "start_exon_idx", "end_exon_idx"]).pipe(
            lambda x: x.chunk(x.sizes)
        ).pipe(_load_zerodim_variables).to_zarr(
            path,
            mode=mode,
            group=constants.NC_GENEJUNCTIONS,
            consolidated=consolidated,
        )
        return

    @classmethod
    def from_arrays(
        cls,
        genes: Genes,
        gene_idx: npt._ArrayLikeInt_co,
        start: npt._ArrayLikeInt_co,
        end: npt._ArrayLikeInt_co,
        denovo: Optional[npt._ArrayLikeBool_co] = None,
        passed_build: Optional[npt._ArrayLikeBool_co] = None,
        simplified: Optional[npt._ArrayLikeBool_co] = None,
        connected_exons: Optional[Exons] = None,
    ) -> "GeneJunctions":
        """Create :class:`GeneJunctions` from :class:`Genes` and input arrays"""
        if denovo is None:
            denovo = np.zeros_like(gene_idx, dtype=np.bool_)
        if passed_build is None:
            passed_build = np.zeros_like(gene_idx, dtype=np.bool_)
        if simplified is None:
            simplified = np.zeros_like(gene_idx, dtype=np.bool_)
        return GeneJunctions(
            _GeneJunctions(
                genes._genes,
                gene_idx,
                start,
                end,
                denovo,
                passed_build,
                simplified,
                connected_exons=None
                if connected_exons is None
                else connected_exons._exons,
            )
        )

    @staticmethod
    def from_zarr(
        path: Union[str, Path],
        genes: Optional[Genes] = None,
    ) -> "GeneJunctions":
        """Load :py:class:`GeneJunctions` from specified path

        Parameters
        ----------
        path: Union[str, Path]
            path to zarr file
        genes: Optional[Genes]
            genes on which the junctions are defined. If None, try loading from
            zarr file. Note that new_majiq checks if objects refer to the
            same genes (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        if genes is None:
            genes = Genes.from_zarr(path)
        with xr.open_zarr(path, group=constants.NC_GENEJUNCTIONS) as df:
            df.load()
            return GeneJunctions.from_arrays(
                genes,
                df.gene_idx.values,
                df.start.values,
                df.end.values,
                denovo=df.denovo.values,
                passed_build=df.passed_build.values,
                simplified=df.simplified.values,
            )
