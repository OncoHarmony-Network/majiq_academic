"""
GeneJunctions.py

Container of gene junctions that are open intervals over genes

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import numpy as np
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

    @staticmethod
    def load_dataset(path: Union[str, Path]) -> xr.Dataset:
        """Load junctions from zarr file as :py:class:`xr.Dataset`

        Load junctions from zarr file as :py:class:`xr.Dataset`. This is a more
        lightweight representation of junctions when it is not necessary to
        connect to gene or contig information

        Parameters
        ----------
        path: Union[str, Path]
            Path to where :py:class:`GeneJunctions` are saved (this includes
            splicegraphs)

        Returns
        -------
        xr.Dataset
        """
        with xr.open_zarr(path, group=constants.NC_GENEJUNCTIONS) as df:
            return df.load()

    @staticmethod
    def from_dataset_and_genes(df: xr.Dataset, genes: Genes) -> "GeneJunctions":
        """Create :py:class:`GeneJunctions` from junction dataset and :py:class:`Genes`

        Parameters
        ----------
        df: xr.Dataset
            Variables:
                gene_idx(gj_idx) uint64
                start(gj_idx) int64
                end(gj_idx) int64
                denovo(gj_idx) bool
                passed_build(gj_idx) bool
                simplified(gj_idx) bool
        genes: Genes
            Genes matched to gene_idx in input dataset

        Returns
        -------
        GeneJunctions
        """
        return GeneJunctions(
            _GeneJunctions(
                genes._genes,
                df.gene_idx.values,
                df.start.values,
                df.end.values,
                df.denovo.values,
                df.passed_build.values,
                df.simplified.values,
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
        return GeneJunctions.from_dataset_and_genes(
            GeneJunctions.load_dataset(path), genes
        )
