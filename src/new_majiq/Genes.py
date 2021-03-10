"""
Genes.py

Container of genes that are closed intervals over different contigs

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from typing import (
    List,
    Optional,
    Union,
)
from new_majiq.internals import Genes as _Genes
from new_majiq.Contigs import Contigs
from new_majiq.ContigRegions import ContigRegions
from pathlib import Path


class Genes(ContigRegions):
    """Genes defined on contigs used by a splicegraph"""

    def __init__(self, genes: _Genes):
        super().__init__(genes)
        return

    @property
    def _genes(self) -> _Genes:
        """expose underlying internals representation of Genes"""
        return self._contig_regions

    @property
    def gene_idx(self) -> np.ndarray:
        return self._region_idx

    def __getitem__(self, gene_id: str) -> Optional[int]:
        """Get gene_idx for specified gene_id"""
        raise NotImplementedError("Need to expose safe_idx")

    @property
    def gene_id(self) -> List[str]:
        return self._genes.gene_id

    @property
    def gene_name(self) -> List[str]:
        return self._genes.gene_name

    def annotate_gene_idx(self, df: xr.Dataset) -> xr.Dataset:
        """For now, just add gene_id to df using df.gene_idx"""
        return df.assign_coords(
            gene_id=(df.gene_idx.dims, np.array(self.gene_id)[df.gene_idx]),
        )

    @property
    def df(self) -> xr.Dataset:
        """view of underlying genes object as xarray Dataset"""
        return xr.Dataset(
            {},
            {
                "gene_idx": self.gene_idx,
                "contig_idx": ("gene_idx", self.contig_idx),
                "start": ("gene_idx", self.start),
                "end": ("gene_idx", self.end),
                "strand": ("gene_idx", self.strand),
                "gene_id": ("gene_idx", self.gene_id),
                "gene_name": ("gene_idx", self.gene_name),
            },
        )

    def to_netcdf(self, path: Union[str, Path], mode: str) -> None:
        """Serialize to netcdf format. Note contigs need to be saved separately"""
        self.df.to_netcdf(path, mode, group=constants.NC_GENES)
        return

    @classmethod
    def from_netcdf(
        self,
        path: Union[str, Path],
        contigs: Optional[Contigs] = None,
    ) -> "Genes":
        """Read genes from netcdf file.

        Parameters
        ----------
        path: Union[str, Path]
            path to netcdf file
        contigs: Optional[Contigs]
            contigs on which the genes are defined. If None, try loading from
            netcdf file. Note that new_majiq checks if objects refer to the
            same contigs (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        df = xr.open_dataset(path, group=constants.NC_GENES)
        if contigs is None:
            contigs = Contigs.from_netcdf(path)
        return Genes(
            _Genes(
                contigs._contigs,
                df.contig_idx,
                df.start,
                df.end,
                df.strand,
                df.gene_id.values.tolist(),
                df.gene_name.values.tolist(),
            )
        )
