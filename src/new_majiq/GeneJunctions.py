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
    Sequence,
    TYPE_CHECKING,
    Union,
)
from new_majiq.internals import GeneJunctions as _GeneJunctions
from new_majiq.GeneConnections import GeneConnections
from new_majiq.Genes import Genes
from new_majiq.Exons import Exons
from pathlib import Path

if TYPE_CHECKING:
    from new_majiq.PassedJunctions import (
        GroupJunctionsGenerator,
        PassedJunctionsGenerator,
    )


class GeneJunctions(GeneConnections):
    def __init__(self, gene_junctions: _GeneJunctions):
        super().__init__(gene_junctions)
        return

    def build_group(self, exons: Exons) -> "GroupJunctionsGenerator":
        """Create accumulator of per-experiment passed junctions for build group"""
        from new_majiq.PassedJunctions import GroupJunctionsGenerator

        return GroupJunctionsGenerator(self, exons)

    def builder(self) -> "PassedJunctionsGenerator":
        """Create accumulator of passed junctions from build groups"""
        from new_majiq.PassedJunctions import PassedJunctionsGenerator

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

    def to_zarr(self, path: Union[str, Path], mode: str) -> None:
        """Serialize to zarr format. Note genes need to be saved separately"""
        self.df.drop_vars(["gj_idx", "start_exon_idx", "end_exon_idx"]).to_zarr(
            path,
            mode=mode,
            group=constants.NC_GENEJUNCTIONS,
        )
        return

    @staticmethod
    def load_dataset(path: Union[str, Path]) -> xr.Dataset:
        """load junctions table from file"""
        with xr.open_zarr(path, group=constants.NC_GENEJUNCTIONS) as df:
            return df.load()

    @staticmethod
    def _combine_datasets(dfs: Sequence[xr.Dataset]) -> xr.Dataset:
        """combine junctions found in multiple splicegraphs

        Assumes that they all share the same genes
        """
        return (
            xr.concat(
                [
                    df
                    # work with not_passed so that can combine using ALL
                    .assign_coords(not_passed=~df.passed_build)
                    # set index so can see matches
                    .set_index(gj_idx=["gene_idx", "start", "end"])
                    # select variables we want to use
                    .reset_coords()[["denovo", "not_passed", "simplified"]]
                    for df in dfs
                ],
                # concatenate over new dimension inputs
                dim="inputs",
                # take all unique indexes (gene_idx, start, end)
                join="outer",
                # missing values ~ is denovo, not passed, and simplified
                fill_value=True,
                # don't retain any attributes here
                combine_attrs="drop",
            )
            # aggregate over inputs using ALL rule
            .all("inputs")
            # get coordinates we want back ~ from load_dataset
            .assign_coords(passed_build=lambda df: ~df.not_passed)
            .drop_vars(["not_passed"])
            .set_coords(["denovo", "simplified"])
            .reset_index("gj_idx")
        )

    @staticmethod
    def from_dataset_and_genes(df: xr.Dataset, genes: Genes) -> "GeneJunctions":
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
        """Read exons from zarr file

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
