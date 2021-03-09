"""
SpliceGraph.py

Combined regions of contigs containing genes containing exons connected by
introns and junctions

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.Contigs import Contigs
from new_majiq.Genes import Genes
from new_majiq.Exons import Exons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.internals import SpliceGraph as _SpliceGraph
from new_majiq.internals import GFF3Types

from typing import (
    Final,
    Union,
)
from pathlib import Path


class SpliceGraph(object):
    def __init__(self, sg: _SpliceGraph):
        self._sg: Final[_SpliceGraph] = sg
        return

    @classmethod
    def from_gff3(
        cls,
        path: Union[str, Path],
        process_ir: bool = constants.DEFAULT_BUILD_PROCESS_IR,
        gff3_types: GFF3Types = constants.DEFAULT_BUILD_GFF3TYPES,
    ) -> "SpliceGraph":
        return SpliceGraph(_SpliceGraph.from_gff3(str(path), process_ir, gff3_types))

    @property
    def _contigs(self) -> Contigs:
        return Contigs(self._sg._contigs)

    @property
    def _genes(self) -> Genes:
        return Genes(self._sg._genes)

    @property
    def _exons(self) -> Exons:
        return Exons(self._sg._exons)

    @property
    def _introns(self) -> GeneIntrons:
        return GeneIntrons(self._sg._introns)

    @property
    def _junctions(self) -> GeneJunctions:
        return GeneJunctions(self._sg._junctions)

    def to_netcdf(self, path: Union[str, Path]) -> None:
        """Serialize splicegraph to netcdf format"""
        if Path(path).exists():
            raise ValueError(
                f"Will not save splicegraph to already existing file {path}."
                " Please delete and try again if you want it there, otherwise"
                " pick a different output path"
            )
        self._contigs.to_netcdf(path, "w")
        self._genes.to_netcdf(path, "a")
        self._exons.to_netcdf(path, "a")
        self._introns.to_netcdf(path, "a")
        self._junctions.to_netcdf(path, "a")
        return

    @classmethod
    def from_netcdf(cls, path: Union[str, Path]) -> "SpliceGraph":
        """Load SpliceGraph from specified path"""
        contigs = Contigs.from_netcdf(path)
        genes = Genes.from_netcdf(path, contigs)
        exons = Exons.from_netcdf(path, genes)
        introns = GeneIntrons.from_netcdf(path, genes)
        junctions = GeneJunctions.from_netcdf(path, genes)
        return SpliceGraph(
            _SpliceGraph(
                contigs._contigs,
                genes._genes,
                exons._exons,
                junctions._gene_junctions,
                introns._gene_introns,
            )
        )
