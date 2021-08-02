"""
SpliceGraph.py

Combined regions of contigs containing genes containing exons connected by
introns and junctions

Author: Joseph K Aicher
"""

import new_majiq.constants as constants

from new_majiq.Contigs import Contigs
from new_majiq.Genes import Genes
from new_majiq.Exons import Exons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.ExonConnections import ExonConnections
from new_majiq.internals import SpliceGraph as _SpliceGraph
from new_majiq.internals import GFF3Types

from typing import (
    Final,
    Optional,
    Union,
)
from pathlib import Path


class SpliceGraph(object):
    def __init__(self, sg: _SpliceGraph):
        self._sg: Final[_SpliceGraph] = sg
        return

    @classmethod
    def from_components(
        cls,
        contigs: Contigs,
        genes: Genes,
        exons: Exons,
        junctions: GeneJunctions,
        introns: GeneIntrons,
    ) -> "SpliceGraph":
        return SpliceGraph(
            _SpliceGraph(
                contigs._contigs,
                genes._genes,
                exons._exons,
                junctions._gene_junctions,
                introns._gene_introns,
            )
        )

    def with_updated_exon_connections(
        self, exon_connections: ExonConnections
    ) -> "SpliceGraph":
        """convenience function to build splicegraph with same contigs/genes"""
        return SpliceGraph.from_components(
            self.contigs,
            self.genes,
            exon_connections.exons,
            exon_connections.junctions,
            exon_connections.introns,
        )

    @classmethod
    def from_gff3(
        cls,
        path: Union[str, Path],
        process_ir: bool = constants.DEFAULT_BUILD_PROCESS_IR,
        gff3_types: GFF3Types = constants.DEFAULT_BUILD_GFF3TYPES,
    ) -> "SpliceGraph":
        return SpliceGraph(_SpliceGraph.from_gff3(str(path), process_ir, gff3_types))

    @property
    def contigs(self) -> Contigs:
        return Contigs(self._sg._contigs)

    @property
    def genes(self) -> Genes:
        return Genes(self._sg._genes)

    @property
    def exons(self) -> Exons:
        return Exons(self._sg._exons)

    @property
    def introns(self) -> GeneIntrons:
        return GeneIntrons(self._sg._introns)

    @property
    def junctions(self) -> GeneJunctions:
        return GeneJunctions(self._sg._junctions)

    @property
    def exon_connections(self) -> ExonConnections:
        return ExonConnections(self._sg._exon_connections)

    def to_zarr(self, path: Union[str, Path]) -> None:
        """Serialize splicegraph to zarr format"""
        if Path(path).exists():
            raise ValueError(
                f"Will not save splicegraph to already existing file {path}."
                " Please delete and try again if you want it there, otherwise"
                " pick a different output path"
            )
        self.contigs.to_zarr(path, "w")
        self.genes.to_zarr(path, "a")
        self.exons.to_zarr(path, "a")
        self.introns.to_zarr(path, "a")
        self.junctions.to_zarr(path, "a")
        return

    @classmethod
    def from_zarr(
        cls, path: Union[str, Path], genes: Optional[Genes] = None
    ) -> "SpliceGraph":
        """Load SpliceGraph from specified path"""
        if genes is None:
            contigs = Contigs.from_zarr(path)
            genes = Genes.from_zarr(path, contigs)
        else:
            # make sure that genes match what are in file
            if genes != Genes.from_zarr(path):
                raise ValueError(
                    "Optionally provied genes do not match those from path"
                )
            contigs = genes.contigs
        exons = Exons.from_zarr(path, genes)
        introns = GeneIntrons.from_zarr(path, genes)
        junctions = GeneJunctions.from_zarr(path, genes)
        return SpliceGraph(
            _SpliceGraph(
                contigs._contigs,
                genes._genes,
                exons._exons,
                junctions._gene_junctions,
                introns._gene_introns,
            )
        )
