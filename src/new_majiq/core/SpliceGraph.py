"""
SpliceGraph.py

Combined regions of contigs containing genes containing exons connected by
introns and junctions

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Any, Callable, Final, Optional, Union

import new_majiq.constants as constants
from new_majiq.internals import SpliceGraph as _SpliceGraph

from .Contigs import Contigs
from .ExonConnections import ExonConnections
from .Exons import Exons
from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions
from .Genes import Genes
from .GFF3TypesMap import GFF3TypesMap


class SpliceGraph(object):
    """Representation of all possible splicing changes in each gene.

    Representation of splicing in each gene as exons connected by spliced
    junctions and retained introns.
    This representation is composed of:

    - the collection of all :class:`Contigs` (e.g. chromosomes) on which genes
      can be defined (:attr:`SpliceGraph.contigs`),
    - the collection of all :class:`Genes` and their coordinates on the
      splicegraph contigs (:attr:`SpliceGraph.genes`),
    - the collection of all :class:`Exons` for each gene
      (:attr:`SpliceGraph.exons`),
    - the collection of all :class:`GeneIntrons` for each gene
      (:attr:`SpliceGraph.introns`),
    - the collection of all :class:`GeneJunctions` for each gene
      (:attr:`SpliceGraph.junctions`),
    - :class:`ExonConnections` associating introns and junctions to source and
      target exons, enabling the identification of splicing events, e.g. LSVs
      (:attr:`SpliceGraph.exon_connections`).

    Parameters
    ----------
    sg: _SpliceGraph
        Underlying object binding the internal C++ API

    See Also
    --------
    SpliceGraph.from_gff3 : Initialize :class:`SpliceGraph` from GFF3 file
    SpliceGraph.from_components : Construct :class:`SpliceGraph` from components
    SpliceGraph.from_zarr : Load :class:`SpliceGraph` saved in Zarr format
    """

    def __init__(self, sg: _SpliceGraph):
        """Construct :class:`SpliceGraph` using object from internal C++ API

        Parameters
        ----------
        sg: _SpliceGraph
            Underlying object binding the internal C++ API
        """
        self._sg: Final[_SpliceGraph] = sg
        return

    def __repr__(self) -> str:
        return (
            f"SpliceGraph["
            f"{len(self.contigs)} contigs,"
            f" {len(self.genes)} genes,"
            f" {len(self.exons)}/{len(self.introns)}/{len(self.junctions)}"
            " exons/introns/junctions]"
        )

    @classmethod
    def from_components(
        cls,
        contigs: Contigs,
        genes: Genes,
        exons: Exons,
        junctions: GeneJunctions,
        introns: GeneIntrons,
    ) -> "SpliceGraph":
        """Construct :py:class:`SpliceGraph` with given components

        Parameters
        ----------
        contigs: Contigs
        genes: Genes
        exons: Exons
        junctions: GeneJunctions
        introns: GeneIntrons

        Returns
        -------
        SpliceGraph
        """
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
        """Create :py:class:`SpliceGraph` from exon connections with same genes

        Parameters
        ----------
        exon_connections: ExonConnections

        Returns
        -------
        SpliceGraph
        """
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
        gff3_types: Optional[GFF3TypesMap] = None,
        log_function: Optional[Callable[[str, str, int], Any]] = None,
    ) -> "SpliceGraph":
        """Create :py:class:`SpliceGraph` from GFF3 transcriptome annotations

        Parameters
        ----------
        path: Union[str, Path]
            Path to GFF3 file (can be gzipped)
        process_ir: bool
            Identify annotated introns. This should generally be True
        gff3_types: Optional[GFF3TypesMap]
            How GFF3 lines will be parsed hierarchically for genes, transcript,
            and exon definitions. If None, use default initialization of
            :class:`GFF3TypesMap`.
        log_function: Optional[Callable[str, str, int]]
            GFF3 types that were not accepted or explicitly ignored will be
            reported per unique type/part of hierarchy where they were missing,
            along with their counts. This function will be called for each
            unique type: log_function(type, location_in_hierarchy, count).
        """
        if gff3_types is None:
            gff3_types = GFF3TypesMap()
        return SpliceGraph(
            _SpliceGraph.from_gff3(
                str(path), process_ir, gff3_types.current_map, log_function
            )
        )

    @property
    def contigs(self) -> Contigs:
        """The collection of all :class:`Contigs` on which genes can be defined"""
        return Contigs(self._sg._contigs)

    @property
    def genes(self) -> Genes:
        """The collection of all :class:`Genes` and their coordinates on contigs"""
        return Genes(self._sg._genes)

    @property
    def exons(self) -> Exons:
        """The collection of all :class:`Exons` for each gene"""
        return Exons(self._sg._exons)

    @property
    def introns(self) -> GeneIntrons:
        """The collection of all :class:`GeneIntrons` for each gene"""
        return GeneIntrons(self._sg._introns)

    @property
    def junctions(self) -> GeneJunctions:
        """The collection of all :class:`GeneJunctions` for each gene"""
        return GeneJunctions(self._sg._junctions)

    @property
    def exon_connections(self) -> ExonConnections:
        """:class:`ExonConnections` for the splicegraph"""
        return ExonConnections(self._sg._exon_connections)

    def to_zarr(self, path: Union[str, Path]) -> None:
        """Save :py:class:`SpliceGraph` to specified path

        Parameters
        ----------
        path: Union[str, Path]
            Path for zarr store with SpliceGraph components
        """
        if Path(path).exists():
            raise ValueError(
                f"Will not save splicegraph to already existing file {path}."
                " Please delete and try again if you want it there, otherwise"
                " pick a different output path"
            )
        self.contigs.to_zarr(path, "w", consolidated=False)
        self.genes.to_zarr(path, "a", consolidated=False)
        self.exons.to_zarr(path, "a", consolidated=False)
        self.introns.to_zarr(path, "a", consolidated=False)
        self.junctions.to_zarr(path, "a", consolidated=True)
        return

    @classmethod
    def from_zarr(
        cls, path: Union[str, Path], genes: Optional[Genes] = None
    ) -> "SpliceGraph":
        """Load :py:class:`SpliceGraph` from specified path

        Parameters
        ----------
        path: Union[str, Path]
            Path where splicegraph is stored in zarr format
        genes: Optional[Genes]
            If specified, :py:class:`Genes` that has already been loaded.
            Used when multiple objects refer to the same set of genes.
            Otherwise, load from path.
        """
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

    def to_sqlite(
        self, path: Union[str, Path], genome_name: str = "default_genome"
    ) -> None:
        """Save splicegraph to legacy format"""
        import sqlite3

        import numpy as np
        import pandas as pd

        # open up connection
        conn = sqlite3.connect(path)

        # we will want to index into gene_ids
        gene_id = np.array(self.genes.gene_id)

        # save table of genes
        pd.DataFrame(
            {
                "name": self.genes.gene_name,
                "strand": [x.decode() for x in self.genes.strand],
                "chromosome": np.array(self.contigs.seqid)[self.genes.contig_idx],
            },
            index=pd.Index(gene_id, name="id"),
        ).to_sql(
            "gene",
            conn,
            dtype={
                "id": "VARCHAR NOT NULL",
                "name": "VARCHAR",
                "strand": "VARCHAR",
                "chromosome": "VARCHAR",
            },
            schema=None,
        )
        # save table of exons
        pd.DataFrame(
            {
                "gene_id": gene_id[self.exons.gene_idx],
                "start": self.exons.start,
                "end": self.exons.end,
                "annotated_start": np.where(
                    self.exons.is_denovo(), self.exons.start, self.exons.annotated_start
                ),
                "annotated_end": np.where(
                    self.exons.is_denovo(), self.exons.end, self.exons.annotated_end
                ),
                "annotated": ~self.exons.is_denovo(),
            }
        ).set_index(["gene_id", "start", "end"]).to_sql(
            "exon",
            conn,
            dtype={
                "gene_id": "VARCHAR NOT NULL",
                "start": "INTEGER NOT NULL",
                "end": "INTEGER NOT NULL",
                "annotated_start": "INTEGER",
                "annotated_end": "INTEGER",
                "annotated": "BOOLEAN",
            },
        )
        # save table of junctions
        pd.DataFrame(
            {
                "gene_id": gene_id[self.junctions.gene_idx],
                "start": self.junctions.start,
                "end": self.junctions.end,
                # we don't track has_reads, so not exact
                "has_reads": self.junctions.passed_build,
                "annotated": ~self.junctions.denovo,
                "is_simplified": self.junctions.simplified,
                "is_constitutive": self.exon_connections.is_constitutive(
                    self.junctions.src_exon_idx(), np.array(b"s")
                ),
                "has_flag": self.junctions.passed_build,
            }
        ).set_index(["gene_id", "start", "end"]).to_sql(
            "junction",
            conn,
            dtype={
                "gene_id": "VARCHAR NOT NULL",
                "start": "INTEGER NOT NULL",
                "end": "INTEGER NOT NULL",
                "has_reads": "BOOLEAN",
                "annotated": "BOOLEAN",
                "is_simplified": "BOOLEAN",
                "is_constitutive": "BOOLEAN",
                "has_flag": "BOOLEAN",
            },
            schema=None,
        )
        # save table of introns
        pd.DataFrame(
            {
                "gene_id": gene_id[self.introns.gene_idx],
                "start": self.introns.start,
                "end": self.introns.end,
                # we don't track has_reads, so not exact
                "has_reads": self.introns.passed_build,
                "annotated": ~self.introns.denovo,
                "is_simplified": self.introns.simplified,
                "is_constitutive": self.exon_connections.is_constitutive(
                    self.introns.src_exon_idx(), np.array(b"s")
                ),
                "has_flag": self.introns.passed_build,
            }
        ).set_index(["gene_id", "start", "end"]).to_sql(
            "intron_retention",
            conn,
            dtype={
                "gene_id": "VARCHAR NOT NULL",
                "start": "INTEGER NOT NULL",
                "end": "INTEGER NOT NULL",
                "has_reads": "BOOLEAN",
                "annotated": "BOOLEAN",
                "is_simplified": "BOOLEAN",
                "is_constitutive": "BOOLEAN",
                "has_flag": "BOOLEAN",
            },
            schema=None,
        )
        # save empty tables of alt_starts, alt_ends, gene overlap, file version,
        for alt_table in ("alt_start", "alt_end"):
            pd.DataFrame({"gene_id": [], "coordinate": []}).astype(
                {"gene_id": str, "coordinate": int}
            ).set_index(["gene_id", "coordinate"]).to_sql(
                alt_table,
                conn,
                dtype={"gene_id": "VARCHAR NOT NULL", "coordinate": "INTEGER NOT NULL"},
            )
        pd.DataFrame().reindex(columns=["gene_id_1", "gene_id_2"]).astype(
            str
        ).set_index(["gene_id_1", "gene_id_2"]).to_sql(
            "gene_overlap",
            conn,
            dtype={"gene_id_1": "VARCHAR NOT NULL", "gene_id_2": "VARCHAR NOT NULL"},
        )
        pd.DataFrame().reindex(columns=["id", "value"]).astype(int).set_index(
            "id"
        ).to_sql(
            "file_version", conn, dtype={"id": "INTEGER NOT NULL", "value": "INTEGER"}
        )
        # genome name
        pd.DataFrame(
            {"name": [genome_name]}, index=pd.Index([1], name="id", dtype=int)
        ).to_sql("genome", conn, dtype={"id": "INTEGER NOT NULL", "name": "VARCHAR"})
        # empty experiment information
        pd.DataFrame({}, index=pd.Index([], name="name", dtype=str)).to_sql(
            "experiment", conn, dtype={"name": "VARCHAR NOT NULL"}
        )
        pd.DataFrame().reindex(
            columns=[
                "experiment_name",
                "junction_gene_id",
                "junction_start",
                "junction_end",
                "reads",
            ]
        ).astype(
            {
                "experiment_name": str,
                "junction_gene_id": str,
                "junction_start": int,
                "junction_end": int,
                "reads": int,
            }
        ).set_index(
            ["experiment_name", "junction_gene_id", "junction_start", "junction_end"]
        ).to_sql(
            "junction_reads",
            conn,
            dtype={
                "experiment_name": "VARCHAR NOT NULL",
                "junction_gene_id": "VARCHAR NOT NULL",
                "junction_start": "INTEGER NOT NULL",
                "junction_end": "INTEGER NOT NULL",
                "reads": "INTEGER NOT NULL",
            },
        )
        pd.DataFrame().reindex(
            columns=[
                "experiment_name",
                "intron_retention_gene_id",
                "intron_retention_start",
                "intron_retention_end",
                "reads",
            ]
        ).astype(
            {
                "experiment_name": str,
                "intron_retention_gene_id": str,
                "intron_retention_start": int,
                "intron_retention_end": int,
                "reads": int,
            }
        ).set_index(
            [
                "experiment_name",
                "intron_retention_gene_id",
                "intron_retention_start",
                "intron_retention_end",
            ]
        ).to_sql(
            "intron_retention_reads",
            conn,
            dtype={
                "experiment_name": "VARCHAR NOT NULL",
                "intron_retention_gene_id": "VARCHAR NOT NULL",
                "intron_retention_start": "INTEGER NOT NULL",
                "intron_retention_end": "INTEGER NOT NULL",
                "reads": "INTEGER NOT NULL",
            },
        )

        # close connection
        conn.close()
        return
