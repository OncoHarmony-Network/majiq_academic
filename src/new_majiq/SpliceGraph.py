"""
SpliceGraph.py

Combined regions of contigs containing genes containing exons connected by
introns and junctions

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Dict, Final, Optional, Union

import new_majiq.constants as constants
from new_majiq.Contigs import Contigs
from new_majiq.ExonConnections import ExonConnections
from new_majiq.Exons import Exons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.Genes import Genes
from new_majiq.internals import GFF3FeatureType
from new_majiq.internals import SpliceGraph as _SpliceGraph


class SpliceGraph(object):
    """Representation of all genes as exons connected by introns and junctions

    Parameters
    ----------
    sg: _SpliceGraph
        Underlying object binding the internal C++ API

    See Also
    --------
    SpliceGraph.from_gff3
    SpliceGraph.from_components
    SpliceGraph.from_zarr
    """

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
        gff3_types: Dict[str, GFF3FeatureType] = constants.DEFAULT_BUILD_GFF3TYPES,
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
