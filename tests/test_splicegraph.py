"""
test_splicegraph.py

Unit tests of expected functionality with respect to splicegraphs

Author: Joseph K Aicher
"""

from pathlib import Path

import pytest
import xarray as xr

import new_majiq as nm


@pytest.fixture
def base_contigs() -> nm.Contigs:
    """contigs C, Ax, B"""
    return nm.Contigs.from_list(["C", "Ax", "B"])


@pytest.fixture
def base_genes(base_contigs: nm.Contigs) -> nm.Genes:
    return nm.Genes.from_arrays(
        base_contigs,
        # each contig has 2 genes: nonoverlapping, overlapping opposite
        # strands, overlapping same strand
        [0, 0, 1, 1, 2, 2],
        [100, 5000, 100, 5000, 100, 5000],
        [4500, 10000, 10000, 10000, 20000, 10000],
        ["+", "+", "+", "-", "+", "+"],
        # *no* overlap, *opp*osite strands, *same* strands
        ["no1", "no2", "opp1", "opp2", "same1", "same2"],
        ["no1", "no2", "opp1", "opp2", "same1", "same2"],
    )


def get_starts(coords):
    """get starts from list[(start, end)]"""
    return [x[0] for x in coords]


def get_ends(coords):
    """get ends from list[(start, end)]"""
    return [x[1] for x in coords]


@pytest.fixture
def base_genejunctions(base_genes: nm.Genes) -> nm.GeneJunctions:
    # one shared junction (6200, 8000) (ignoring strand)
    even_junctions_short = sorted(
        [
            (300, 600),
            (300, 4000),
            (800, 4000),
        ]
    )
    even_junctions_long = sorted(
        [
            *even_junctions_short,
            (4500, 6000),
            (6200, 8000),
            (8100, 8300),
            (8400, 9800),
        ]
    )
    odd_junctions = sorted(
        [
            (5500, 6000),
            (6200, 8000),
            (8100, 8600),
            (8100, 8650),
            (8800, 9800),
        ]
    )
    # construct indexes/coordinates for junctions
    gene_idx = (
        [0] * len(even_junctions_short)
        + [1] * len(odd_junctions)
        + [2] * len(even_junctions_long)
        + [3] * len(odd_junctions)
        + [4] * len(even_junctions_long)
        + [5] * len(odd_junctions)
    )
    start = (
        get_starts(even_junctions_short)
        + get_starts(odd_junctions)
        + get_starts(even_junctions_long)
        + get_starts(odd_junctions)
        + get_starts(even_junctions_long)
        + get_starts(odd_junctions)
    )
    end = (
        get_ends(even_junctions_short)
        + get_ends(odd_junctions)
        + get_ends(even_junctions_long)
        + get_ends(odd_junctions)
        + get_ends(even_junctions_long)
        + get_ends(odd_junctions)
    )
    return nm.GeneJunctions.from_arrays(base_genes, gene_idx, start, end)


@pytest.fixture
def base_geneintrons(base_genes: nm.Genes) -> nm.GeneIntrons:
    return nm.GeneIntrons.from_arrays(base_genes, [0, 2, 4], [301] * 3, [599] * 3)


@pytest.fixture
def base_exons(base_genes: nm.Genes) -> nm.Exons:
    # 0, 2, 4 genes are designed to share same region 100, 4500 with alternatives
    # 2 and 4 have additional constitutive annotated exons
    even_exons_short = sorted(
        [
            (100, 300),
            # there should be annotated intron [301, 599]
            (600, 800),
            (4000, 4500),
        ]
    )
    shared_exons = [
        (6000, 6200),
        (8000, 8100),
        (9800, 10000),
    ]
    odd_exons = sorted([(5000, 5500), (8600, 8800), *shared_exons])
    even_exons_long = sorted([*even_exons_short, *shared_exons, (8300, 8400)])

    # construct indexes/coordinates for exons
    gene_idx = (
        [0] * len(even_exons_short)
        + [1] * len(odd_exons)
        + [2] * len(even_exons_long)
        + [3] * len(odd_exons)
        + [4] * len(even_exons_long)
        + [5] * len(odd_exons)
    )
    start = (
        get_starts(even_exons_short)
        + get_starts(odd_exons)
        + get_starts(even_exons_long)
        + get_starts(odd_exons)
        + get_starts(even_exons_long)
        + get_starts(odd_exons)
    )
    end = (
        get_ends(even_exons_short)
        + get_ends(odd_exons)
        + get_ends(even_exons_long)
        + get_ends(odd_exons)
        + get_ends(even_exons_long)
        + get_ends(odd_exons)
    )
    # construct exons
    return nm.Exons.from_arrays(base_genes, gene_idx, start, end)


@pytest.fixture
def base_splicegraph(
    base_contigs: nm.Contigs,
    base_genes: nm.Genes,
    base_exons: nm.Exons,
    base_geneintrons: nm.GeneIntrons,
    base_genejunctions: nm.GeneJunctions,
) -> nm.SpliceGraph:
    return nm.SpliceGraph.from_components(
        base_contigs, base_genes, base_exons, base_genejunctions, base_geneintrons
    )


def test_contigs_nonunique():
    with pytest.raises(ValueError):
        contigs = nm.Contigs.from_list(["A", "A"])  # noqa: F841
    return


def test_contigs_attributes(base_contigs: nm.Contigs):
    base_contigs.checksum()  # checksum should work
    assert base_contigs["C"] == 0
    assert base_contigs["Ax"] == 1
    assert base_contigs["B"] == 2
    with pytest.raises(KeyError):
        base_contigs["invalid contig"]
    base_contigs.df
    return


def test_genes_constructor(base_contigs: nm.Contigs):
    # this should work
    nm.Genes.from_arrays(base_contigs, [0], [1], [3], ["+"], ["id"], ["name"])
    # nonunique id constraint
    with pytest.raises(ValueError):
        nm.Genes.from_arrays(
            base_contigs, [0] * 2, [1] * 2, [3] * 2, ["+"] * 2, ["id"] * 2, ["name"] * 2
        )
    # nonunique names are allowed, and genes can have exactly the same coordinates
    nm.Genes.from_arrays(
        base_contigs, [0] * 2, [1] * 2, [3] * 2, ["+"] * 2, ["id1", "id2"], ["name"] * 2
    )
    # but regions have to be ordered
    with pytest.raises(ValueError):
        nm.Genes.from_arrays(
            base_contigs,
            [1, 0],
            [1] * 2,
            [3] * 2,
            ["+"] * 2,
            ["id1", "id2"],
            ["name"] * 2,
        )
    with pytest.raises(ValueError):
        nm.Genes.from_arrays(
            base_contigs,
            [0] * 2,
            [1, 0],
            [3] * 2,
            ["+"] * 2,
            ["id1", "id2"],
            ["name"] * 2,
        )
    with pytest.raises(ValueError):
        nm.Genes.from_arrays(
            base_contigs,
            [0] * 2,
            [1] * 2,
            [3, 2],
            ["+"] * 2,
            ["id1", "id2"],
            ["name"] * 2,
        )
    with pytest.raises(ValueError):
        nm.Genes.from_arrays(
            base_contigs,
            [0] * 2,
            [1] * 2,
            [3] * 2,
            ["-", "+"],
            ["id1", "id2"],
            ["name"] * 2,
        )
    return


def test_genes_attributes(base_genes: nm.Genes):
    base_genes.checksum()  # checksum should work
    assert base_genes["no1"] == 0
    assert base_genes["opp2"] == 3
    assert base_genes["same1"] == 4
    assert base_genes.slice_for_contig(0) == slice(0, 2)
    assert base_genes.slice_for_contig(2) == slice(4, 6)
    with pytest.raises(KeyError):
        base_genes["invalid gene"]
    base_genes.df
    return


def test_splicegraph_roundtrip(base_splicegraph: nm.SpliceGraph, tmpdir: Path):
    path = tmpdir / "sg.zarr"
    base_splicegraph.to_zarr(path)
    roundtrip_splicegraph = nm.SpliceGraph.from_zarr(path)
    xr.testing.assert_equal(
        base_splicegraph.contigs.df, roundtrip_splicegraph.contigs.df
    )
    xr.testing.assert_equal(base_splicegraph.genes.df, roundtrip_splicegraph.genes.df)
    xr.testing.assert_equal(base_splicegraph.exons.df, roundtrip_splicegraph.exons.df)
    xr.testing.assert_equal(
        base_splicegraph.introns.df, roundtrip_splicegraph.introns.df
    )
    xr.testing.assert_equal(
        base_splicegraph.junctions.df, roundtrip_splicegraph.junctions.df
    )
