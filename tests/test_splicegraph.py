"""
test_splicegraph.py

Unit tests of expected functionality with respect to splicegraphs

Author: Joseph K Aicher
"""

from pathlib import Path

import pytest
import xarray as xr

import new_majiq as nm


@pytest.fixture(scope="module")
def base_contigs() -> nm.Contigs:
    """contigs C, Ax, B"""
    return nm.Contigs.from_list(["C", "Ax", "B"])


@pytest.fixture(scope="module")
def base_genes(base_contigs: nm.Contigs) -> nm.Genes:
    return nm.Genes.from_arrays(
        base_contigs,
        # each contig has 2 genes: nonoverlapping, overlapping opposite
        # strands, overlapping same strand
        [0, 0, 1, 1, 2, 2],
        [100, 10100, 100, 10100, 100, 10100],
        [10000, 20000, 20000, 20000, 20000, 20000],
        ["+", "+", "+", "-", "+", "+"],
        # *no* overlap, *opp*osite strands, *same* strands
        ["no1", "no2", "opp1", "opp2", "same1", "same2"],
        ["no1", "no2", "opp1", "opp2", "same1", "same2"],
    )


def test_contigs_nonunique():
    with pytest.raises(ValueError):
        contigs = nm.Contigs.from_list(["A", "A"])  # noqa: F841
    return


def test_contigs_roundtrip(base_contigs: nm.Contigs, tmpdir: Path):
    path = tmpdir / "contigs.zarr"
    base_contigs.to_zarr(path, mode="w")
    roundtrip_contigs = nm.Contigs.from_zarr(path)
    xr.testing.assert_equal(base_contigs.df, roundtrip_contigs.df)
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


def test_genes_roundtrip(base_genes: nm.Genes, tmpdir: Path):
    path = tmpdir / "genes.zarr"
    base_genes.contigs.to_zarr(path, mode="w", consolidated=False)
    base_genes.to_zarr(path, mode="a", consolidated=True)
    roundtrip_genes = nm.Genes.from_zarr(path)
    xr.testing.assert_equal(base_genes.df, roundtrip_genes.df)
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
