"""
splicegraph.py (new_majiq.mock)

Create mocked splicegraphs

Author: Joseph K Aicher
"""

from typing import List, NamedTuple, Tuple

import numpy as np

from new_majiq.core.SpliceGraph import (
    Contigs,
    Exons,
    GeneIntrons,
    GeneJunctions,
    Genes,
    SpliceGraph,
)

DEFAULT_EXON_LENGTH = 100
DEFAULT_INTRON_LENGTH = 300
DEFAULT_EXTENSION_LENGTH = 24


class _SpliceGraphCoordinates(NamedTuple):
    exons: List[Tuple[int, int]]
    introns: List[Tuple[int, int]]
    junctions: List[Tuple[int, int]]


def mock_splicegraph(
    num_modules: int,
    start: int = 1,
    exon_length: int = DEFAULT_EXON_LENGTH,
    intron_length: int = DEFAULT_INTRON_LENGTH,
    extension_length: int = DEFAULT_EXTENSION_LENGTH,
) -> SpliceGraph:
    """create mock splicegraph"""
    return _coordinates_to_splicegraph(
        _mock_modules(
            num_modules,
            start=start,
            exon_length=exon_length,
            intron_length=intron_length,
            extension_length=extension_length,
        )
    )


def _coordinates_to_splicegraph(x: _SpliceGraphCoordinates) -> SpliceGraph:
    """convert coordinates for single gene to SpliceGraph"""
    contigs = Contigs.from_list(["mock_contig"])
    genes = Genes.from_arrays(
        contigs, [0], [x.exons[0][0]], [x.exons[-1][-1]], ["+"], ["mock_gene"]
    )
    exons = Exons.from_arrays(
        genes, [0] * len(x.exons), [y for y, _ in x.exons], [y for _, y in x.exons]
    )
    introns = GeneIntrons.from_arrays(
        genes,
        [0] * len(x.introns),
        [y for y, _ in x.introns],
        [y for _, y in x.introns],
        passed_build=[True] * len(x.introns),
    )
    junctions = GeneJunctions.from_arrays(
        genes,
        [0] * len(x.junctions),
        [y for y, _ in x.junctions],
        [y for _, y in x.junctions],
        passed_build=[True] * len(x.junctions),
    )
    return SpliceGraph.from_components(contigs, genes, exons, junctions, introns)


def _mock_modules(
    num_modules: int,
    start: int = 1,
    exon_length: int = DEFAULT_EXON_LENGTH,
    intron_length: int = DEFAULT_INTRON_LENGTH,
    extension_length: int = DEFAULT_EXTENSION_LENGTH,
) -> _SpliceGraphCoordinates:
    """Get coordinates for an multiple modules starting at start"""
    assert num_modules > 0
    assert start > 0
    assert exon_length > 2 * extension_length
    assert intron_length > 0
    assert extension_length > 0
    result = _SpliceGraphCoordinates([], [], [])
    module_events = np.random.randint(2, size=(num_modules, 4), dtype=bool)
    for i, events_i in enumerate(module_events):
        es, alt_start, alt_end, ir = events_i
        module_i = _mock_single_module(
            es,
            alt_start,
            alt_end,
            ir,
            start=start,
            exon_length=exon_length,
            intron_length=intron_length,
            extension_length=extension_length,
        )
        result.exons.extend(module_i.exons)
        result.introns.extend(module_i.introns)
        result.junctions.extend(module_i.junctions)
        if i < num_modules - 1:  # there are still modules left
            start = result.exons[-1][-1] + intron_length + 1
            result.junctions.append((result.exons[-1][-1], start))
    return result


def _mock_single_module(
    es: bool,
    alt_start: bool,
    alt_end: bool,
    ir: bool,
    start: int = 1,
    exon_length: int = DEFAULT_EXON_LENGTH,
    intron_length: int = DEFAULT_INTRON_LENGTH,
    extension_length: int = DEFAULT_EXTENSION_LENGTH,
) -> _SpliceGraphCoordinates:
    """Get coordinates for an individual module starting at start"""
    # assert exon_length > 2 * extension_length
    exon_idx = np.arange(3 if es else 2)
    exon_start = start + (exon_length + intron_length) * exon_idx
    exon_end = (start + exon_length - 1) + (exon_length + intron_length) * exon_idx
    intron_start = 1 + exon_end[:-1]
    intron_end = exon_start[1:] - 1
    junction_start = list(intron_start - 1)
    junction_end = list(1 + intron_end)
    if es:
        junction_start.append(exon_end[0])
        junction_end.append(exon_start[-1])
    if alt_start:
        junction_start.append(exon_end[0] - extension_length)
        junction_end.append(np.random.choice(exon_start[1:]))
    if alt_end:
        junction_start.append(np.random.choice(exon_end[:-1]))
        junction_end.append(exon_start[-1] + extension_length)
    if (not ir) and (es or alt_start or alt_end):
        intron_start = intron_start[[]]
        intron_end = intron_end[[]]
    return _SpliceGraphCoordinates(
        exons=sorted(zip(exon_start, exon_end)),
        introns=sorted(zip(intron_start, intron_end)),
        junctions=sorted(zip(junction_start, junction_end)),
    )
