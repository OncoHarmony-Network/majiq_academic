"""
StrandDetection.py

Approach for automatically detecting and adjusting strandness of SJJunctionsBins

Author: Joseph K Aicher
"""

import numpy as np

import new_majiq.constants as constants
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.internals import ExperimentStrandness
from new_majiq.logger import get_logger
from new_majiq.SJIntrons import SJIntrons
from new_majiq.SJIntronsBins import SJIntronsBins
from new_majiq.SJJunctionsBins import SJJunctionsBins
from new_majiq.SpliceGraph import SpliceGraph
from new_majiq.SpliceGraphReads import SpliceGraphReads


def detect_strand(
    sjbins: SJJunctionsBins,
    sg: SpliceGraph,
    minreads: int = constants.DEFAULT_BAM_STRAND_MINREADS,
    minjunctions: int = constants.DEFAULT_BAM_STRAND_MINJUNCTIONS,
    mindeviation: float = constants.DEFAULT_BAM_STRAND_MINDEVIATION,
) -> SJJunctionsBins:
    """Use splicegraph to get updated sjbins, strandness"""
    log = get_logger()
    if sjbins.strandness == ExperimentStrandness.NONE:
        # there's nothing to do
        return sjbins
    # otherwise, we are going to detect strand by matching to splicegraph
    empty_sg_introns = GeneIntrons.from_genes(sg.genes)
    empty_sibins = SJIntronsBins.from_regions(
        SJIntrons.from_contigs(sjbins.regions.contigs),
        sjbins.total_bins,
        original_path=sjbins.original_path,
        original_version=sjbins.original_version,
        original_time=sjbins.original_time,
    )
    # get reads in original and reversed strand directions
    reads = SpliceGraphReads.from_connections_and_sj(
        empty_sg_introns, sg.junctions, empty_sibins, sjbins
    ).junctions_reads
    sjbins_flipped = sjbins.flip_strand()
    reads_flipped = SpliceGraphReads.from_connections_and_sj(
        empty_sg_introns, sg.junctions, empty_sibins, sjbins_flipped
    ).junctions_reads
    # compute ratios for junctions with at least minreads
    reads_total = reads + reads_flipped
    minreads_mask = reads_total >= minreads
    ratios = reads[minreads_mask] / reads_total[minreads_mask]
    # check that we had enough junctions
    if len(ratios) < minjunctions:
        log.info(
            f"Only {len(ratios)} junctions with at least {minreads} reads"
            f" vs minimum {minjunctions} junctions for inferring strandedness"
        )
        return sjbins.to_unstranded()
    # compute ratio
    median_ratio = np.median(ratios)
    deviation = np.abs(median_ratio - 0.5)
    log.info(
        f"Median ratio of original stranded reads vs total is {median_ratio:.1%}"
        f" (deviates by {deviation:.1%} from unstranded expectation)"
        f" from {len(ratios)} junctions with at least {minreads} reads"
    )
    if deviation < mindeviation:
        # not far enough from 0.5 to justify strandedness
        return sjbins.to_unstranded()
    elif median_ratio < 0.5:
        # we need to flip it
        return sjbins_flipped
    else:
        # current strandedness is appropriate
        return sjbins
