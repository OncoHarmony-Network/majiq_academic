"""
StrandDetection.py

Approach for automatically detecting and adjusting strandness of SJJunctionsBins

Author: Joseph K Aicher
"""

import numpy as np
import new_majiq.constants as constants

from new_majiq.logger import get_logger

from new_majiq.SpliceGraph import SpliceGraph
from new_majiq.SpliceGraphReads import SpliceGraphReads
from new_majiq.SJJunctionsBins import SJJunctionsBins
from new_majiq.SJIntronsBins import SJIntronsBins
from new_majiq.SJIntrons import SJIntrons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.internals import ExperimentStrandness


from typing import (
    NamedTuple,
)


class DetectedSJStrand(NamedTuple):
    updated_sjbins: SJJunctionsBins
    updated_strandness: ExperimentStrandness


def detect_strand(
    sjbins: SJJunctionsBins,
    strandness: ExperimentStrandness,
    sg: SpliceGraph,
    minreads: int = constants.DEFAULT_BAM_STRAND_MINREADS,
    minjunctions: int = constants.DEFAULT_BAM_STRAND_MINJUNCTIONS,
    mindeviation: float = constants.DEFAULT_BAM_STRAND_MINDEVIATION,
) -> DetectedSJStrand:
    """Use splicegraph to get updated sjbins, strandness"""
    log = get_logger()
    if strandness == ExperimentStrandness.NONE:
        # there's nothing to do
        return DetectedSJStrand(
            updated_sjbins=sjbins,
            updated_strandness=strandness,
        )
    log.info("Detecting strandedness of SJJunctionsBins with SpliceGraph")
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
        return DetectedSJStrand(sjbins.to_unstranded(), ExperimentStrandness.NONE)
    # compute ratio
    median_ratio = np.median(ratios)
    if np.abs(median_ratio - 0.5) < mindeviation:
        log.info(
            f"Median ratio of original stranded reads vs total is {median_ratio},"
            f" which is within {mindeviation} of 0.5"
        )
        return DetectedSJStrand(sjbins.to_unstranded(), ExperimentStrandness.NONE)
    else:
        log.info(
            f"Median ratio of original stranded reads vs total is {median_ratio},"
            f" which is more than {mindeviation} from 0.5"
        )
    # at this point, we know it's stranded, so we flip if median_ratio < 0.5
    if median_ratio < 0.5:
        # we need to flip it
        return DetectedSJStrand(
            updated_sjbins=sjbins_flipped,
            updated_strandness=(
                ExperimentStrandness.REVERSE
                if strandness == ExperimentStrandness.FORWARD
                else ExperimentStrandness.FORWARD
            ),
        )
    else:
        # no flipping it
        return DetectedSJStrand(sjbins, strandness)
