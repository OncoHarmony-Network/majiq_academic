"""
SJExperiment.py

Convenience class for representing SJIntronsBins and SJJunctionsBins for the
same experiment

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Final, Union

import numpy as np

import new_majiq.constants as constants
from new_majiq.Exons import Exons
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.internals import ExperimentStrandness
from new_majiq.logger import get_logger
from new_majiq.SJIntrons import SJIntrons
from new_majiq.SJIntronsBins import SJIntronsBins
from new_majiq.SJJunctionsBins import SJJunctionsBins
from new_majiq.SpliceGraph import SpliceGraph
from new_majiq.SpliceGraphReads import SpliceGraphReads


class SJExperiment(object):
    """Per-bin read coverage over introns and junctions in the same experiment

    Parameters
    ----------
    introns: SJIntronsBins
    junctions: SJJunctionsBins

    Notes
    -----
    Enforces requirement that introns and junctions share the same original
    (bam) path and version of majiq used to do parsing
    """

    def __init__(self, introns: SJIntronsBins, junctions: SJJunctionsBins):
        """Per bin reads for introns and junctions for the same experiment

        Parameters
        ----------
        introns: SJIntronsBins
        junctions: SJJunctionsBins
        """
        if introns.original_path != junctions.original_path:
            raise ValueError(
                "SJExperiment introns and junctions must have same original path"
            )
        if introns.original_version != junctions.original_version:
            raise ValueError(
                "SJExperiment introns and junctions must have same original version"
            )
        self._introns: Final[SJIntronsBins] = introns
        self._junctions: Final[SJJunctionsBins] = junctions
        return

    @property
    def introns(self) -> SJIntronsBins:
        """Bin reads for introns for this experiment"""
        return self._introns

    @property
    def junctions(self) -> SJJunctionsBins:
        """Bin reads for junctions for this experiment"""
        return self._junctions

    @property
    def original_path(self) -> str:
        return self.junctions.original_path

    @property
    def original_version(self) -> str:
        return self.junctions.original_version

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "SJExperiment":
        """Load introns and junctions reads from specified path"""
        return SJExperiment(
            SJIntronsBins.from_zarr(path), SJJunctionsBins.from_zarr(path)
        )

    def to_zarr(self, path: Union[str, Path], consolidated: bool = True) -> None:
        """Serialize to zarr format"""
        self.junctions.to_zarr(path, consolidated=False)
        self.introns.to_zarr(path, consolidated=True, check_experiment_if_exists=False)
        return

    @classmethod
    def from_bam(
        cls,
        path: Union[str, Path],
        sg: SpliceGraph,
        strandness: Union[str, ExperimentStrandness] = "auto",
        update_exons: bool = False,
        nthreads: int = constants.DEFAULT_BAM_NTHREADS,
        allow_disjoint_contigs: bool = constants.DEFAULT_BAM_ALLOW_DISJOINT_CONTIGS,
        auto_minreads: int = constants.DEFAULT_BAM_STRAND_MINREADS,
        auto_minjunctions: int = constants.DEFAULT_BAM_STRAND_MINJUNCTIONS,
        auto_mediantolerance: float = constants.DEFAULT_BAM_STRAND_MINDEVIATION,
    ) -> "SJExperiment":
        """Load SJ information from path using given splicegraph, strandness

        Parameters
        ----------
        path: Union[str, Path]
            Path for input BAM file
        sg: SpliceGraph
            Used to determine regions for contig intron coverage
        strandness: Union[str, ExperimentStrandness]
            strandness to parse BAM files with. If str, one of "auto",
            "forward", "reverse", or "none" (not case-sensitive).
            If auto, automatically detects strand using median ratio of forward
            vs reverse stranded reads at annotated junctions
        update_exons: bool
            Experimental -- use junction coverage to definitively ignore
            intronic coverage in potential denovo exons (or exon extension)
        nthreads: int
            Number of threads to parse BAM with
        allow_disjoint_contigs: bool
            If true, don't raise error if BAM doesn't have any contigs that
            overlap sg
        auto_minreads: int
            Only consider evidence from splicegraph junctions with at least
            this many total (unstranded) reads
        auto_minjunctions: int
            Infer unstranded if the number of splicegraph junctions with
            sufficient reads is less than this argument
        auto_mediantolerance: float
            Infer unstranded if the median proportion of junctions of forward
            strand vs all reads deviates from 0.5 by at most this amount
        """
        log = get_logger()
        # manage strand
        autostrand: bool = False
        if not isinstance(strandness, ExperimentStrandness):
            strandness = strandness.upper()  # make uppercase
            if strandness == "AUTO":
                autostrand = True
                strandness = ExperimentStrandness.FORWARD
            elif strandness == "FORWARD":
                strandness = ExperimentStrandness.FORWARD
            elif strandness == "REVERSE":
                strandness = ExperimentStrandness.REVERSE
            elif strandness == "NONE":
                strandness = ExperimentStrandness.NONE
            else:
                raise ValueError(
                    f"Invalid strandness {strandness},"
                    " must be AUTO, FORWARD, REVERSE, or NONE"
                )
        # load junctions
        log.info(f"Parsing alignments from {path} for junctions")
        junctions: SJJunctionsBins = SJJunctionsBins.from_bam(
            path, strandness=strandness, nthreads=nthreads
        )
        if not (set(sg.contigs.seqid) & set(junctions.regions.contigs.seqid)):
            if allow_disjoint_contigs:
                log.warning("Contigs from splicegraph and BAM are disjoint!")
            else:
                log.error(
                    "Contigs from splicegraph and BAM are disjoint!\n"
                    f" sg contigs = {sg.contigs.seqid},\n"
                    f" bam contigs = {junctions.regions.contigs.seqid},\n"
                    "Set allow_disjoint_contigs if this is okay."
                )
                raise RuntimeError("Contigs from splicegraph and BAM are disjoint")
        if autostrand:
            log.info("Inferring strandness comparing counts on splicegraph junctions")
            junctions = cls.detect_strand(
                junctions,
                sg,
                minreads=auto_minreads,
                minjunctions=auto_minjunctions,
                mindeviation=auto_mediantolerance,
            )
            log.info(f"Inferred strandness: {junctions.strandness.name}")
        # determine introns
        log.info("Using gene introns/exons to define regions for intronic coverage")
        gene_introns: GeneIntrons = sg.introns
        exons: Exons = sg.exons
        if update_exons:
            log.info("Identifying potential denovo exons from input junctions")
            # TODO (change parameters used for reliable updated junctions?)
            updated_junctions = (
                sg.junctions.builder()
                .add_group(sg.junctions.build_group(exons).add_experiment(junctions))
                .get_passed()
            )
            exons = exons.infer_with_junctions(updated_junctions)
        log.info(f"Parsing alignments from {path} for introns")
        introns = SJIntronsBins.from_bam(
            path,
            junctions.total_bins,
            exons,
            gene_introns,
            strandness=junctions.strandness,
            nthreads=nthreads,
        )
        return SJExperiment(introns, junctions)

    @staticmethod
    def detect_strand(
        sj_junctions: SJJunctionsBins,
        sg: SpliceGraph,
        minreads: int = constants.DEFAULT_BAM_STRAND_MINREADS,
        minjunctions: int = constants.DEFAULT_BAM_STRAND_MINJUNCTIONS,
        mindeviation: float = constants.DEFAULT_BAM_STRAND_MINDEVIATION,
    ) -> SJJunctionsBins:
        """Use splicegraph to get correct strand orientation of junctions"""
        log = get_logger()
        if sj_junctions.strandness == ExperimentStrandness.NONE:
            # there's nothing to do
            return sj_junctions
        # otherwise, we are going to detect strand by matching to splicegraph

        # ignore introns
        empty_sg_introns = GeneIntrons.from_genes(sg.genes)
        sj = SJExperiment.from_SJJunctionsBins(sj_junctions)

        # get reads in original and reversed strand directions
        reads = SpliceGraphReads._internals_from_connections_and_sj(
            empty_sg_introns, sg.junctions, sj
        ).junctions_reads
        sj_flipped = SJExperiment(sj.introns, sj_junctions.flip_strand())
        reads_flipped = SpliceGraphReads._internals_from_connections_and_sj(
            empty_sg_introns, sg.junctions, sj_flipped
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
            return sj_junctions.to_unstranded()

        # compute ratio
        median_ratio = np.median(ratios)
        deviation = np.abs(median_ratio - 0.5)
        log.info(
            f"Median ratio of original stranded reads vs total is {median_ratio:.1%}"
            f" (deviates by {deviation:.1%} from unstranded expectation)"
            f" from {len(ratios)} junctions with at least {minreads} reads"
        )

        # return sj junctions indicated by deviation from 50%
        if deviation < mindeviation:
            # not far enough from 0.5 to justify strandedness
            return sj_junctions.to_unstranded()
        elif median_ratio < 0.5:
            # we need to flip it
            return sj_flipped.junctions
        else:
            # current strandedness is appropriate
            return sj_junctions

    @classmethod
    def from_SJJunctionsBins(cls, junctions: SJJunctionsBins) -> "SJExperiment":
        """SJExperiment with given junction coverage, zero intron coverage"""
        return SJExperiment(
            SJIntronsBins.from_regions(
                SJIntrons.from_contigs(junctions.regions.contigs),
                junctions.total_bins,
                original_path=junctions.original_path,
                original_version=junctions.original_version,
                original_time=junctions.original_time,
            ),
            junctions,
        )
