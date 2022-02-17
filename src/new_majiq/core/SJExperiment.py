"""
SJExperiment.py

Convenience class for representing SJIntronsBins and SJJunctionsBins for the
same experiment

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Final, Optional, Union

import numpy as np

import new_majiq.constants as constants
from new_majiq.internals import ExperimentStrandness, ExperimentThresholds
from new_majiq.logger import get_logger

from .Exons import Exons
from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions
from .SJIntrons import SJIntrons
from .SJIntronsBins import SJIntronsBins
from .SJJunctionsBins import SJJunctionsBins
from .SpliceGraph import SpliceGraph
from .SpliceGraphReads import SpliceGraphReads


class SJExperiment(object):
    """Spliced junction and retained intron coverage for the same experiment.

    Junctions and retained introns are defined over contigs, not genes
    (i.e. :class:`ContigIntrons`, :class:`ContigJunctions`, not
    :class:`GeneIntrons`, :class:`GeneJunctions`).
    Coverage assigned to these regions is later mapped to their gene
    equivalents (matching contig and strand; matching coordinates for
    junctions, matching annotated status and overlapping coordinates for
    introns).

    Coverage is obtained over unique bins.
    Bins represent positions at which the start of the feature would intersect
    the aligned read.
    In the case of junctions, there is a single bin for each possible position,
    so the total number of bins/positions is determined by the maximum read
    length.
    For introns, nonzero length increases the total number of possible
    positions.
    In order to treat introns and junctions similarly, MAJIQ groups together
    contiguous positions into the same number of bins as junctions.

    This is represented by:

    - :class:`SJIntronsBins` with intron coverage and assigned regions
      (:attr:`SJExperiment.introns`),
    - :class:`SJJunctionsBins` with junction coverage and assigned regions
      (:attr:`SJExperiment.junctions`).

    Parameters
    ----------
    introns: SJIntronsBins
    junctions: SJJunctionsBins

    See Also
    --------
    SJExperiment.from_bam : Load :class:`SJExperiment` from BAM file
    SJExperiment.from_zarr : Load :class:`SJExperiment` saved to Zarr file
    """

    def __init__(self, introns: SJIntronsBins, junctions: SJJunctionsBins):
        """Initialize :class:`SJExperiment` with intron and junction coverage

        Parameters
        ----------
        introns: SJIntronsBins
        junctions: SJJunctionsBins

        Notes
        -----
        Enforces requirement that introns and junctions share the same original
        (bam) path and version of majiq used to do parsing
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

    def __repr__(self) -> str:
        """String representation of SJExperiment"""
        return f"SJExperiment(original_path=[...]/{Path(self.original_path).name})"

    @property
    def introns(self) -> SJIntronsBins:
        """:py:class:`SJIntronsBins` with intron coverage for this experiment"""
        return self._introns

    @property
    def junctions(self) -> SJJunctionsBins:
        """:py:class:`SJJunctionsBins` with junctioncoverage for this experiment"""
        return self._junctions

    @property
    def original_path(self) -> str:
        """Path to BAM file that was processed for coverage"""
        return self.junctions.original_path

    @property
    def original_version(self) -> str:
        """Version of MAJIQ that was used to process the coverage"""
        return self.junctions.original_version

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "SJExperiment":
        """Load :py:class:`SJExperiment` from specified path

        Parameters
        ----------
        path: Union[str, Path]
            Path where coverage is stored in zarr format

        Returns
        -------
        SJExperiment
        """
        return SJExperiment(
            SJIntronsBins.from_zarr(path), SJJunctionsBins.from_zarr(path)
        )

    @classmethod
    def original_path_from_zarr(cls, path: Union[str, Path]) -> str:
        """Extract bam path from Zarr without loading entire file

        Parameters
        ----------
        path: Union[str, Path]
            Path where coverage is stored in zarr format

        Returns
        -------
        str
            Original BAM path used to load the :class:`SJExperiment` saved at
            path
        """
        return SJJunctionsBins.original_path_from_zarr(path)

    def to_zarr(self, path: Union[str, Path], consolidated: bool = True) -> None:
        """Save :py:class:`SJExperiment` to specified path

        Parameters
        ----------
        path: Union[str, Path]
            Path for zarr store with intron and junction coverage
        consolidated: bool
            Finalize metadata for zarr file
        """
        self.junctions.to_zarr(path, consolidated=False)
        self.introns.to_zarr(path, consolidated=True, check_experiment_if_exists=False)
        return

    @classmethod
    def from_bam(
        cls,
        path: Union[str, Path],
        sg: SpliceGraph,
        strandness: Union[str, ExperimentStrandness] = "auto",
        nthreads: int = constants.DEFAULT_BAM_NTHREADS,
        allow_disjoint_contigs: bool = constants.DEFAULT_BAM_ALLOW_DISJOINT_CONTIGS,
        auto_minreads: int = constants.DEFAULT_BAM_STRAND_MINREADS,
        auto_minjunctions: int = constants.DEFAULT_BAM_STRAND_MINJUNCTIONS,
        auto_mediantolerance: float = constants.DEFAULT_BAM_STRAND_MINDEVIATION,
        update_exons_thresholds: Optional[
            ExperimentThresholds
        ] = constants.DEFAULT_BUILD_EXP_THRESHOLDS,
    ) -> "SJExperiment":
        """Load :class:`SJExperiment` from BAM file

        Load intron and junction coverage from BAM file using splicegraph to
        define intron regions for coverage. If strandness is "auto", use
        coverage assigned to each strand on known junctions to infer
        experimental strandness.

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
        update_exons_thresholds: Optional[ExperimentThresholds]
            If specified, use denovo junction coverage passing thresholds to
            split intronic coverage by novel splice sites. If None, do not
            update junctions from splicegraph using input junction coverage.

        Returns
        -------
        SJExperiment
            intron and junction coverage from specified BAM file
        """
        log = get_logger()
        log.info("Loading SJExperiment from aligned reads in %s", path)
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
        log.debug("Parsing SJJunctions from %s", path)
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
            log.debug("Inferring strandness comparing counts on splicegraph junctions")
            junctions = cls.detect_strand(
                junctions,
                sg,
                minreads=auto_minreads,
                minjunctions=auto_minjunctions,
                mindeviation=auto_mediantolerance,
            )
            log.info(f"Inferred strandness: {junctions.strandness.name}")
        # determine introns
        log.debug("Using gene introns/exons to define regions for intronic coverage")
        gene_introns: GeneIntrons = sg.introns
        gene_junctions: GeneJunctions = sg.junctions
        exons: Exons = sg.exons
        if update_exons_thresholds:
            log.debug("Updating splicegraph junctions using SJ coverage")
            gene_junctions = (
                sg.junctions.builder()
                .add_group(
                    sg.junctions.build_group(exons).add_experiment(
                        junctions, thresholds=update_exons_thresholds, add_denovo=True
                    )
                )
                .get_passed()
            )
        # update exons using annotated boundaries, splicesites from gene_junctions
        exons = exons.infer_with_junctions(gene_junctions, full_inference=False)
        log.debug("Parsing SJIntrons from %s", path)
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
        """Return :class:`SJJunctionsBins` with appropriate strandness

        Given `sj_junctions`, flip or remove strandness as appropriate to match
        junctions in `sg`.
        This is done by comparing the number of reads assigned to each junction
        in `sg` when junctions are treated as forward vs reverse stranded.
        MAJIQ restricts analysis to junctions with more than `minreads` reads.
        Then, it calls an experiment as stranded if the median ratio of stranded
        vs total reads is within `mindeviation` of 0.5.
        If stranded, it chooses the strand direction consistent with the median
        ratio.

        Parameters
        ----------
        sj_junctions: SJJunctionsBins
            junction coverage to align strandness for relative to splicegraph
        sg: SpliceGraph
            splicegraph with junctions which are checked against junction
            coverage
        """
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
        ).junctions_values
        sj_flipped = SJExperiment(sj.introns, sj_junctions.flip_strand())
        reads_flipped = SpliceGraphReads._internals_from_connections_and_sj(
            empty_sg_introns, sg.junctions, sj_flipped
        ).junctions_values

        # compute ratios for junctions with at least minreads
        reads_total = reads + reads_flipped
        minreads_mask = reads_total >= minreads
        ratios = reads[minreads_mask] / reads_total[minreads_mask]

        # check that we had enough junctions
        if len(ratios) < minjunctions:
            log.info(
                "Too few junctions with enough reads to infer strandness"
                " (%d < minjunctions=%d; with minreads=%d)",
                len(ratios),
                minjunctions,
                minreads,
            )
            return sj_junctions.to_unstranded()

        # compute ratio
        median_ratio = np.median(ratios)
        deviation = np.abs(median_ratio - 0.5)
        log.debug(
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
        """Return :class:`SJExperiment` with junction coverage, no introns"""
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
