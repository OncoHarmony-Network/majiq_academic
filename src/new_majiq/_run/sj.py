"""
sj.py

Translate input BAM file to SJ file

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    check_nonnegative_factory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = (
    "Translate RNA-seq alignments to raw bin reads for junctions and intronic regions"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "bam", type=ExistingResolvedPath, help="Path to RNA-seq alignments as BAM"
    )
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        help="Path to splicegraph file to determine intronic regions that do"
        " not overlap exons",
    )
    parser.add_argument(
        "sj",
        type=NewResolvedPath,
        help="Path for SJ file with raw bin reads for junctions and introns",
    )
    strandness = parser.add_argument_group("strandness arguments")
    strandness.add_argument(
        "--strandness",
        type=str,
        default="AUTO",
        choices=(
            "AUTO",
            nm.ExperimentStrandness.NONE.name,
            nm.ExperimentStrandness.FORWARD.name,
            nm.ExperimentStrandness.REVERSE.name,
        ),
        help="Strandness of input BAM."
        " AUTO = automatically detect strand (use median ratio of forward vs"
        " reverse stranded reads at annotated junctions). (default: %(default)s)",
    )
    strandness.add_argument(
        "--auto-minreads",
        metavar="N",
        type=check_nonnegative_factory(int, reject_zero=True),
        default=nm.constants.DEFAULT_BAM_STRAND_MINREADS,
        help="For automatic detection of strand. Only consider evidence from"
        " splicegraph junctions with at least this many total (unstranded) reads"
        " (default: %(default)s)",
    )
    strandness.add_argument(
        "--auto-minjunctions",
        metavar="N",
        type=check_nonnegative_factory(int, reject_zero=True),
        default=nm.constants.DEFAULT_BAM_STRAND_MINJUNCTIONS,
        help="For automatic detection of strand. Infer unstranded if the number"
        " of splicegraph junctions with sufficient reads is less than this argument"
        " (default: %(default)s)",
    )
    strandness.add_argument(
        "--auto-mediantolerance",
        metavar="X",
        type=check_nonnegative_factory(float, reject_zero=True),
        default=nm.constants.DEFAULT_BAM_STRAND_MINDEVIATION,
        help="For automatic detection of strand. Infer unstranded if the median"
        " proportion over junctions of forward strand reads vs all reads"
        " deviates from 0.5 by only this amount (default: %(default)s)",
    )
    parser.add_argument(
        "--update-exons",
        action="store_true",
        default=False,
        help="Experimental: Use junction coverage to definitively ignore"
        " intronic coverage in potential denovo exons (or exon extension)",
    )
    # fail if no overlapping contigs?
    disjoint_contigs = parser.add_argument_group("disjoint contigs arguments")
    disjoint_contigs_ex = disjoint_contigs.add_mutually_exclusive_group()
    disjoint_contigs_ex.add_argument(
        "--allow-disjoint-contigs",
        action="store_true",
        dest="allow_disjoint_contigs",
        default=nm.constants.DEFAULT_BAM_ALLOW_DISJOINT_CONTIGS,
        help="Warn, but do not fail, when BAM has different contigs than"
        " splicegraph (default allow_disjoint_contigs = %(default)s)",
    )
    disjoint_contigs_ex.add_argument(
        "--reject-disjoint-contigs",
        action="store_false",
        dest="allow_disjoint_contigs",
        default=nm.constants.DEFAULT_BAM_ALLOW_DISJOINT_CONTIGS,
        help="Fail when BAM has different contigs than splicegraph"
        " (default allow_disjoint_contigs = %(default)s)",
    )
    resources_args(parser, use_dask=False)
    return


def run(args: argparse.Namespace) -> None:
    log = get_logger()

    log.info(f"Loading splicegraph ({args.splicegraph})")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    # load junctions, introns
    sj = nm.SJExperiment.from_bam(
        args.bam,
        sg,
        args.strandness,
        update_exons=args.update_exons,
        nthreads=args.nthreads,
        allow_disjoint_contigs=args.allow_disjoint_contigs,
        auto_minreads=args.auto_minreads,
        auto_minjunctions=args.auto_minjunctions,
        auto_mediantolerance=args.auto_mediantolerance,
    )
    log.info(f"Saving junction and intron coverage to {args.sj}")
    sj.to_zarr(args.sj)
    return


def main(sys_args: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args(sys_args)
    run(args)
    return


subcommand = GenericSubcommand(DESCRIPTION, add_args, run)


if __name__ == "__main__":
    main()
