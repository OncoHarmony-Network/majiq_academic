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
from new_majiq._run.build_args import sj_strandness_args
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
    update_exons = parser.add_argument_group("update exons arguments")
    update_exons.add_argument(
        "--no-update-exons",
        action="store_true",
        default=False,
        help="Do not use junction coverage to update splicegraph/regions"
        " used for intronic coverage."
        " Use splicegraph as passed in in order to determine regions for"
        " intronic coverage"
        " (default: update splicegraph using --update-minreads/--update-minpos)",
    )
    update_exons.add_argument(
        "--update-minreads",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_BUILD_MINDENOVO,
        metavar="N",
        help="Set the minimum number of reads to pass a novel junction."
        " Novel splicesites used to split regions for intronic coverage."
        " (default: %(default)s)",
    )
    update_exons.add_argument(
        "--update-minpos",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_BUILD_MINPOS,
        metavar="N",
        help="Set the minimum number of nonzero positions to pass a novel"
        " junction."
        " Novel splicesites used to split regions for intronic coverage."
        " (default: %(default)s)",
    )
    sj_strandness_args(parser)
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
    update_exons_thresholds: Optional[nm.ExperimentThresholds] = None
    if args.no_update_exons:
        log.debug("Not updating splicegraph with junction coverage")
    else:
        log.debug(
            "Update splicegraph with junction coverage; add novel junctions if"
            " number of reads >= %d, number of nonzero positions >= %d",
            args.update_minreads,
            args.update_minpos,
        )
        update_exons_thresholds = nm.ExperimentThresholds(
            minreads=args.update_minreads,
            mindenovo=args.update_minreads,
            minpos=args.update_minpos,
        )
    # load junctions, introns
    sj = nm.SJExperiment.from_bam(
        args.bam,
        sg,
        args.strandness,
        nthreads=args.nthreads,
        allow_disjoint_contigs=args.allow_disjoint_contigs,
        auto_minreads=args.auto_minreads,
        auto_minjunctions=args.auto_minjunctions,
        auto_mediantolerance=args.auto_mediantolerance,
        update_exons_thresholds=update_exons_thresholds,
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
