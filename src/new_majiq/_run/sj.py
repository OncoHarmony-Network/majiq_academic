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
    sj_strandness_args(parser)
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
