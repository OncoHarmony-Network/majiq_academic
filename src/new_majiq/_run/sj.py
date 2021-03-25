"""
sj.py

Translate input BAM file to SJ file

Author: Joseph K Aicher
"""

import argparse
import new_majiq as nm
import new_majiq.constants as constants

from new_majiq.logger import get_logger
from new_majiq._run._majiq_args import check_nonnegative_factory
from pathlib import Path
from new_majiq._run._run import GenericSubcommand
from typing import (
    List,
    Optional,
)


DESCRIPTION = (
    "Translate RNA-seq alignments to raw bin reads for junctions and intronic regions"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument("bam", type=Path, help="Path to RNA-seq alignments as BAM")
    parser.add_argument(
        "splicegraph",
        type=Path,
        help="Path to splicegraph file to determine intronic regions that do"
        " not overlap exons",
    )
    parser.add_argument(
        "sj",
        type=Path,
        help="Path for SJ file with raw bin reads for junctions and introns",
    )
    parser.add_argument(
        "--strandness",
        type=str,
        default="NONE",
        choices=("NONE", "FORWARD", "REVERSE"),
        help="Strandness of input BAM (default: %(default)s)",
    )
    # TODO (enable skipping of parsing introns)
    parser.add_argument(
        "--nthreads",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_BAM_NTHREADS,
        help="Number of threads used for parsing input alignments. It is "
        " highly recommended to use multiple threads (default: %(default)s)",
    )
    parser.add_argument(
        "--update-exons",
        action="store_true",
        default=False,
        help="Experimental: Use junction coverage to definitively ignore"
        " intronic coverage in potential denovo exons (or exon extension)",
    )
    # fail if no overlapping contigs?
    disjoint_contigs_ex = parser.add_mutually_exclusive_group()
    disjoint_contigs_ex.add_argument(
        "--allow-disjoint-contigs",
        action="store_true",
        dest="allow_disjoint_contigs",
        default=constants.DEFAULT_BAM_ALLOW_DISJOINT_CONTIGS,
        help="Warn, but do not fail, when BAM has different contigs than"
        " splicegraph (default allow_disjoint_contigs = %(default)s)",
    )
    disjoint_contigs_ex.add_argument(
        "--reject-disjoint-contigs",
        action="store_false",
        dest="allow_disjoint_contigs",
        default=constants.DEFAULT_BAM_ALLOW_DISJOINT_CONTIGS,
        help="Fail when BAM has different contigs than splicegraph"
        " (default allow_disjoint_contigs = %(default)s)",
    )
    return


def run(args: argparse.Namespace) -> None:
    if not args.bam.exists():
        raise ValueError(f"Was unable to find input alignments at {args.bam}")
    if not args.splicegraph.exists():
        raise ValueError(f"Was unable to find input splicegraph at {args.splicegraph}")
    if args.sj.exists():
        raise ValueError(f"Output path {args.sj} already exists")
    log = get_logger()

    strandness = nm.ExperimentStrandness(ord(args.strandness[0]))
    log.info(f"Loading splicegraph ({args.splicegraph.resolve()})")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    log.info(f"Parsing alignments from {args.bam.resolve()} for junctions")
    sj_junctions = nm.SJJunctionsBins.from_bam(
        args.bam,
        strandness=strandness,
        nthreads=args.nthreads,
    )
    if not (set(sg.contigs.seqid) & set(sj_junctions.regions.contigs.seqid)):
        # disjoint sets of contigs from bam vs contigs
        if args.allow_disjoint_contigs:
            log.warning("Contigs from splicegraph and BAM are disjoint!")
        else:
            log.error(
                "Contigs from splicegraph and BAM are disjoint!"
                f"\n\tSplicegraph contigs = {sg.contigs.seqid}"
                f"\n\tBAM contigs = {sj_junctions.regions.contigs.seqid}"
                "\nAAdd flag `--allow-disjoint-contigs` if this is what you"
                " really want"
            )
            raise RuntimeError("Contigs from splicegraph and BAM are disjoint")
    log.info("Using gene introns/exons to define regions for intronic coverage")
    gene_introns: nm.GeneIntrons = sg.introns
    exons: nm.Exons
    if args.update_exons:
        log.info("Identifying potential denovo exons from input junctions")
        # TODO (change parameters used for reliable updated junctions?)
        updated_junctions = (
            sg.junctions.builder()
            .add_group(sg.junctions.build_group(sg.exons).add_experiment(sj_junctions))
            .get_passed()
        )
        exons = sg.exons.infer_with_junctions(updated_junctions)
    else:
        exons = sg.exons
    log.info(f"Parsing alignments from {args.bam.resolve()} for introns")
    sj_introns = nm.SJIntronsBins.from_bam(
        args.bam,
        total_bins=sj_junctions.total_bins,
        exons=exons,
        gene_introns=gene_introns,
        strandness=strandness,
        nthreads=args.nthreads,
    )
    log.info(f"Saving junction and intron coverage to {args.sj.resolve()}")
    sj_junctions.to_zarr(args.sj)
    sj_introns.to_zarr(args.sj)
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
