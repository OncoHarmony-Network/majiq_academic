"""
gff3.py

Translate input GFF3 file to base annotated splicegraph file

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = "Translate input GFF3 file to base splicegraph file"


def add_args(parser: argparse.ArgumentParser) -> None:
    """add arguments to parser"""
    parser.add_argument(
        "gff3",
        type=ExistingResolvedPath,
        help="Path to GFF3 file (uncompressed or gzipped) to convert",
    )
    parser.add_argument(
        "splicegraph",
        type=NewResolvedPath,
        help="Path to save resulting splicegraph. Fails if path already exists.",
    )
    # annotate ir options
    introns_ex = parser.add_mutually_exclusive_group()
    introns_ex.add_argument(
        "--ignore_ir",
        action="store_false",
        dest="process_ir",
        default=nm.constants.DEFAULT_BUILD_PROCESS_IR,
        help="Ignore annotated introns when processing annotations"
        " (default process_ir=%(default)s)",
    )
    introns_ex.add_argument(
        "--process_ir",
        action="store_true",
        dest="process_ir",
        default=nm.constants.DEFAULT_BUILD_PROCESS_IR,
        help="Ensure annotated introns are processed (default process_ir=%(default)s)",
    )
    # enable configuration of GFF3 parsing
    gff3_options = parser.add_argument_group("GFF3 parsing options arguments")
    gff3_base_types = gff3_options.add_mutually_exclusive_group()
    gff3_base_types.add_argument(
        "--features-default",
        action="store_true",
        default=True,
        help="Initialize parsing of GFF3 transcripts using MAJIQ defaults"
        " for select GFF3 types."
        " Use --types-* to map additional GFF3 types. (default)",
    )
    gff3_base_types.add_argument(
        "--features-none",
        dest="features_default",
        action="store_false",
        default=True,
        help="Initialize parsing of GFF3 transcripts with no GFF3 types."
        " Use --types-* to map additional GFF3 types"
        " (default: features-default)",
    )
    gff3_base_types.add_argument(
        "--features-tsv",
        metavar="tsv",
        default=None,
        type=ExistingResolvedPath,
        help="Initialize parsing of GFF3 transcripts with first two columns of"
        " tab delimited file (first column: GFF3 type to map, second column:"
        f" action for that GFF3 type, from {nm.GFF3TypesMap.VALID_ACTIONS()})"
        " (default: features-default)",
    )
    StoreGFF3Types = StoreRequiredUniqueActionFactory()
    gff3_options.add_argument(
        "--types-genes",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Map specified types from GFF3 as genes if found as ancestor of exon",
    )
    gff3_options.add_argument(
        "--types-transcripts",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Map specified types from GFF3 as transcripts if found as parent of exon",
    )
    gff3_options.add_argument(
        "--types-exons",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Map specified types from GFF3 as exons",
    )
    gff3_options.add_argument(
        "--types-silent",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Specified types from GFF3 will be ignored silently if found as"
        " parent (ignored transcript) or top-level ancestor (ignored gene)",
    )
    gff3_options.add_argument(
        "--types-hard-skip",
        metavar="T",
        nargs="*",
        action=StoreGFF3Types,
        type=str,
        default=list(),
        help="Specified types from GFF3 should never be an ancestor of exons"
        " and can be completely skipped for potential performance improvements",
    )
    return


def run(args: argparse.Namespace) -> None:
    log = get_logger()
    log.info(f"Processing GFF3 {args.gff3} to create annotated splicegraph")

    def log_function(gff3_type: str, missing_reason: str, count: int) -> None:
        """log function for :meth:`SpliceGraph.from_gff3`

        log function for :meth:`SpliceGraph.from_gff3` as warnings for
        non-silent skipped parents/top-level ancestores
        """
        log.warning(
            "GFF3 type '%s' skipped %d times (%s)",
            gff3_type,
            count,
            missing_reason,
        )
        return

    # get GFF3 types for parsing
    gff3_types: nm.GFF3TypesMap
    if args.features_default and not args.features_tsv:
        gff3_types = nm.GFF3TypesMap()  # use default values
    else:
        gff3_types = nm.GFF3TypesMap({})  # start with empty values
        if args.features_tsv:
            # process each line in args.features_tsv
            for x in open(args.features_tsv, "r"):
                gff3_key, gff3_action, *_ = x.strip().split("\t")
                if gff3_key in gff3_types:
                    raise ValueError(
                        f"{args.features_tsv} specifies {gff3_key} more than once"
                    )
                gff3_types[gff3_key] = gff3_action
    gff3_types = nm.GFF3TypesMap.from_types_sets(
        gene_types=gff3_types.gene_types() | set(args.types_genes),
        transcript_types=gff3_types.transcript_types() | set(args.types_transcripts),
        exon_types=gff3_types.exon_types() | set(args.types_exons),
        silent_types=gff3_types.silent_types() | set(args.types_silent),
        hard_skip_types=gff3_types.hard_skip_types() | set(args.types_hard_skip),
    )
    log.info("Parsing GFF3 with %s", gff3_types)

    sg = nm.SpliceGraph.from_gff3(
        args.gff3,
        process_ir=args.process_ir,
        gff3_types=gff3_types,
        log_function=log_function,
    )
    log.info(f"Saving annotated splicegraph to {args.splicegraph}")
    sg.to_zarr(args.splicegraph)
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
