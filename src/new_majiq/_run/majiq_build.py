#!/usr/bin/env python
"""
majiq_build.py

Provides entry-point to different majiq-build scripts

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict

from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build import subcommand as build
from new_majiq._run.combine import subcommand as combine
from new_majiq._run.gff3 import subcommand as gff3
from new_majiq._run.psi_coverage import subcommand as psi_coverage
from new_majiq._run.sg_coverage import subcommand as sg_coverage
from new_majiq._run.sj import subcommand as sj

SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "gff3": gff3,
    "sj": sj,
    "update": build,
    "combine": combine,
    "psi-coverage": psi_coverage,
    "sg-coverage": sg_coverage,
}


def main() -> None:
    """Entry-point into multiple tools using subcommands"""
    # build parser
    parser = argparse.ArgumentParser(
        description="Tools to analyze RNA splicing by building splicegraphs"
        " over genes and assessing coverage over junctions and retained introns"
    )
    # add subparsers
    subparsers = parser.add_subparsers(required=True, help="")
    for src_name, src_module in SUBPARSER_SOURCES.items():
        src_parser = subparsers.add_parser(
            src_name,
            help=src_module.DESCRIPTION,
            description=src_module.DESCRIPTION,
        )
        src_parser.set_defaults(func=src_module.run)
        src_module.add_args(src_parser)

    # check length of input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # parse arguments now
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()