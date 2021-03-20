#!/usr/bin/env python
"""
tools.py

Provides entry-point to different majiq-tools scripts

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict
from new_majiq._run._run import GenericSubcommand
from new_majiq._run.cite import subcommand as cite
from new_majiq._run.gff3 import subcommand as gff3
from new_majiq._run.sj import subcommand as sj
from new_majiq._run.build import subcommand as build
from new_majiq._run.build_group import subcommand as build_group
from new_majiq._run.simplify import subcommand as simplify
from new_majiq._run.lsv_coverage import subcommand as lsv_coverage
from new_majiq._run.sg_coverage import subcommand as sg_coverage
from new_majiq._run.sg_coverage_summarize import subcommand as sg_coverage_summarize
from new_majiq._run.quantify import subcommand as quantify
from new_majiq._run.legacy_psi import subcommand as legacy_psi


SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "gff3": gff3,
    "sj": sj,
    "build": build,
    "build-group": build_group,
    "simplify": simplify,
    "lsv-coverage": lsv_coverage,
    "sg-coverage": sg_coverage,
    "sg-coverage-summary": sg_coverage_summarize,
    "quantify": quantify,
    "legacy-psi": legacy_psi,
    "cite": cite,
}


def main() -> None:
    """Entry-point into multiple tools using subcommands"""
    # build parser
    parser = argparse.ArgumentParser(
        description="Tools to detect, quantify, and analyze RNA splicing"
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
