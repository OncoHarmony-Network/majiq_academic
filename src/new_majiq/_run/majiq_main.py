#!/usr/bin/env python
"""
majiq_main.py

Provides entry-point to primary majiq scripts

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict

from new_majiq._run._run import GenericSubcommand
from new_majiq._run.build_pipeline import subcommand as build_pipeline
from new_majiq._run.cite import subcommand as cite
from new_majiq._run.deltapsi import subcommand as deltapsi
from new_majiq._run.heterogen import subcommand as heterogen
from new_majiq._run.moccasin import subcommand_pipeline as moccasin_pipeline
from new_majiq._run.quantify import subcommand as psi

SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "build": build_pipeline,
    "moccasin": moccasin_pipeline,
    "psi": psi,
    "deltapsi": deltapsi,
    "heterogen": heterogen,
    "cite": cite,
}


def main() -> None:
    """Entry-point into multiple tools using subcommands"""
    # build parser
    parser = argparse.ArgumentParser(
        description="Tools to detect, quantify, and analyze RNA splicing",
        epilog="More with majiq-{build,moccasin,quantify,mendelian} or new-majiq",
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
