#!/usr/bin/env python
"""
majiq_moccasin.py

Provides entry-point to different majiq scripts for MOCCASIN-powered batch
correction

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict

from new_majiq._run._run import GenericSubcommand
from new_majiq._run.moccasin import subcommand_coverage_infer as moccasin_coverage_infer
from new_majiq._run.moccasin import subcommand_coverage_model as moccasin_coverage_model
from new_majiq._run.moccasin import subcommand_factors_infer as moccasin_factors_infer
from new_majiq._run.moccasin import subcommand_factors_model as moccasin_factors_model

SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "factors-model": moccasin_factors_model,
    "factors-infer": moccasin_factors_infer,
    "coverage-model": moccasin_coverage_model,
    "coverage-infer": moccasin_coverage_infer,
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
