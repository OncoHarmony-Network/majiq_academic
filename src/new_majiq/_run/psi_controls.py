"""
psi_controls.py

Compute summary of input PsiCoverage with given quantiles

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_range_factory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand

DESCRIPTION = "Summarize PSI controls for repeated comparison to cases"


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "controls",
        metavar="controls_summary",
        type=NewResolvedPath,
        help="Path for output PsiControlsSummary file",
    )
    parser.add_argument(
        "psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths to PsiCoverage files used as controls",
    )
    parser.add_argument(
        "--alpha",
        metavar="A",
        # floats on [0, 1]
        type=check_range_factory(float, 0, 1, True, True),
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=nm.constants.DEFAULT_OUTLIERS_ALPHA,
        help="Threshold for control quantiles in both directions (two-sided"
        " comparison). Quantiles are 0.5 * alpha and 1.0 - 0.5 * alpha."
        " (default: %(default)s)",
    )

    # resources
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    log = nm.logger.get_logger()
    log.info(f"Joining {len(args.psicov)} input PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    log.info(f"Summarizing {psicov}")
    nm.PsiControlsSummary.from_psicov(psicov, args.alpha).to_zarr(
        args.controls, show_progress=args.show_progress
    )
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
