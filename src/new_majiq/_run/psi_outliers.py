"""
psi_outliers.py

Identify aberrant splicing events of interests in cases vs controls

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional, cast

from dask.delayed import Delayed
from dask.distributed import progress

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_range_factory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.core._workarounds import _load_zerodim_variables

DESCRIPTION = "Identify outliers from PSI cases (N ~ 1) vs PSI controls (N >> 1)"


def add_args(parser: argparse.ArgumentParser) -> None:
    # input information
    parser.add_argument(
        "cases",
        metavar="cases_psicov",
        type=ExistingResolvedPath,
        help="Path to PsiCoverage file with case experiments",
    )
    parser.add_argument(
        "controls",
        metavar="controls_summary",
        type=ExistingResolvedPath,
        help="Path to PsiControlsSummary file for comparison to controls",
    )
    # output information
    parser.add_argument(
        "outliers",
        metavar="out_outliers",
        type=NewResolvedPath,
        help="Path for output zarr with PsiOutliers computations for specified"
        " parameters",
    )
    # thresholds for case (posterior) quantiles
    parser.add_argument(
        "--alpha",
        metavar="A",
        # floats on [0, 1]
        type=check_range_factory(float, 0, 1, True, True),
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="Threshold for case quantiles in both directions (two-sided"
        " comparison). Quantiles are 0.5 * alpha and 1.0 - 0.5 * alpha."
        " If not specified, use same values as from input controls"
        " (default: %(default)s)",
    )

    # resources
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    log = nm.logger.get_logger()
    log.info(f"Loading cases from {args.cases}")
    cases = nm.PsiCoverage.from_zarr(args.cases)
    log.info(f"Loading controls summary from {args.controls}")
    controls = nm.PsiControlsSummary.from_zarr(args.controls)
    log.info(f"Comparing cases {cases} vs {controls}")
    outliers = nm.PsiOutliers(cases, controls)
    ds = outliers.dataset(cases_alpha=args.alpha)
    log.info(f"Saving comparisons to {args.outliers}")
    save_ds_future = cast(
        Delayed,
        ds.pipe(_load_zerodim_variables).to_zarr(
            args.outliers,
            mode="w",
            group=nm.constants.NC_PSIOUTLIERS,
            consolidated=False,
            compute=False,
        ),
    )
    if args.show_progress:
        save_ds_future = save_ds_future.persist()
        progress(save_ds_future)
    else:
        save_ds_future.compute()
    outliers.events.chunk(outliers.events.sizes).to_zarr(
        args.outliers, mode="a", group=nm.constants.NC_EVENTS, consolidated=True
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
