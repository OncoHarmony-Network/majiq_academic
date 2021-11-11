"""
psi_outliers.py

Identify aberrant splicing events of interests in cases vs controls

Author: Joseph K Aicher

Parameters
----------
splicegraph
-cases PSICOV [PSICOV ...]
-controls PSICOV [PSICOV ...]
(if shared-events-with, use that splicegraph's annotations for output table)
(--shared-events-with SG | --exclude-events-from SG | --all-events)
--control-min-experiments (0.9 default)
--alpha Q [Q ...] ([0.05] default). Used for quantiles, similar to two-sided...
--controls-alpha Q [Q ...] (None -> use cases-alpha, then alpha)
--cases-alpha Q [Q ...] (None -> use controls-alpha, then alpha)
(dask resources)

What do I want to get back?
(cases_prefix) psi_mean, psi_std
(cases_prefix, cases_quantiles) case_psi_quantile
controls_psi_median
(controls_quantiles) controls_psi_quantile
(cases_prefix, cases_alpha, controls_alpha) psi_gap

Note that {cases,controls}_quantiles are actually just for computation, and
really a stacked/unstacked version of ({cases,controls}_alpha, lower[bool]).
psi_gap is computed by taking the maximum of the lower - higher in either
direction, and of 0 (if they overlap, both differences will be negative.

Easy enough to get all of these values.
xr.Dataset(dict(), dict(alpha=[0.05, 0.1, 0.3], is_lb=[False, True])).pipe(lambda ds: (q := 0.5 * ds.alpha).where(ds.is_lb, 1 - q)).rename("q").to_dataset().set_coords("q")
Should check that input values for alpha are between 0 and 1

We can make psi_gap unitless by scaling by the within-group difference in psi:
psi_gap / ((case_quantiles[lb=False] - case_quantiles[lb=True]) + (control_quantiles[lb=False] - case_quantiles[lb=True])).
This is a derived quantity that should help account/filter out junctions with
high variability.

Something which we would have to do more implementation for would be filtering
out low coverage LSVs in the context of their splicing modules.
We would have to set up approach for identifying the modules in new-majiq.
We would have to define coverage in splicing module relative to LSV.
Easiest thing to do would probably be to summarize coverage at LSVs at start/end
of each module as reference (as the source/sink of flow of coverage through
splicing decisions).
We would have to be careful about how we defined things, etc.
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
from new_majiq._workarounds import _load_zerodim_variables

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
