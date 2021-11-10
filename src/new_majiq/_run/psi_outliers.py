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
from typing import List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_nonnegative_factory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand

DESCRIPTION = "Identify outliers from PSI cases (N ~ 1) vs PSI controls (N >> 1)"


def add_args(parser: argparse.ArgumentParser) -> None:
    comparison_req = parser.add_argument_group(
        "required quantification group arguments"
    )
    StorePSICovPaths = StoreRequiredUniqueActionFactory()
    comparison_req.add_argument(
        "--cases",
        type=ExistingResolvedPath,
        action=StorePSICovPaths,
        nargs="+",
        required=True,
        help="Paths to PsiCoverage files for case experiments",
    )
    comparison_req.add_argument(
        "--controls",
        type=ExistingResolvedPath,
        action=StorePSICovPaths,
        nargs="+",
        required=True,
        help="Paths to PsiCoverage files for control experiments",
    )
    # splicegraph information
    parser.add_argument(
        "splicegraph",
        metavar="SG",
        type=ExistingResolvedPath,
        help="Splicegraph on which cases, controls coverage defined",
    )
    # output information
    parser.add_argument(
        "output_zarr",
        metavar="OUT",
        type=NewResolvedPath,
        help="Path for output zarr describing results",
    )

    # compare with other splicegraphs
    sg_comparison = parser.add_argument_group("splicegraph comparison arguments")
    sg_comparison_ex = sg_comparison.add_mutually_exclusive_group(required=True)
    sg_comparison_ex.add_argument(
        "--shared-events-with",
        metavar="SG",
        type=ExistingResolvedPath,
        default=None,
        help="Only consider events which are also found in this splicegraph."
        " Annotations from this table will override those from the positional"
        " splicegraph argument. This is generally used to pass in second-pass"
        " splicegraph when analyzing quantifications from the first pass.",
    )
    sg_comparison_ex.add_argument(
        "--exclude-events-from",
        metavar="SG",
        type=ExistingResolvedPath,
        default=None,
        help="Only consider events found in the positional splicegraph argument"
        " that are not also found in this splicegraph. This is generally used"
        " to pass in the first-pass splicegraph when analyzing quantifications"
        " from the second pass.",
    )
    sg_comparison_ex.add_argument(
        "--all-events",
        default=None,
        action="store_true",
        help="Consider all events from input coverage files, matching input"
        " splicegraph. This is generally used when evaluating for outliers with"
        " only a single pass.",
    )

    # filtering criteria
    filtering = parser.add_argument_group("Outlier threshold arguments")
    filtering.add_argument(
        "--min-experiments-controls",
        metavar="X",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_OUTLIERS_MINEXPERIMENTS,
        help="min_experiments used for masking controls. If less than 1,"
        " fraction of the total number of experiments, >= 1 the absolute number"
        " that must be passed. Generally set to high percentage so that"
        " reported quantiles are reasonable (default: %(default)s)",
    )
    filtering.add_argument(
        "--alpha",
        "-a",
        dest="alpha",
        type=float,
        nargs="+",
        default=nm.constants.DEFAULT_OUTLIERS_ALPHA,
        help="Threshold for case and control quantiles in both directions"
        " (two-sided comparison)."
        " Quantiles are 0.5 * alpha and 1.0 - 0.5 * alpha. Superseded by"
        " --alpha-cases and --alpha-controls when both specified"
        " (default: %(default)s)",
    )
    filtering.add_argument(
        "--alpha-cases",
        dest="alpha_cases",
        type=float,
        nargs="+",
        default=None,
        help="Threshold for case quantiles in both directions"
        " (two-sided comparison)."
        " Quantiles are 0.5 * alpha and 1.0 - 0.5 * alpha. If specified,"
        " --alpha-controls must also be specified"
        " (default: %(default)s)",
    )
    filtering.add_argument(
        "--alpha-controls",
        dest="alpha_controls",
        type=float,
        nargs="+",
        default=None,
        help="Threshold for case quantiles in both directions"
        " (two-sided comparison)."
        " Quantiles are 0.5 * alpha and 1.0 - 0.5 * alpha. If specified,"
        " --alpha-cases must also be specified"
        " (default: %(default)s)",
    )

    # resources
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    pass


def main(sys_args: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args(sys_args)
    run(args)
    return


subcommand = GenericSubcommand(DESCRIPTION, add_args, run)


if __name__ == "__main__":
    main()
