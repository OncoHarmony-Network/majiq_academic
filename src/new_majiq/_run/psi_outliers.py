"""
psi_outliers.py

Identify aberrant splicing events of interests in cases vs controls

Author: Joseph K Aicher
"""

import argparse
import json
import sys
from typing import Any, Dict, List, Optional

import pandas as pd

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_range_factory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand

DESCRIPTION = "Identify outliers from PSI cases (N ~ 1) vs PSI controls (N >> 1)"


def add_args(parser: argparse.ArgumentParser) -> None:
    StoreUniqueInputPaths = StoreRequiredUniqueActionFactory()
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        action=StoreUniqueInputPaths,
        help="Path for input SpliceGraph on which comparisons will be made",
    )
    parser.add_argument(
        "controls",
        type=ExistingResolvedPath,
        action=StoreUniqueInputPaths,
        help="Path for PsiControlsSummary defined over splicegraph",
    )
    parser.add_argument(
        "cases",
        type=ExistingResolvedPath,
        action=StoreUniqueInputPaths,
        help="Path for PsiCoverage for cases defined over splicegraph",
    )
    parser.add_argument(
        "--output-tsv",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path for output TSV file (default: stdout)",
    )
    parser.add_argument(
        "--events-summary",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=None,
        help="Path for summary per splicing event (default: None)",
    )
    parser.add_argument(
        "--genes-summary",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=None,
        help="Path for summary per splicing gene (default: None)",
    )
    parser.add_argument(
        "--pass-previous",
        metavar=("base-splicegraph", "base-controls", "base-cases"),
        nargs=3,
        type=ExistingResolvedPath,
        action=StoreUniqueInputPaths,
        default=None,
        help="Paths for input SpliceGraph, PsiControlsSummary, and PsiCoverage"
        " for previous pass. If specified, comparisons from both passes"
        " will be merged on events defined `splicegraph`."
        " Splicegraphs are required to have the same genes."
        " Controls are required to be defined relative to the same"
        " experiments."
        " For cases, experiments in `cases` must also be in `base-cases`."
        " This helps combine results from two passes of majiq build, one with"
        " just controls that can be shared between analyses, and another with"
        " the cases added."
        " (default: only process one set of quantifications)",
    )
    parser.add_argument(
        "--min-experiments",
        metavar="X",
        type=check_range_factory(float, minval=0, mininclude=False),
        default=nm.constants.DEFAULT_OUTLIERS_MINEXPERIMENTS,
        help="The fraction (value < 1) or absolute number (value >= 1) of"
        " experiments that must pass individually for an event in the controls"
        " for a comparison to be made (default: %(default)s)",
    )
    parser.add_argument(
        "--select-alpha",
        metavar="A",
        type=check_range_factory(float, 0, 1, True, True),
        default=None,
        help="Select single threshold on controls/case quantils in both"
        " directions (two-sided comparison)."
        " Quantiles are 0.5*alpha and 1.0-0.5*alpha."
        " Selected threshold must have been used when creating controls"
        " (default: use all values of alpha available in controls)",
    )
    # resources
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    log = nm.logger.get_logger()
    metadata: Dict[str, Any] = dict()
    metadata["command"] = " ".join(sys.argv)
    metadata["version"] = nm.__version__

    log.debug("Loading input splicegraphs")
    sg_list: List[nm.SpliceGraph] = []
    use_genes: Optional[nm.Genes] = None
    if args.pass_previous:
        log.info("Loading previous pass splicegraph from %s", args.pass_previous[0])
        sg_list.append(nm.SpliceGraph.from_zarr(args.pass_previous[0]))
        use_genes = sg_list[0].genes
    log.info("Loading input splicegraph from %s", args.splicegraph)
    sg_list.append(nm.SpliceGraph.from_zarr(args.splicegraph, genes=use_genes))

    log.debug("Loading input controls")
    controls_list: List[nm.PsiControlsSummary] = []
    if args.pass_previous:
        log.info("Loading previous pass controls from %s", args.pass_previous[1])
        controls_list.append(nm.PsiControlsSummary.from_zarr(args.pass_previous[1]))
    log.info("Loading input controls from %s", args.controls)
    controls_list.append(nm.PsiControlsSummary.from_zarr(args.controls))
    if args.pass_previous:
        log.debug("Verifying that controls were defined with the same experiments")
        if missing := set(controls_list[0].prefixes).symmetric_difference(
            set(controls_list[1].prefixes)
        ):
            raise ValueError(
                "Input controls are defined with different experiments"
                f" between passes (not shared: {missing})"
            )
        log.debug("Verifying that controls were defined with the same values of alpha")
        if set(controls_list[0].alpha.values) != set(controls_list[1].alpha.values):
            raise ValueError(
                "Input controls are defined with different values of alpha"
                f" (pass1: {set(controls_list[0].alpha.values)},"
                f" pass2: {set(controls_list[1].alpha.values)})"
            )
    if args.select_alpha and args.select_alpha not in controls_list[0].alpha.values:
        raise ValueError(
            f"Selected value of alpha ({args.select_alpha}) is not defined for"
            f" controls (defined: {controls_list[0].alpha.values.tolist()})"
        )
    metadata["n_controls"] = controls_list[0].num_prefixes
    log.info("Quantifying potential outliers relative to controls %s", controls_list)

    log.debug("Loading input cases")
    cases_list: List[nm.PsiCoverage] = []
    if args.pass_previous:
        log.info("Loading previous pass cases from %s", args.pass_previous[2])
        cases_list.append(nm.PsiCoverage.from_zarr(args.pass_previous[2]))
    log.info("Loading input cases from %s", args.cases)
    cases_list.append(nm.PsiCoverage.from_zarr(args.cases))
    if args.pass_previous:
        log.debug("Verifying that previous pass cases has experiments from input cases")
        if missing := set(cases_list[1].prefixes).difference(cases_list[0].prefixes):
            raise ValueError(
                "Previous pass cases does not have coverage for same"
                f" experiments as input cases ({missing = })"
            )
        log.debug("Taking subset of experiments in base cases")
        cases_list[0] = cases_list[0][cases_list[1].prefixes]
    metadata["n_cases"] = cases_list[0].num_prefixes
    log.info("Quantifying potential outliers from cases %s", cases_list)

    log.info("Quantifying differences between cases vs controls")
    df_list: List[pd.DataFrame] = [
        nm.PsiOutliers(cases, controls).to_dataframe(
            controls_min_experiments=args.min_experiments,
            alpha=args.select_alpha,
            show_progress=args.show_progress,
        )
        for cases, controls in zip(cases_list, controls_list)
    ]
    log.info("Annotating quantifications with splicegraph information")
    df: pd.DataFrame
    if args.pass_previous:
        # keep any information compatible with current splicegraph, even if it
        # isn't a strict LSV in current splicegraph
        events = sg_list[-1].exon_connections.lsvs(
            nm.constants.SelectLSVs.PERMISSIVE_LSVS
        )
        # get events found in controls
        events_list = [
            controls.get_events(sg.introns, sg.junctions)
            for controls, sg in zip(controls_list, sg_list)
        ]
        # merge df_list onto events
        df = events.merge_dataframes(df_list, events_list)
    else:
        df = (
            controls_list[0]
            .get_events(sg_list[0].introns, sg_list[0].junctions)
            .ec_dataframe.join(df_list[0], how="inner", on="ec_idx")
        )

    try:
        output_name = args.output_tsv.name
    except AttributeError:
        output_name = args.output_tsv
    log.info("Saving metadata and table of comparisons to %s", output_name)
    metadata_json = json.dumps(metadata, sort_keys=True, indent=4)
    args.output_tsv.write("# {}\n".format(metadata_json.replace("\n", "\n# ")))
    (
        df
        # manually format probability columns
        .pipe(
            lambda df: df.assign(
                **{
                    col: df[col].apply(lambda x: f"{x:.3e}")
                    for col in df.columns
                    if "probability" in col
                }
            )
        )
        # numeric columns need at most 4 digits precision
        .round(4)
        # save result in TSV format
        .to_csv(args.output_tsv, sep="\t", index=not args.splicegraph)
    )

    if args.events_summary or args.genes_summary:
        df_events = nm.PsiOutliers.summarize_df_events(df)
        if args.events_summary:
            try:
                output_name = args.events_summary.name
            except AttributeError:
                output_name = args.events_summary
            log.info(
                "Saving metadata and summary per splicing event to %s", output_name
            )
            args.events_summary.write(
                "# {}\n".format(metadata_json.replace("\n", "\n# "))
            )
            (
                df_events
                # manually format probability columns
                .pipe(
                    lambda df: df.assign(
                        **{
                            col: df[col].apply(lambda x: f"{x:.3e}")
                            for col in df.columns
                            if "probability" in col
                        }
                    )
                )
                # numeric columns need at most 4 digits precision
                .round(4)
                # save result in TSV format
                .to_csv(args.events_summary, sep="\t", index=True)
            )
        if args.genes_summary:
            try:
                output_name = args.genes_summary.name
            except AttributeError:
                output_name = args.genes_summary
            log.info("Saving metadata and summary per gene to %s", output_name)
            args.genes_summary.write(
                "# {}\n".format(metadata_json.replace("\n", "\n# "))
            )
            (
                nm.PsiOutliers.summarize_df_genes(df_events)
                # manually format probability columns
                .pipe(
                    lambda df: df.assign(
                        **{
                            col: df[col].apply(lambda x: f"{x:.3e}")
                            for col in df.columns
                            if "probability" in col
                        }
                    )
                )
                # numeric columns need at most 4 digits precision
                .round(4)
                # save result in TSV format
                .to_csv(args.genes_summary, sep="\t", index=True)
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
