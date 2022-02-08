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
    check_range_factory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand

DESCRIPTION = "Identify outliers from PSI cases (N ~ 1) vs PSI controls (N >> 1)"


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--pass",
        dest="pass_args_list",
        metavar=("splicegraph", "controls", "cases"),
        nargs=3,
        type=ExistingResolvedPath,
        action="append",
        required=True,
        help="Specify an analysis pass of SpliceGraph, PsiControlsSummary,"
        " and PsiCoverage for cases for outliers analysis."
        " Use this flag multiple times to merge analysis over multiple"
        " build/quantification passes."
        " Final results will be relative to the last specified pass."
        " Each pass' PsiControlsSummary and PsiCoverage must share events"
        " defined on the pass' SpliceGraph."
        " Each SpliceGraph must share genes."
        " Controls are required to share the same experiments."
        " The cases from each pass must include the same experiments as the"
        " cases from the last specified pass.",
    )
    annotated = parser.add_argument_group("annotated features arguments")
    annotated_ex = annotated.add_mutually_exclusive_group(required=False)
    annotated_ex.add_argument(
        "--annotated",
        dest="annotated",
        metavar="splicegraph",
        type=ExistingResolvedPath,
        default=None,
        help="Identify novel events/exons/introns/junctions relative to this"
        " splicegraph."
        " Default: use the splicegraph from the second-last analysis pass"
        " if more than one pass was specified (otherwise '--no-annotated').",
    )
    annotated_ex.add_argument(
        "--no-annotated",
        dest="annotated",
        action="store_false",
        default=None,
        help="Identify novel exons/introns/junctions using definitions from"
        " last specified SpliceGraph."
        " Do not identify novel features relative to another splicegraph."
        " Default: use the splicegraph from the second-last analysis pass"
        " if more than one pass was specified",
    )
    tsv = parser.add_argument_group("output TSV locations arguments")
    tsv.add_argument(
        "--output-tsv",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path for output TSV file (default: stdout)",
    )
    tsv.add_argument(
        "--events-summary",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=None,
        help="Path for summary per splicing event (default: None)",
    )
    tsv.add_argument(
        "--genes-summary",
        metavar="TSV",
        type=argparse.FileType("w"),
        default=None,
        help="Path for summary per splicing gene (default: None)",
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
        help="Select single threshold on controls/case quantiles in both"
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

    # args.pass_args_list: List[Tuple[sg, controls, cases]]

    log.debug("Loading input splicegraphs")
    use_genes = nm.Genes.from_zarr(args.pass_args_list[-1][0])
    sg_list: List[nm.SpliceGraph] = [
        nm.SpliceGraph.from_zarr(sg_path, genes=use_genes)
        for sg_path, _, _ in args.pass_args_list
    ]
    log.info(
        "Outliers analysis on %d passes with splicegraphs %s", len(sg_list), sg_list
    )

    log.debug("Loading input controls")
    controls_list: List[nm.PsiControlsSummary] = [
        nm.PsiControlsSummary.from_zarr(controls_path)
        for _, controls_path, _ in args.pass_args_list
    ]
    if len(controls_list) > 1:
        log.debug("Verifying that controls were defined with the same experiments")
        for pass_ct, compare_controls in enumerate(controls_list[:-1], 1):
            if missing := set(compare_controls.prefixes).symmetric_difference(
                set(controls_list[-1].prefixes)
            ):
                raise ValueError(
                    f"Pass {pass_ct} controls were not defined with same"
                    f" experiments as pass {len(controls_list)}"
                    f" (not shared: {missing})."
                )
        if args.select_alpha is None:
            log.debug(
                "Verifying that controls were defined with the same values of alpha"
            )
            for pass_ct, compare_controls in enumerate(controls_list[:-1], 1):
                if set(compare_controls.alpha.values) != set(
                    controls_list[-1].alpha.values
                ):
                    raise ValueError(
                        "Input controls are defined with different values of alpha"
                        f" (pass {pass_ct}: {set(compare_controls.alpha.values)},"
                        f" pass {len(controls_list)}:"
                        f" {set(controls_list[-1].alpha.values)})."
                    )
    if args.select_alpha is not None:
        for pass_ct, compare_controls in enumerate(controls_list, 1):
            if args.select_alpha not in compare_controls.alpha.values:
                raise ValueError(
                    f"Selected value of alpha ({args.select_alpha}) is not"
                    f" defined for pass {pass_ct} controls from"
                    f" {args.pass_args_list[pass_ct][1]}"
                    f" (defined: {compare_controls.alpha.values})."
                )
    metadata["n_controls"] = controls_list[-1].num_prefixes
    log.info("Quantifying potential outliers relative to controls %s", controls_list)

    log.debug("Loading input cases")
    cases_list: List[nm.PsiCoverage] = [
        nm.PsiCoverage.from_zarr(cases_path) for _, _, cases_path in args.pass_args_list
    ]
    if len(cases_list) > 1:
        for pass_ct, compare_cases in enumerate(cases_list[:-1], 1):
            if missing := set(cases_list[-1].prefixes).difference(
                compare_cases.prefixes
            ):
                raise ValueError(
                    f"Pass {pass_ct} cases does not have coverage for same"
                    f" experiments as pass {len(cases_list)} ({missing = })"
                )
            else:
                # taking subset of experiments in previous passes
                cases_list[pass_ct - 1] = compare_cases[cases_list[-1].prefixes]
    metadata["n_cases"] = cases_list[-1].num_prefixes
    log.info("Quantifying potential outliers from cases %s", cases_list)

    log.debug("Defining how novel features will be identified")
    annotated: Optional[nm.SpliceGraph]
    if args.annotated is None and len(sg_list) > 1:
        # second last splicegraph as default
        args.annotated = args.pass_args_list[-2][0]
    if args.annotated:
        try:
            annotated = sg_list[
                [sg_path for sg_path, _, _ in args.pass_args_list].index(args.annotated)
            ]
        except ValueError:  # args.annotated not from one of the input passes
            annotated = nm.SpliceGraph.from_zarr(args.annotated, genes=use_genes)
        log.info(
            "Novel features will be defined relative to %s (%s)",
            args.annotated,
            annotated,
        )
    else:  # explicitly used --no-annotated or only one pass specified
        annotated = None
        log.info("Novel features will be defined as in %s", args.pass_args_list[-1][0])

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
    if len(df_list) > 1:
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
        df = events.merge_dataframes(df_list, events_list, annotated=annotated)
    else:
        df = (
            controls_list[0]
            .get_events(sg_list[0].introns, sg_list[0].junctions)
            .ec_dataframe(annotated=annotated)
            .join(df_list[0], how="inner", on="ec_idx")
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
        .to_csv(args.output_tsv, sep="\t", index=False)
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
