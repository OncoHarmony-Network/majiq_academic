"""
heterogen.py

Analysis of groups of independent experiments

Author: Joseph K Aicher
"""

import argparse
import json
import sys
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from dask.distributed import Client

import new_majiq as nm
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_characters_factory,
    check_nonnegative_factory,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = "Test differences in PSI for two groups of independent experiments"


def add_args(parser: argparse.ArgumentParser) -> None:
    comparison_req = parser.add_argument_group("Required specification of groups")
    StorePSICovPaths = StoreRequiredUniqueActionFactory()
    comparison_req.add_argument(
        "-grp1",
        type=ExistingResolvedPath,
        action=StorePSICovPaths,
        nargs="+",
        dest="psi1",
        required=True,
        help="Paths to PsiCoverage files for experiments in first group",
    )
    comparison_req.add_argument(
        "-grp2",
        type=ExistingResolvedPath,
        action=StorePSICovPaths,
        nargs="+",
        dest="psi2",
        required=True,
        help="Paths to PsiCoverage files for experiments in second group",
    )
    StoreGroupNames = StoreRequiredUniqueActionFactory()
    check_group_chars = check_characters_factory(
        nm.constants.ALLOWED_GROUP_NAME_CHARS, "alphanumeric or underscore characters"
    )
    comparison_req.add_argument(
        "-n",
        "--names",
        nargs=2,
        metavar=("NAME_GRP1", "NAME_GRP2"),
        required=True,
        action=StoreGroupNames,
        type=check_group_chars,
        help="The names that identify the groups being compared.",
    )

    parser.add_argument(
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        help="Threshold for group filters. This specifies the fraction"
        " (value < 1) or absolute number (value >= 1) of experiments that must"
        " pass individually in each group for an LSV to be quantified"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--stats",
        type=str,
        nargs="+",
        choices=set(nm.constants.STATS_AVAILABLE.keys()),
        default=nm.constants.DEFAULT_HET_USESTATS,
        help="Specify which stats to use for testing differences between groups"
        f" (available: {set(nm.constants.STATS_AVAILABLE.keys())})"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--splicegraph",
        type=ExistingResolvedPath,
        default=None,
        help="If specified, annotate quantifications with splicegraph information",
    )
    parser.add_argument(
        "--output-tsv",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path for output TSV file (default: stdout)",
    )
    parser.add_argument(
        "--population-quantiles",
        type=float,
        nargs="*",
        default=nm.constants.DEFAULT_HET_POPULATION_QUANTILES,
        help="Quantiles of PSI (besides median) per group to report"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--nthreads",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_QUANTIFY_NTHREADS,
        help="Number of threads used by Dask scheduler to quantify in chunks"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--memory-limit",
        type=str,
        default="auto",
        help="Memory limit to pass to dask cluster (default: %(default)s)",
    )

    psisample_settings = parser.add_argument_group(
        "Settings for repeated testing on posterior samples"
    )
    psisample_settings.add_argument(
        "--psisamples",
        type=check_nonnegative_factory(int, False),
        default=nm.constants.DEFAULT_HET_PSISAMPLES,
        help="Number of samples/tests to run and summarize. If set to zero,"
        " testing will only be done on posterior means (default: %(default)s).",
    )
    psisample_settings.add_argument(
        "--pvalue-quantiles",
        type=float,
        nargs="*",
        default=nm.constants.DEFAULT_HET_PVALUE_QUANTILES,
        help="Report these quantiles of pvalues on psisamples (default: %(default)s)",
    )
    return


def run(args: argparse.Namespace) -> None:
    population_quantiles = sorted(set(np.round(args.population_quantiles, 3)))
    pvalue_quantiles = sorted(set(np.round(args.pvalue_quantiles, 3)))
    use_stats = sorted(set(args.stats))

    log = get_logger()
    client = Client(
        n_workers=1,
        threads_per_worker=args.nthreads,
        dashboard_address=None,
        memory_limit=args.memory_limit,
    )
    log.info(client)
    nm.rng_resize(args.nthreads)
    metadata: Dict[str, Any] = dict()
    metadata["command"] = " ".join(sys.argv)
    metadata["version"] = nm.__version__
    metadata["population_quantiles"] = population_quantiles
    metadata["stats"] = use_stats
    metadata["pvalue_quantiles"] = pvalue_quantiles
    metadata["psisamples"] = args.psisamples

    # load psi1, psi2
    log.info("Loading input coverage files")
    psi1 = nm.PsiCoverage.from_zarr(args.psi1)
    psi2 = nm.PsiCoverage.from_zarr(args.psi2)
    if not psi1.events.equals(psi2.events):
        raise ValueError("Events from psi1 do not match events from psi2")
    group_sizes = {
        args.names[0]: psi1.num_prefixes,
        args.names[1]: psi2.num_prefixes,
    }
    metadata["group_sizes"] = group_sizes
    log.info(f"Analyzing {args.names[0]}({psi1}) vs {args.names[1]}({psi2})")

    # initialize list of data frames that will be concatenated (on columns)
    concat_df: List[pd.DataFrame] = list()

    # did we have a splicegraph to work with? add annotations...
    if args.splicegraph:
        log.info(f"Loading splicegraph from {args.splicegraph}")
        sg = nm.SpliceGraph.from_zarr(args.splicegraph)
        events = psi1.get_events(sg.introns, sg.junctions)
        concat_df.append(events.ec_dataframe)

    # dataset of quantifications I want
    log.info(f"Performing quantification on means and {args.psisamples} psisamples")
    heterogen = nm.Heterogen(
        psi1,
        psi2,
        min_experiments_f=args.min_experiments,
        name1=args.names[0],
        name2=args.names[1],
    )
    # which tests? which psi? which quantiles (pop and pval)? what psisamples?
    ds_quant = heterogen.dataset(
        population_quantiles=population_quantiles,
        pvalue_quantiles=pvalue_quantiles,
        psisamples=args.psisamples,
        use_stats=use_stats,
    ).load()

    log.info("Reshaping resulting quantifications to table")
    # pvalue columns that have dims (ec_idx, stats)
    df_pvalues = (
        ds_quant[
            [
                name
                for name, v in ds_quant.data_vars.items()
                if set(v.dims) == {"ec_idx", "stats"}
            ]
        ]
        .to_dataframe()
        .unstack(["stats"])
        .sort_index(axis=1)
    )
    df_pvalues.columns = [f"{stat}-{var}" for var, stat in df_pvalues.columns]
    concat_df.append(df_pvalues)
    if args.psisamples > 0 and pvalue_quantiles:
        # have dims (ec_idx, stats, pval_quantile)
        df_psisamples = (
            ds_quant[
                [
                    name
                    for name, v in ds_quant.data_vars.items()
                    if set(v.dims) == {"ec_idx", "stats", "pval_quantile"}
                ]
            ]
            .to_dataframe()
            .unstack(["stats", "pval_quantile"])
            .sort_index(axis=1)
        )
        df_psisamples.columns = [
            f"{stat}-{var}_{q:0.3f}" for var, stat, q in df_psisamples.columns
        ]
        concat_df.append(df_psisamples)
    # columns that have dims (grp, ec_idx) (always)
    df_medians = (
        ds_quant[
            [
                name
                for name, v in ds_quant.data_vars.items()
                if set(v.dims) == {"ec_idx", "grp"}
            ]
        ]
        .to_dataframe()
        .unstack(["grp"])
        .sort_index(axis=1)
    )
    df_medians.columns = [f"{grp}-{var}" for var, grp in df_medians.columns]
    concat_df.append(df_medians)
    if population_quantiles:
        df_pq = (
            ds_quant[
                [
                    name
                    for name, v in ds_quant.data_vars.items()
                    if set(v.dims) == {"ec_idx", "grp", "population_quantile"}
                ]
            ]
            .to_dataframe()
            .unstack(["grp", "population_quantile"])
            .sort_index(axis=1)
        )
        df_pq.columns = [f"{grp}-{var}_{q:0.3f}" for var, grp, q in df_pq.columns]
        concat_df.append(df_pq)

    log.info(f"Writing metadata to {args.output_tsv.name}")
    metadata_json = json.dumps(metadata, sort_keys=True, indent=4)
    args.output_tsv.write("# {}\n".format(metadata_json.replace("\n", "\n# ")))
    log.info(f"Writing table to {args.output_tsv.name}")
    (
        # concatenate columns together
        pd.concat(concat_df, axis=1, join="inner")
        # remove rows where no input passed
        .loc[ds_quant["passed"].values]
        # any column with pvalue in it needs to be manually formatted
        .pipe(
            lambda df: df.assign(
                **{
                    col: df[col].apply(lambda x: f"{x:.3e}")
                    for col in df.columns
                    if "pvalue" in col
                }
            )
        )
        # other numeric columns need at most 4 digits precision (psi)
        .round(4)
        # save result in TSV format
        .to_csv(args.output_tsv, sep="\t", index=not args.splicegraph)
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
