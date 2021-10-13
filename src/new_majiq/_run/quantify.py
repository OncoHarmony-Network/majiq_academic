"""
quantify.py

Quantify input PsiCoverage files

Author: Joseph K Aicher
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import xarray as xr

import new_majiq as nm
from new_majiq._run._majiq_args import check_nonnegative_factory
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = "Quantify PSI from PsiCoverage files"


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "psicov",
        type=Path,
        nargs="+",
        help="Paths to PsiCoverage files. All must have been generated using same splicegraph",
    )
    parser.add_argument(
        "--splicegraph",
        type=Path,
        default=None,
        help="If specified, annotate quantifications with splicegraph information",
    )
    parser.add_argument(
        "--quantiles",
        type=float,
        nargs="+",
        default=None,
        help="If specified, calculate/report PSI posterior quantiles",
    )
    parser.add_argument(
        "--min-experiments",
        type=check_nonnegative_factory(float, True),
        default=None,
        help="If specified, treat samples as replicates and quantify combined coverage for events that passed in at least min experiments (proportion of experiments if < 1)",
    )
    parser.add_argument(
        "--output-tsv",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path for output TSV file (default: stdout)",
    )
    return


def run(args: argparse.Namespace) -> None:
    if not all(p.exists() for p in args.psicov):
        missing = sorted(p for p in args.psicov if not p.exists())
        raise ValueError(f"Unable to find input coverage ({missing = })")
    if args.splicegraph and not args.splicegraph.exists():
        raise ValueError(f"Unable to find input splicegraph {args.splicegraph}")

    log = get_logger()
    metadata: Dict[str, Any] = dict()
    metadata["command"] = " ".join(sys.argv)
    metadata["version"] = nm.__version__

    log.info(f"Joining {len(args.psicov)} input PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    metadata["n_experiments"] = psicov.df.sizes["prefix"]
    log.info(f"Using coverage from {metadata['n_experiments']} input experiments")
    if args.min_experiments is not None:
        log.info(f"Aggregating coverage using min-experiments = {args.min_experiments}")
        psicov = psicov.sum("aggregate", min_experiments_f=args.min_experiments)

    # initialize list of data frames that will be concatenated (on columns)
    concat_df: List[pd.DataFrame] = list()

    # did we have a splicegraph to work with? add annotations...
    if args.splicegraph:
        log.info(f"Loading splicegraph from {args.splicegraph}")
        sg = nm.SpliceGraph.from_zarr(args.splicegraph)
        events = psicov.get_events(sg.introns, sg.junctions)
        concat_df.append(events.ec_dataframe)

    # get desired quantifications
    quantify_vars = {
        "any_passed": psicov.event_passed.any("prefix"),
        "raw_psi_mean": psicov.raw_posterior_mean,
        "raw_psi_std": np.sqrt(psicov.raw_posterior_variance),
        "bootstrap_psi_mean": psicov.bootstrap_posterior_mean,
        "bootstrap_psi_std": np.sqrt(psicov.bootstrap_posterior_variance),
        "raw_coverage": psicov.raw_coverage,
    }
    if args.quantiles:
        quantiles = sorted(set(np.round(args.quantiles, 3)))
        log.info("Will compute the following posterior quantiles: {quantiles}")
        quantify_vars["approx_psi_quantile"] = psicov.approximate_quantile(quantiles)
    log.info("Performing quantification")
    ds_quant = xr.Dataset(quantify_vars).reset_coords(drop=True).load()  # type: ignore[arg-type]

    log.info("Reshaping resulting quantifications to table")
    # all quantifications but quantiles
    df_quant = (
        ds_quant.drop_dims("quantiles", errors="ignore")
        .drop_vars("any_passed")
        .to_dataframe()
        .unstack("prefix")
        .reorder_levels([1, 0], axis=1)
        .sort_index(axis=1)
    )
    df_quant.columns = [
        f"{prefix} {var}" if ds_quant.sizes["prefix"] > 1 else var
        for (prefix, var) in df_quant.columns.values
    ]
    concat_df.append(df_quant)
    if args.quantiles:
        df_quantiles = (
            ds_quant["approx_psi_quantile"]
            .to_series()
            .unstack(["prefix", "quantiles"])
            .sort_index(axis=1)
        )
        df_quantiles.columns = [
            f"{prefix} approx_psi_quantile_{q:0.3f}"
            if ds_quant.sizes["prefix"] > 1
            else f"approx_psi_quantile_{q:0.3f}"
            for (prefix, q) in df_quantiles.columns.values
        ]
        concat_df.append(df_quantiles)
    log.info(f"Writing metadata to {args.output_tsv.name}")
    metadata_json = json.dumps(metadata, sort_keys=True, indent=4)
    args.output_tsv.write("# {}\n".format(metadata_json.replace("\n", "\n# ")))
    log.info(f"Writing table to {args.output_tsv.name}")
    (
        # concatenate columns together
        pd.concat(concat_df, axis=1, join="inner")
        # remove rows where no input passed
        .loc[ds_quant["any_passed"].values]
        # numeric columns need at most 4 digits precision
        .round(4)
        # save result in TSV format
        .to_csv(args.output_tsv, sep="\t", index=False)
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
