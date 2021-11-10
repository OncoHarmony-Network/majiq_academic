"""
quantify.py

Quantify input PsiCoverage files

Author: Joseph K Aicher
"""

import argparse
import json
import sys
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from dask.distributed import progress

import new_majiq as nm
from new_majiq._run._majiq_args import quantify_nocomparison_args, resources_args
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = "Quantify PSI from PsiCoverage files"


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--quantiles",
        metavar="Q",
        type=float,
        nargs="*",
        default=list(),
        help="If specified, calculate/report PSI posterior quantiles",
    )
    quantify_nocomparison_args(parser)
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
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

    log.info("Performing quantification")
    ds_quant = psicov.dataset(quantiles=sorted(set(np.round(args.quantiles, 3))))
    if args.show_progress:
        ds_quant = ds_quant.persist()
        progress(*(x.data for x in ds_quant.variables.values() if x.chunks))
    ds_quant = ds_quant.load()

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
            ds_quant[[name for name, v in ds_quant.items() if "quantiles" in v.dims]]
            .to_dataframe()
            .unstack(["prefix", "quantiles"])
            .sort_index(axis=1)
        )
        df_quantiles.columns = [
            f"{prefix} {var}_{q:0.3f}"
            if ds_quant.sizes["prefix"] > 1
            else f"{var}_{q:0.3f}"
            for (var, prefix, q) in df_quantiles.columns.values
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
