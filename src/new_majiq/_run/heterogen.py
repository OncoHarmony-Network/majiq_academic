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

import new_majiq as nm
from new_majiq._run._majiq_args import (
    check_nonnegative_factory,
    quantify_comparison_args,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = "Test differences in PSI for two groups of independent experiments"


def add_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--population-quantiles",
        metavar="Q",
        type=float,
        nargs="*",
        default=nm.constants.DEFAULT_HET_POPULATION_QUANTILES,
        help="Quantiles of PSI (besides median) per group to report"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--stats",
        metavar="S",
        type=str,
        nargs="+",
        choices=set(nm.constants.STATS_AVAILABLE.keys()),
        default=nm.constants.DEFAULT_HET_USESTATS,
        help="Specify which stats to use for testing differences between groups"
        f" (available: {set(nm.constants.STATS_AVAILABLE.keys())})"
        " (default: %(default)s)",
    )
    quantify_comparison_args(parser)

    psisample_settings = parser.add_argument_group(
        "psisamples (posterior samples) testing arguments"
    )
    psisample_settings.add_argument(
        "--psisamples",
        metavar="J",
        type=check_nonnegative_factory(int, False),
        default=nm.constants.DEFAULT_HET_PSISAMPLES,
        help="Number of samples/tests to run and summarize. If set to zero,"
        " testing will only be done on posterior means (default: %(default)s).",
    )
    psisample_settings.add_argument(
        "--pvalue-quantiles",
        metavar="Q",
        type=float,
        nargs="*",
        default=nm.constants.DEFAULT_HET_PVALUE_QUANTILES,
        help="Report these quantiles of pvalues on psisamples (default: %(default)s)",
    )

    distribution_settings = parser.add_argument_group(
        "posteriors for testing arguments"
    )
    distribution_raw = distribution_settings.add_mutually_exclusive_group()
    distribution_raw.add_argument(
        "--test-raw",
        dest="raw_stats",
        action="store_true",
        default=nm.constants.DEFAULT_HET_RAWSTATS,
        help="Test for differences on raw_psi_mean (default: %(default)s)",
    )
    distribution_raw.add_argument(
        "--ignore-raw",
        dest="raw_stats",
        action="store_false",
        default=nm.constants.DEFAULT_HET_RAWSTATS,
        help="Do not test for differences on raw_psi_mean (default: %(default)s)",
    )
    distribution_approximate = distribution_settings.add_mutually_exclusive_group()
    distribution_approximate.add_argument(
        "--test-bootstrap",
        dest="approximate_stats",
        action="store_true",
        default=nm.constants.DEFAULT_HET_APPROXSTATS,
        help="Test for differences on bootstrap_psi_mean and psisamples"
        " samples from smooth approximation of mixture of bootstrap posteriors"
        " (default: %(default)s)",
    )
    distribution_approximate.add_argument(
        "--ignore-bootstrap",
        dest="approximate_stats",
        action="store_false",
        default=nm.constants.DEFAULT_HET_APPROXSTATS,
        help="Do not test for differences with approximation of bootstrap"
        " posteriors (default: %(default)s)",
    )
    resources_args(parser, use_dask=True)
    parser.add_argument(
        "--psibins",
        metavar="B",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_QUANTIFY_PSIBINS,
        help="Number of bins for discretizing PSI distribution for visualization"
        " (default: %(default)s)",
    )
    return


def run(args: argparse.Namespace) -> None:
    population_quantiles = sorted(set(np.round(args.population_quantiles, 3)))
    pvalue_quantiles = sorted(set(np.round(args.pvalue_quantiles, 3)))
    use_stats = sorted(set(args.stats))

    log = get_logger()
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
    log.info(f"Comparing {args.names[0]}({psi1}) vs {args.names[1]}({psi2})")

    # dataset of quantifications I want
    log.info(f"Performing quantification on means and {args.psisamples} psisamples")
    heterogen = nm.Heterogen(
        psi1,
        psi2,
        min_experiments_f=args.min_experiments,
        name1=args.names[0],
        name2=args.names[1],
    )
    heterogen_voila = heterogen.dataset(
        pvalue_quantiles=pvalue_quantiles,
        use_stats=use_stats,
        psisamples=args.psisamples,
        psibins=nm.constants.DEFAULT_QUANTIFY_PSIBINS,
    )

    if args.output_voila:
        log.info("Saving quantifications for VOILA to %s", args.output_voila)
        heterogen_voila.to_zarr(
            args.output_voila,
            ec_chunksize=args.chunksize,
            show_progress=args.show_progress,
        )
        heterogen_voila = nm.HeterogenDataset.from_zarr(args.output_voila)

    sg: Optional[nm.SpliceGraph] = None
    if args.splicegraph:
        log.debug("Loading splicegraph from %s", args.splicegraph)
        sg = nm.SpliceGraph.from_zarr(args.splicegraph)

    log.info("Summarizing population quantiles per group")
    df = heterogen_voila.to_dataframe(
        sg=sg,
        population_quantiles=population_quantiles,
        show_progress=args.show_progress,
    )

    try:
        output_name = args.output_tsv.name
    except AttributeError:
        output_name = args.output_tsv
    log.info(f"Writing metadata to {output_name}")
    metadata_json = json.dumps(metadata, sort_keys=True, indent=4)
    args.output_tsv.write("# {}\n".format(metadata_json.replace("\n", "\n# ")))
    log.info(f"Writing table to {output_name}")
    (
        df
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
