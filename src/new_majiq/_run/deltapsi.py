"""
deltapsi.py

Quantify differences between two groups of replicate experiments (MAJIQ
deltapsi)

Author: Joseph K Aicher
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

import new_majiq as nm
from new_majiq._run._majiq_args import (
    StoreRequiredUniqueActionFactory,
    check_characters_factory,
    check_nonnegative_factory,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.DeltaPsi import DeltaPsi
from new_majiq.DPsiPrior import DPsiPrior
from new_majiq.logger import get_logger

DESCRIPTION = "Quantify dPSI from two groups of replicate experiments"


def add_args(parser: argparse.ArgumentParser) -> None:
    comparison_req = parser.add_argument_group("Required specification of groups")
    comparison_req.add_argument(
        "-psi1",
        type=Path,
        nargs="+",
        dest="psi1",
        required=True,
        help="Paths to PsiCoverage files for experiments contributing to psi1",
    )
    comparison_req.add_argument(
        "-psi2",
        type=Path,
        nargs="+",
        dest="psi2",
        required=True,
        help="Paths to PsiCoverage files for experiments contributing to psi2",
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
        "--splicegraph",
        type=Path,
        default=None,
        help="If specified, annotate quantifications with splicegraph information",
    )
    parser.add_argument(
        "--output-tsv",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Path for output TSV file (default: stdout)",
    )

    quant_settings = parser.add_argument_group("Quantification settings")
    quant_settings.add_argument(
        "--psibins",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_QUANTIFY_PSIBINS,
        help="Number of bins for discretizing PSI distribution (twice this for"
        " dPSI) (default: %(default)s)",
    )
    quant_settings.add_argument(
        "--changing-threshold",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_DPSI_CHANGING_THRESHOLD,
        help="Report posterior probability that abs(dPSI) greater than this"
        " threshold (default: %(default)s)",
    )
    quant_settings.add_argument(
        "--nonchanging-threshold",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_DPSI_NONCHANGING_THRESHOLD,
        help="Report posterior probability that abs(dPSI) less than this"
        " threshold (default: %(default)s)",
    )
    quant_settings.add_argument(
        "--use-posterior",
        type=str,
        default="smooth",
        choices=["smooth", "legacy", "both"],
        help="Perform deltapsi inference with new (smooth) vs legacy approach"
        " (or both) (default: %(default)s",
    )

    prior_args = parser.add_argument_group("Configure prior on deltapsi")
    default_vs_empirical = prior_args.add_mutually_exclusive_group()
    default_vs_empirical.add_argument(
        "--empirical-prior",
        action="store_true",
        dest="empirical_prior",
        default=True,
        help="Update default prior using empirical difference in PSI from"
        " high-confidence binary events (default empirical=%(default)s)",
    )
    default_vs_empirical.add_argument(
        "--default-prior",
        action="store_false",
        dest="empirical_prior",
        default=True,
        help="Use default prior on deltapsi for quantification",
    )
    prior_args.add_argument(
        "--prior-minreads",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_DPSI_PRIOR_MINREADS,
        help="Only use subset of connections with at least this many reads with"
        " min-experiments for empirical estimation of deltapsi prior"
        " (default: %(default)s)",
    )
    prior_args.add_argument(
        "--prior-minevents",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_DPSI_PRIOR_MINLSV,
        help="Require dPSI from at least this many binary splicing events to"
        " perform empirical estimation of deltapsi prior (default: %(default)s)",
    )
    prior_args.add_argument(
        "--prior-iter",
        type=check_nonnegative_factory(int, False),
        default=nm.constants.DEFAULT_DPSI_PRIOR_MAXITER,
        help="Number of EM updates for beta distribution parameters (will do 1"
        " more for mixture probabilities) (default: %(default)s)",
    )

    return


def run(args: argparse.Namespace) -> None:
    if not all(p.exists() for p in args.psi1):
        missing = sorted(p for p in args.psi1 if not p.exists())
        raise ValueError(f"Unable to find input coverage ({missing = })")
    if not all(p.exists() for p in args.psi2):
        missing = sorted(p for p in args.psi2 if not p.exists())
        raise ValueError(f"Unable to find input coverage ({missing = })")
    if args.splicegraph and not args.splicegraph.exists():
        raise ValueError(f"Unable to find input splicegraph {args.splicegraph}")

    log = get_logger()
    metadata: Dict[str, Any] = dict()
    metadata["command"] = " ".join(sys.argv)
    metadata["version"] = nm.__version__
    metadata["changing_threshold"] = args.changing_threshold
    metadata["nonchanging_threshold"] = args.nonchanging_threshold
    metadata["psibins"] = args.psibins

    # load psi1, psi2
    log.info("Loading input coverage files")
    psi1 = nm.PsiCoverage.from_zarr(args.psi1)
    psi2 = nm.PsiCoverage.from_zarr(args.psi2)
    if not psi1.events.equals(psi2.events):
        raise ValueError("Events from psi1 do not match events from psi2")
    group_sizes = {
        args.names[0]: psi1.df.sizes["prefix"],
        args.names[1]: psi2.df.sizes["prefix"],
    }
    metadata["group_sizes"] = group_sizes
    log.info(f"Number of experiments per group: {group_sizes}")

    # initialize list of data frames that will be concatenated (on columns)
    concat_df: List[pd.DataFrame] = list()

    # did we have a splicegraph to work with? add annotations...
    if args.splicegraph:
        log.info(f"Loading splicegraph from {args.splicegraph}")
        sg = nm.SpliceGraph.from_zarr(args.splicegraph)
        events = psi1.get_events(sg.introns, sg.junctions)
        concat_df.append(events.ec_dataframe)

    prior = DPsiPrior()
    if args.empirical_prior:
        log.info(f"Starting from default deltapsi prior {prior}")
        prior = prior.empirical_update(
            psi1,
            psi2,
            minreads=args.prior_minreads,
            min_experiments_f=args.min_experiments,
            min_lsvs=args.prior_minevents,
            n_update_a=args.prior_iter,
        )
        log.info(f"Using deltapsi prior {prior}")
    else:
        log.info(f"Using default deltapsi prior {prior}")

    # dataset of quantifications I want
    log.info("Performing quantification")
    deltapsi = DeltaPsi(
        psi1,
        psi2,
        prior,
        psibins=args.psibins,
        min_experiments_f=args.min_experiments,
        name1=args.names[0],
        name2=args.names[1],
    )
    ds_quant = deltapsi.dataset(
        changing_threshold=args.changing_threshold,
        nonchanging_threshold=args.nonchanging_threshold,
        use_posterior=args.use_posterior,
    ).load()

    log.info("Reshaping resulting quantifications to table")
    concat_df.append(ds_quant.drop_vars("passed").to_dataframe())
    log.info(f"Writing metadata to {args.output_tsv.name}")
    metadata_json = json.dumps(metadata, sort_keys=True, indent=4)
    args.output_tsv.write("# {}\n".format(metadata_json.replace("\n", "\n# ")))
    log.info(f"Writing table to {args.output_tsv.name}")
    (
        pd.concat(concat_df, axis=1, join="inner")
        # remove rows where no input passed
        .loc[ds_quant["passed"].values]
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
        # other numeric columns need at most 4 digits precision
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
