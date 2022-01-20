"""
deltapsi.py

Quantify differences between two groups of replicate experiments (MAJIQ
deltapsi)

Author: Joseph K Aicher
"""

import argparse
import json
import sys
from enum import Enum
from typing import Any, Dict, List, Optional

import new_majiq as nm
from new_majiq._run._majiq_args import (
    check_nonnegative_factory,
    quantify_comparison_args,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger

DESCRIPTION = "Quantify dPSI from two groups of replicate experiments"


class DPsiPriorType(Enum):
    DEFAULT_PRIOR = "default_prior"
    EMPIRICAL_PRIOR = "empirical_prior"


def add_args(parser: argparse.ArgumentParser) -> None:
    quantify_comparison_args(parser)

    quant_settings = parser.add_argument_group("deltapsi inference arguments")
    quant_settings.add_argument(
        "--psibins",
        metavar="B",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_QUANTIFY_PSIBINS,
        help="Number of bins for discretizing PSI distribution (twice this for"
        " dPSI) (default: %(default)s)",
    )
    quant_settings.add_argument(
        "--changing-threshold",
        metavar="CT",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_DPSI_CHANGING_THRESHOLD,
        help="Report posterior probability that abs(dPSI) greater than this"
        " threshold (default: %(default)s)",
    )
    quant_settings.add_argument(
        "--nonchanging-threshold",
        metavar="NCT",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_DPSI_NONCHANGING_THRESHOLD,
        help="Report posterior probability that abs(dPSI) less than this"
        " threshold (default: %(default)s)",
    )

    prior_args = parser.add_argument_group("deltapsi prior arguments")
    prior_type = prior_args.add_mutually_exclusive_group()
    prior_type.add_argument(
        "--empirical-prior",
        dest="prior_type",
        default=DPsiPriorType.EMPIRICAL_PRIOR,
        action="store_const",
        const=DPsiPriorType.EMPIRICAL_PRIOR,
        help="Update default prior using EM on empirical difference in PSI from"
        " high-confidence binary events (default: %(default)s)",
    )
    prior_type.add_argument(
        "--default-prior",
        dest="prior_type",
        default=DPsiPriorType.EMPIRICAL_PRIOR,
        action="store_const",
        const=DPsiPriorType.DEFAULT_PRIOR,
        help="Use default prior on deltapsi for quantification (default: %(default)s)",
    )
    prior_args.add_argument(
        "--prior-minreads",
        metavar="R",
        type=check_nonnegative_factory(float, True),
        default=nm.constants.DEFAULT_DPSI_PRIOR_MINREADS,
        help="Only use subset of connections with at least this many reads with"
        " min-experiments for empirical estimation of deltapsi prior"
        " (default: %(default)s)",
    )
    prior_args.add_argument(
        "--prior-minevents",
        metavar="N",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_DPSI_PRIOR_MINLSV,
        help="Require dPSI from at least this many binary splicing events to"
        " perform empirical estimation of deltapsi prior (default: %(default)s)",
    )
    prior_args.add_argument(
        "--prior-iter",
        metavar="N",
        type=check_nonnegative_factory(int, False),
        default=nm.constants.DEFAULT_DPSI_PRIOR_MAXITER,
        help="Number of EM updates for beta distribution parameters (will do 1"
        " more for mixture probabilities) (default: %(default)s)",
    )

    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
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
        args.names[0]: psi1.num_prefixes,
        args.names[1]: psi2.num_prefixes,
    }
    metadata["group_sizes"] = group_sizes
    log.info(f"Comparing {args.names[0]}({psi1}) vs {args.names[1]}({psi2})")

    prior = nm.DPsiPrior()
    if args.prior_type == DPsiPriorType.DEFAULT_PRIOR:
        log.info(f"Using default deltapsi prior {prior}")
    else:
        log.info(f"Starting from default deltapsi prior {prior}")
        prior = prior.empirical_update(
            psi1,
            psi2,
            minreads=args.prior_minreads,
            min_experiments_f=args.min_experiments,
            min_lsvs=args.prior_minevents,
            n_update_a=args.prior_iter,
            legacy=False,
            show_progress=args.show_progress,
        )
        log.info(f"Using deltapsi prior {prior}")

    # dataset of quantifications I want
    log.debug("Preparing quantification")
    deltapsi = nm.DeltaPsi(
        psi1,
        psi2,
        prior,
        psibins=args.psibins,
        min_experiments_f=args.min_experiments,
        name1=args.names[0],
        name2=args.names[1],
    )
    log.info("Computing quantifications to table")
    deltapsi_voila = deltapsi.dataset

    if args.output_voila:
        log.info("Saving quantifications for VOILA to %s", args.output_voila)
        deltapsi_voila.to_zarr(
            args.output_voila,
            ec_chunksize=args.chunksize,
            show_progress=args.show_progress,
        )
        deltapsi_voila = nm.DeltaPsiDataset.from_zarr(args.output_voila)

    sg: Optional[nm.SpliceGraph] = None
    if args.splicegraph:
        log.debug("Loading splicegraph from %s", args.splicegraph)
        sg = nm.SpliceGraph.from_zarr(args.splicegraph)

    df = deltapsi_voila.to_dataframe(
        sg=sg,
        changing_threshold=args.changing_threshold,
        nonchanging_threshold=args.nonchanging_threshold,
        # no need to show progress if most of the work already done in VOILA file
        show_progress=args.show_progress and args.output_voila is None,
    )

    try:
        output_name = args.output_tsv.name
    except AttributeError:
        output_name = args.output_tsv
    log.info("Writing metadata to %s", output_name)
    metadata_json = json.dumps(metadata, sort_keys=True, indent=4)
    args.output_tsv.write("# {}\n".format(metadata_json.replace("\n", "\n# ")))
    log.info("Writing table to %s", output_name)
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
        # other numeric columns need at most 4 digits precision
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
