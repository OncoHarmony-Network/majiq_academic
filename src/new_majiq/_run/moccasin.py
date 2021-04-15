"""
moccasin.py

Perform batch correction using the moccasin algorithm

Author: Joseph K Aicher
"""

import argparse
import numpy as np
import pandas as pd
import xarray as xr
import new_majiq as nm
import new_majiq.constants as constants
import new_majiq.moccasin as mc

from new_majiq.experiments import bam_experiment_name
from new_majiq._run._majiq_args import check_nonnegative_factory
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger
from pathlib import Path
from typing import (
    List,
    Optional,
    Union,
)


def _open_mf_with_prefix(
    paths: List[Path], group: Optional[str] = None
) -> xr.DataArray:
    """wrap xr.open_mfdataset to load zarr, add prefix from filenames"""
    if len(set(bam_experiment_name(x) for x in paths)) < len(paths):
        raise ValueError("paths have non-unique prefixes")
    return xr.open_mfdataset(
        paths,
        concat_dim="prefix",
        engine="zarr",
        join="exact",
        group=group,
        preprocess=lambda x: x.expand_dims(
            prefix=[bam_experiment_name(x.encoding["source"])]
        ),
    )


def _args_quantifiable_thresholds(parser: argparse.ArgumentParser) -> None:
    """add argument group for quantifier thresholds"""
    thresholds = parser.add_argument_group("Quantifiability thresholds")
    thresholds.add_argument(
        "--minreads",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINREADS,
        help="Minimum readrate per experiment to pass a connection"
        " (default: %(default)s)",
    )
    thresholds.add_argument(
        "--minbins",
        type=check_nonnegative_factory(float, True),
        default=constants.DEFAULT_QUANTIFY_MINBINS,
        help="Minimum number of nonzero bins to pass a connection"
        " (default: %(default)s).",
    )
    return


def _get_quantifiable_thresholds(args: argparse.Namespace) -> nm.QuantifierThresholds:
    """get QuantifierThresholds from parsed arguments"""
    return nm.QuantifierThresholds(
        minreads=args.minreads, minbins=args.minbins, min_experiments_f=1
    )


def _args_factors(parser: argparse.ArgumentParser) -> None:
    """add arguments to specify factors to load"""
    factors_ex = parser.add_mutually_exclusive_group(required=True)
    factors_ex.add_argument(
        "--intercept-only",
        dest="intercept_only",
        action="store_true",
        default=None,
        help="Only use non-confounding intercept term"
        " (only makes sense if augmenting with discovered unknown factors)",
    )
    factors_ex.add_argument(
        "--factors-tsv",
        dest="factors_tsv",
        type=Path,
        default=None,
        help="Path to TSV with model matrix for all prefixes being processed"
        " (required column: prefix)",
    )
    factors_ex.add_argument(
        "--factors-zarr",
        dest="factors_zarr",
        type=Path,
        nargs="+",
        default=None,
        help="Paths to factors per experiment (output by infer-factors)",
    )
    parser.add_argument(
        "--confounding",
        type=str,
        nargs="+",
        default=list(),
        help="Names of confounding variables, only used with --factors-tsv",
    )
    return


def _get_factors(
    prefix: Union[str, List[str]], args: argparse.Namespace
) -> xr.DataArray:
    """get factors matrix for specified prefix(es)"""
    do_squeeze = False  # do we expect a prefix dimension?
    if isinstance(prefix, str):
        prefix = [prefix]  # treat as list of 1
        do_squeeze = True  # squeeze it back out, though
    factors: xr.DataArray
    if args.intercept_only:
        factors = xr.DataArray(
            np.ones((len(prefix), 1)),
            {
                "prefix": prefix,
                "factor": ["intercept"],
                "confounding": ("factor", [False]),
            },
            dims=["prefix", "factor"],
        )
    elif args.factors_tsv:
        df = (
            pd.read_csv(args.factors_tsv, sep="\t")
            # get string prefix index in sorted order
            .astype({"prefix": str}).set_index("prefix", verify_integrity=True)
            # select desired prefixes
            .loc[prefix]
            # convert all other values to float32
            .astype(np.float32)
        )
        if (missing := set(x for x in args.confounding if x not in df.columns)) :
            raise ValueError(
                f"Not all specified confounders were found from TSV ({missing = })"
            )
        factors = xr.DataArray(
            df.values,
            {
                "prefix": df.index,
                "factor": df.columns,
                "confounding": ("factor", [x in args.confounding for x in df.columns]),
            },
            dims=["prefix", "factor"],
        )
    else:  # args.factors_zarr
        prefix_set = set(prefix)
        matched_zarr = [
            x for x in set(args.factors_zarr) if bam_experiment_name(x) in prefix_set
        ]
        if (missing := prefix_set - {bam_experiment_name(x) for x in matched_zarr}) :
            raise ValueError(f"Factors for some prefixes missing ({missing = })")
        factors = _open_mf_with_prefix(matched_zarr).factors
    # do the squeeze, if desired
    if do_squeeze:
        factors = factors.squeeze("prefix")
    return factors


def _args_ruv_parameters(parser: argparse.ArgumentParser) -> None:
    ruv = parser.add_argument_group("Parameters for modeling unknown confounders")
    ruv.add_argument(
        "--ruv-max-new-factors",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_MOCCASIN_RUV_MAX_FACTORS,
        help="Maximum number of new factors the model will add",
    )
    ruv.add_argument(
        "--ruv-max-events",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_MOCCASIN_RUV_MAX_EVENTS,
        help="Maximum number of events to model residuals for",
    )
    return


def args_tmpfile(parser: argparse.ArgumentParser) -> None:
    """arguments for creating tmpfile for individual majiq file"""
    parser.add_argument(
        "coverage",
        type=Path,
        help="Path to LSV coverage (i.e. majiq file) to temporarily"
        " reparameterize for moccasin computation",
    )
    parser.add_argument(
        "tmpfile",
        type=Path,
        help="Path for temporary reparameterization of the majiq file",
    )
    _args_quantifiable_thresholds(parser)
    return


def args_factors_model(parser: argparse.ArgumentParser) -> None:
    """arguments for creating model of unknown factors"""
    _args_factors(parser)  # get factors information
    parser.add_argument(
        "factors_model", type=Path, help="Path for output model for unknown confounders"
    )
    parser.add_argument(
        "tmpfiles",
        metavar="tmpfile",
        type=Path,
        nargs="+",
        help="Paths for input tmpfiles with coverage for unknown confounders model",
    )
    _args_ruv_parameters(parser)
    return


def args_factors_infer(parser: argparse.ArgumentParser) -> None:
    """arguments for inferring factors using factors model"""
    _args_factors(parser)
    parser.add_argument(
        "factors_model",
        type=Path,
        help="Path for input model of unknown confounders",
    )
    parser.add_argument(
        "tmpfile",
        type=Path,
        help="Path to tmpfile for which known/unknown factors will be saved",
    )
    parser.add_argument(
        "output",
        type=Path,
        help="Path for output factors (known and unknown)",
    )
    return


def args_coverage_model(parser: argparse.ArgumentParser) -> None:
    """arguments for modeling coverage"""
    parser.add_argument(
        "coverage_model",
        type=Path,
        help="Output path for coverage model",
    )
    _args_factors(parser)
    parser.add_argument(
        "tmpfiles",
        metavar="tmpfile",
        nargs="+",
        type=Path,
        help="Paths for input tmpfiles with coverage for coverage model",
    )
    return


def args_coverage_infer(parser: argparse.ArgumentParser) -> None:
    """arguments for getting corrected lsv coverage"""
    parser.add_argument(
        "original_coverage",
        type=Path,
        help="Path to original, uncorrected LSV coverage",
    )
    parser.add_argument(
        "tmpfile",
        type=Path,
        help="Path to tmpfile corresponding to original_coverage",
    )
    _args_factors(parser)
    parser.add_argument(
        "coverage_model",
        type=Path,
        help="Input path for coverage model",
    )
    parser.add_argument(
        "corrected_coverage",
        type=Path,
        help="Path for corrected LSV coverage",
    )
    return


def run_tmpfile(args: argparse.Namespace) -> None:
    log = get_logger()
    if bam_experiment_name(args.coverage) != bam_experiment_name(args.tmpfile):
        raise ValueError(
            "Coverage file and temporary reparameterization have different prefixes"
        )
    log.info(f"Loading and reparameterizing coverage in {args.coverage.resolve()}")
    df_tmp = mc.tmp_psi_for_moccasin(args.coverage, _get_quantifiable_thresholds(args))
    log.info(f"Saving reparameterization to {args.tmpfile.resolve()}")
    df_tmp.to_zarr(args.tmpfile, mode="w")
    return


def run_factors_model(args: argparse.Namespace) -> None:
    """model unknown factors using input coverage"""
    prefix = [bam_experiment_name(x) for x in args.tmpfiles]
    if len(prefix) > len(set(prefix)):
        raise ValueError("Passed tmpfiles have duplicate prefixes")
    log = get_logger()
    log.info("Setting up model matrix of known factors")
    factors = _get_factors(prefix, args)
    log.info(f"Opening coverage from {len(args.tmpfiles)} tmpfiles")
    coverage = _open_mf_with_prefix(args.tmpfiles)
    log.info("Learning model for unknown confounding factors")
    model = mc.ModelUnknownConfounders.train(
        coverage.psi,
        coverage._offsets.values,
        factors,
        max_new_factors=args.ruv_max_new_factors,
        max_events=args.ruv_max_events,
    )
    log.info(
        f"Learned {model.num_factors} factors explaining"
        f" {model.explained_variance.sum().values[()]} of residual variance"
    )
    log.info(f"Saving model to {args.factors_model.resolve()}")
    model.to_zarr(args.factors_model)
    return


def run_factors_infer(args: argparse.Namespace):
    """compute unknown factors using input coverage"""
    prefix = bam_experiment_name(args.tmpfile)
    if prefix != bam_experiment_name(args.output):
        raise ValueError("tmpfile and output factors file have different prefixes")
    log = get_logger()
    log.info("Setting up known factors")
    factors = _get_factors(prefix, args).load()
    log.info(f"Opening coverage from {args.tmpfile.resolve()}")
    coverage = xr.open_zarr(args.tmpfile)
    log.info(f"Loading model for unknown factors from {args.factors_model.resolve()}")
    model = mc.ModelUnknownConfounders.from_zarr(args.factors_model)
    log.info("Solving for unknown confounders")
    unknown_factors = model.predict(coverage.psi, factors).load()
    log.info(f"Saving combined factors to {args.output.resolve()}")
    (
        xr.concat([factors, unknown_factors.rename(new_factor="factor")], dim="factor")
        .rename("factors")
        .to_dataset()
        .to_zarr(args.output, mode="w")
    )
    return


def run_coverage_model(args: argparse.Namespace) -> None:
    """Determine parameters for coverage model given factors"""
    prefix = [bam_experiment_name(x) for x in args.tmpfiles]
    if len(prefix) > len(set(prefix)):
        raise ValueError("Passed tmpfiles have duplicate prefixes")
    log = get_logger()
    log.info("Setting up model matrix of all factors")
    factors = _get_factors(prefix, args)
    log.info(f"Opening coverage from {len(args.tmpfiles)} tmpfiles")
    coverage = _open_mf_with_prefix(args.tmpfiles)
    log.info("Solving for model parameters")
    coverage_model = mc.infer_model_params(
        coverage.psi, factors, extra_core_dims=["bootstrap_replicate"]
    )
    log.info(f"Saving model parameters to {args.coverage_model.resolve()}")
    (
        coverage_model.rename("coverage_model")
        .to_dataset()
        .to_zarr(args.coverage_model, mode="w")
    )
    return


def run_coverage_infer(args: argparse.Namespace) -> None:
    """Create corrected LSV coverage file"""
    prefix = bam_experiment_name(args.original_coverage)
    if prefix != bam_experiment_name(args.tmpfile):
        raise ValueError("Original coverage and tmpfile prefixes disagree")
    elif prefix != bam_experiment_name(args.corrected_coverage):
        raise ValueError("Original and corrected coverage prefixes disagree")
    log = get_logger()
    log.info(f"Setting up factors for {prefix = }")
    factors = _get_factors(prefix, args)
    log.info(f"Opening up model parameters from {args.coverage_model.resolve()}")
    coverage_model = xr.open_zarr(args.coverage_model).coverage_model
    df_ec = xr.open_zarr(args.original_coverage, group="events_coverage")
    df_e = xr.open_zarr(args.original_coverage, group="events")
    df_tmp = xr.open_zarr(args.tmpfile)
    corrected_bootstraps = mc.bootstraps_from_tmp_psi(
        mc.correct_with_model(df_tmp.psi, coverage_model, factors),
        df_tmp.total_coverage,
        df_tmp._offsets,
        df_ec.bootstraps,
    )
    df_e.to_zarr(args.corrected_coverage, mode="w", group="events")
    df_ec.assign(bootstraps=corrected_bootstraps).to_zarr(
        args.corrected_coverage, mode="a", group="events_coverage"
    )
    return


subcommand_tmpfile = GenericSubcommand(
    "(Advanced) Create temporary reparameterization of coverage for batch correction",
    args_tmpfile,
    run_tmpfile,
)
subcommand_factors_model = GenericSubcommand(
    "(Advanced) Create model of unknown confounding factors using tmpfiles and known factors",
    args_factors_model,
    run_factors_model,
)
subcommand_factors_infer = GenericSubcommand(
    "(Advanced) Save known/unknown factors for input experiment",
    args_factors_infer,
    run_factors_infer,
)
subcommand_coverage_model = GenericSubcommand(
    "(Advanced) Create model of coverage ratios using tmpfiles and factors",
    args_coverage_model,
    run_coverage_model,
)
subcommand_coverage_infer = GenericSubcommand(
    "(Advanced) Save corrected lsv coverage for input experiment",
    args_coverage_infer,
    run_coverage_infer,
)
