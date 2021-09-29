"""
moccasin.py

Perform batch correction using the moccasin algorithm

Author: Joseph K Aicher
"""

import argparse
from pathlib import Path
from typing import List, Union

import moccasin.moccasin as mc
import numpy as np
import pandas as pd
import xarray as xr
from dask.distributed import Client

import new_majiq as nm
import new_majiq.constants as constants
from new_majiq._offsets import clip_and_normalize_strict
from new_majiq._run._majiq_args import check_nonnegative_factory
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger


def _args_dask(parser: argparse.ArgumentParser) -> None:
    """arguments to pass to Dask scheduler"""
    parser.add_argument(
        "--nthreads",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_QUANTIFY_NTHREADS,
        help="Number of threads used by Dask scheduler to fit models in chunks"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--memory-limit",
        type=str,
        default="auto",
        help="Memory limit to pass to dask cluster (default: %(default)s)",
    )
    # doesn't appear we need to manually limit memory at this time, can do so
    # in the future if necessary
    return


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
        help="Paths to factors matrices including desired prefixes. If a prefix"
        " is in multiple files, the first listed factors file with that prefix"
        " is used",
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
        if missing := set(x for x in args.confounding if x not in df.columns):
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
        factors = xr.open_mfdataset(
            args.factors_zarr,
            engine="zarr",
            combine="nested",
            concat_dim="prefix",
            join="override",
            compat="override",  # use first file if overlapping prefixes
            coords="minimal",
            data_vars="minimal",
        )["factors"]
        # load factors for desired prefixes into memory (as other methods do)
        factors = factors.sel(prefix=prefix).load()
    # do the squeeze, if desired
    if do_squeeze:
        factors = factors.squeeze("prefix")
    log = get_logger()
    log.info(
        f"Using {factors.sizes['factor']} factors "
        f" ({factors['confounding'].sum().values[()]} confounding):"
        + "".join(
            [
                f"\nFactor {name}\t"
                + ("confounding" if is_confounding else "noncounfounding")
                for name, is_confounding in zip(
                    factors["factor"].values, factors["confounding"].values
                )
            ]
        )
    )
    return factors


def _args_ruv_parameters(parser: argparse.ArgumentParser) -> None:
    ruv = parser.add_argument_group("Parameters for modeling unknown confounders")
    ruv.add_argument(
        "--ruv-max-new-factors",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_MOCCASIN_RUV_MAX_FACTORS,
        help="Maximum number of new factors the model will add (default: %(default)s)",
    )
    ruv.add_argument(
        "--ruv-max-events",
        type=check_nonnegative_factory(int, True),
        default=constants.DEFAULT_MOCCASIN_RUV_MAX_EVENTS,
        help="Maximum number of events to model residuals for (default: %(default)s)",
    )
    return


def args_factors_model(parser: argparse.ArgumentParser) -> None:
    """arguments for creating model of unknown factors"""
    _args_factors(parser)  # get factors information
    parser.add_argument(
        "factors_model", type=Path, help="Path for output model for unknown confounders"
    )
    parser.add_argument(
        "psicov",
        type=Path,
        nargs="+",
        help="Paths for input psi coverage files for unknown confounders model",
    )
    _args_ruv_parameters(parser)
    _args_dask(parser)
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
        "output",
        type=Path,
        help="Path for output factors (known and unknown)",
    )
    parser.add_argument(
        "psicov",
        type=Path,
        nargs="+",
        help="Path to psi coverage files for which known/unknown factors will be saved",
    )
    _args_dask(parser)
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
        "psicov",
        nargs="+",
        type=Path,
        help="Paths for input psi coverage files for coverage model",
    )
    _args_dask(parser)
    return


def args_coverage_infer(parser: argparse.ArgumentParser) -> None:
    """arguments for getting corrected lsv coverage"""
    _args_factors(parser)
    parser.add_argument(
        "coverage_model",
        type=Path,
        help="Input path for coverage model",
    )
    parser.add_argument(
        "corrected_psicov",
        type=Path,
        help="Path for output corrected psi coverage",
    )
    parser.add_argument(
        "original_psicov",
        type=Path,
        nargs="+",
        help="Paths to original, uncorrected psi coverage",
    )
    _args_dask(parser)
    return


def run_factors_model(args: argparse.Namespace) -> None:
    """model unknown factors using input coverage"""
    log = get_logger()
    client = Client(
        n_workers=1,
        threads_per_worker=args.nthreads,
        dashboard_address=None,
        memory_limit=args.memory_limit,
    )
    log.info(client)
    log.info(f"Opening coverage from {len(args.psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    log.info("Setting up model matrix of known factors")
    factors = _get_factors(psicov.prefixes, args)
    log.info("Learning model for unknown confounding factors")
    model = mc.ModelUnknownConfounders.train(
        psicov.bootstrap_psi,
        psicov.event_passed,
        psicov.lsv_offsets.values,
        factors,
        max_new_factors=args.ruv_max_new_factors,
        max_events=args.ruv_max_events,
        dim_prefix="prefix",
        dim_factor="factor",
        dim_ecidx="ec_idx",
        log=log,
    )
    # report about explained variance
    explained_variance = model.explained_variance.values
    if len(explained_variance):
        log.info(
            f"Learned {model.num_factors} factors explaining"
            f" {explained_variance.sum():.3%} of residual variance"
            " (after accounting for known factors):"
            + "".join(
                [
                    f"\nfactor {i}\t{x:.3%} (cumulative: {cum_x:.3%})"
                    for i, (x, cum_x) in enumerate(
                        zip(explained_variance, explained_variance.cumsum()), 1
                    )
                ]
            )
        )
    # save model
    log.info(f"Saving model to {args.factors_model.resolve()}")
    model.to_zarr(args.factors_model)
    return


def run_factors_infer(args: argparse.Namespace):
    """compute unknown factors using input coverage"""
    log = get_logger()
    client = Client(
        n_workers=1,
        threads_per_worker=args.nthreads,
        dashboard_address=None,
        memory_limit=args.memory_limit,
    )
    log.info(client)
    log.info(f"Opening coverage from {len(args.psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    log.info("Setting up model matrix of known factors")
    factors = _get_factors(psicov.prefixes, args).load()
    log.info(f"Loading model for unknown factors from {args.factors_model.resolve()}")
    model = mc.ModelUnknownConfounders.from_zarr(args.factors_model)
    log.info("Solving for unknown confounders")
    unknown_factors = model.predict(
        psicov.bootstrap_psi, psicov.event_passed, factors
    ).load()
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
    log = get_logger()
    client = Client(
        n_workers=1,
        threads_per_worker=args.nthreads,
        dashboard_address=None,
        memory_limit=args.memory_limit,
    )
    log.info(client)
    log.info(f"Opening coverage from {len(args.psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    log.info("Setting up model matrix of all factors")
    factors = _get_factors(psicov.prefixes, args).load()
    log.info("Solving for bootstrap model parameters")
    bootstrap_model = mc.infer_model_params(
        psicov.bootstrap_psi,
        psicov.event_passed,
        factors,
        complete=False,
        dim_prefix="prefix",
        dim_factor="factor",
    )
    log.info("Solving for raw model parameters")
    # TODO could potentially share work if combined bootstrap_psi, raw_psi to
    # solve for parameters together
    raw_model = mc.infer_model_params(
        psicov.raw_psi,
        psicov.event_passed,
        factors,
        complete=False,
        dim_prefix="prefix",
        dim_factor="factor",
    )
    log.info(f"Saving model parameters to {args.coverage_model.resolve()}")
    xr.Dataset(
        {
            "bootstrap_model": bootstrap_model,
            "raw_model": raw_model,
        }
    ).to_zarr(args.coverage_model, mode="w")
    return


def run_coverage_infer(args: argparse.Namespace) -> None:
    """Create corrected LSV coverage file"""
    log = get_logger()
    client = Client(
        n_workers=1,
        threads_per_worker=args.nthreads,
        dashboard_address=None,
        memory_limit=args.memory_limit,
    )
    log.info(client)
    log.info(f"Opening coverage from {len(args.original_psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.original_psicov)
    log.info("Setting up model matrix of all factors")
    factors = _get_factors(psicov.prefixes, args).load().chunk({"prefix": 1})
    log.info(f"Opening up model parameters from {args.coverage_model.resolve()}")
    models = xr.open_zarr(args.coverage_model).load()
    log.info("Loading offsets between events")
    offsets = psicov.lsv_offsets.load()

    def clip_and_renormalize_psi(x: xr.DataArray) -> xr.DataArray:
        return xr.apply_ufunc(
            clip_and_normalize_strict,
            x.chunk({"ec_idx": -1}),
            offsets,
            input_core_dims=[["ec_idx"], ["offset_idx"]],
            output_core_dims=[["ec_idx"]],
            output_dtypes=(x.dtype,),
            dask="parallelized",
        )

    log.info("Correcting bootstrap_psi")
    adj_bootstrap_psi = (
        # get linear adjustment of bootstrap_psi
        mc.correct_with_model(
            psicov.bootstrap_psi,
            models.bootstrap_model,
            factors,
            dim_prefix="prefix",
            dim_factor="factor",
        )
        # clip and renormalize
        .pipe(clip_and_renormalize_psi)
        # but must be passed and not null
        .pipe(
            lambda x: x.where(x.notnull() & psicov.event_passed, psicov.bootstrap_psi)
        )
    )
    log.info("Correcting raw_psi")
    adj_raw_psi = (
        # get linear adjustment of raw_psi
        mc.correct_with_model(
            psicov.raw_psi,
            models.raw_model,
            factors,
            dim_prefix="prefix",
            dim_factor="factor",
        )
        # clip and renormalize
        .pipe(clip_and_renormalize_psi)
        # but must be passed and not null
        .pipe(lambda x: x.where(x.notnull() & psicov.event_passed, psicov.raw_psi))
    )
    log.info("Updating PsiCoverage")
    psicov = psicov.updated(
        adj_bootstrap_psi,
        adj_raw_psi,
        model_path=f"{args.coverage_model.resolve()}",
    )
    log.info(f"Saving corrected coverage to {args.corrected_psicov.resolve()}")
    psicov.to_zarr(args.corrected_psicov)
    return


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
    "(Advanced) Save corrected psi coverage for input experiment",
    args_coverage_infer,
    run_coverage_infer,
)
