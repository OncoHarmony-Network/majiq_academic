"""
moccasin.py

Perform batch correction using the moccasin algorithm

Author: Joseph K Aicher
"""

import argparse
import sys
from typing import Dict, List, Union

import moccasin.moccasin as mc
import numpy as np
import pandas as pd
import xarray as xr
from dask.distributed import progress

import new_majiq as nm
from new_majiq._offsets import clip_and_normalize_strict
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_nonnegative_factory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq.logger import get_logger


def _args_factors(parser: argparse.ArgumentParser) -> None:
    """add arguments to specify factors to load"""
    factors = parser.add_argument_group("required input confounders arguments")
    factors_ex = factors.add_mutually_exclusive_group(required=True)
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
        metavar="TSV",
        dest="factors_tsv",
        type=ExistingResolvedPath,
        default=None,
        help="Path to TSV with model matrix for all prefixes being processed"
        " (required column: prefix)",
    )
    factors_ex.add_argument(
        "--factors-zarr",
        metavar="ZARR",
        dest="factors_zarr",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="Paths to factors matrices including desired prefixes. If a prefix"
        " is in multiple files, the first listed factors file with that prefix"
        " is used",
    )
    parser.add_argument(
        "--confounding",
        metavar="VAR",
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


def _args_ruv_parameters(
    parser: argparse.ArgumentParser,
    default_max_new_factors: int = nm.constants.DEFAULT_MOCCASIN_RUV_MAX_FACTORS,
) -> None:
    ruv = parser.add_argument_group("unknown confounders modeling arguments")
    ruv.add_argument(
        "--ruv-max-new-factors",
        metavar="N",
        type=check_nonnegative_factory(int, True),
        default=default_max_new_factors,
        help="Maximum number of new factors the model will add (default: %(default)s)",
    )
    ruv.add_argument(
        "--ruv-max-events",
        metavar="N",
        type=check_nonnegative_factory(int, True),
        default=nm.constants.DEFAULT_MOCCASIN_RUV_MAX_EVENTS,
        help="Maximum number of events to model residuals for (default: %(default)s)",
    )
    return


def args_factors_model(parser: argparse.ArgumentParser) -> None:
    """arguments for creating model of unknown factors"""
    _args_factors(parser)  # get factors information
    parser.add_argument(
        "factors_model",
        type=NewResolvedPath,
        help="Path for output model for unknown confounders",
    )
    parser.add_argument(
        "psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths for input psi coverage files for unknown confounders model",
    )
    _args_ruv_parameters(parser)
    resources_args(parser, use_dask=True)
    return


def args_factors_infer(parser: argparse.ArgumentParser) -> None:
    """arguments for inferring factors using factors model"""
    _args_factors(parser)
    parser.add_argument(
        "factors_model",
        type=ExistingResolvedPath,
        help="Path for input model of unknown confounders",
    )
    parser.add_argument(
        "output",
        type=NewResolvedPath,
        help="Path for output factors (known and unknown)",
    )
    parser.add_argument(
        "psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Path to psi coverage files for which known/unknown factors will be saved",
    )
    resources_args(parser, use_dask=True)
    return


def args_coverage_model(parser: argparse.ArgumentParser) -> None:
    """arguments for modeling coverage"""
    parser.add_argument(
        "coverage_model",
        type=NewResolvedPath,
        help="Output path for coverage model",
    )
    _args_factors(parser)
    parser.add_argument(
        "psicov",
        nargs="+",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        help="Paths for input psi coverage files for coverage model",
    )
    resources_args(parser, use_dask=True)
    return


def args_coverage_infer(parser: argparse.ArgumentParser) -> None:
    """arguments for getting corrected lsv coverage"""
    _args_factors(parser)
    parser.add_argument(
        "coverage_model",
        type=ExistingResolvedPath,
        help="Input path for coverage model",
    )
    parser.add_argument(
        "corrected_psicov",
        type=NewResolvedPath,
        help="Path for output corrected psi coverage",
    )
    parser.add_argument(
        "original_psicov",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        help="Paths to original, uncorrected psi coverage",
    )
    resources_args(parser, use_dask=True)
    return


def run_factors_model(args: argparse.Namespace) -> None:
    """model unknown factors using input coverage"""
    log = get_logger()
    log.info(f"Opening coverage from {len(args.psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    log.info("Setting up model matrix of known factors")
    factors = _get_factors(psicov.prefixes, args)
    log.info(f"Learning model for unknown confounding factors from {psicov}")
    model = mc.ModelUnknownConfounders.train(
        psicov.raw_psi,
        psicov.event_passed,
        psicov.lsv_offsets.values,
        factors,
        max_new_factors=args.ruv_max_new_factors,
        max_events=args.ruv_max_events,
        dim_prefix="prefix",
        dim_factor="factor",
        dim_ecidx="ec_idx",
        log=log,
        show_progress=args.show_progress,
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
    log.info(f"Saving model to {args.factors_model}")
    model.to_zarr(args.factors_model)
    return


def run_factors_infer(args: argparse.Namespace):
    """compute unknown factors using input coverage"""
    log = get_logger()
    log.info(f"Opening coverage from {len(args.psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    log.info("Setting up model matrix of known factors")
    factors = _get_factors(psicov.prefixes, args).load()
    log.info(f"Loading model for unknown factors from {args.factors_model}")
    model = mc.ModelUnknownConfounders.from_zarr(args.factors_model)
    log.info(f"Solving for unknown confounders for {psicov}")
    unknown_factors = model.predict(psicov.raw_psi, psicov.event_passed, factors)
    if args.show_progress:
        unknown_factors = unknown_factors.persist()
        progress(unknown_factors.data)
    unknown_factors = unknown_factors.load()

    log.info(f"Saving combined factors to {args.output}")
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
    log.info(f"Opening coverage from {len(args.psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    log.info("Setting up model matrix of all factors")
    factors = _get_factors(psicov.prefixes, args).load()
    log.info(f"Solving for model parameters from {psicov}")
    bootstrap_model = mc.infer_model_params(
        psicov.bootstrap_psi,
        psicov.event_passed,
        factors,
        complete=False,
        dim_prefix="prefix",
        dim_factor="factor",
    )
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
    models = xr.Dataset(dict(bootstrap_model=bootstrap_model, raw_model=raw_model))
    if args.show_progress:
        models = models.persist()
        progress(*(x.data for x in models.variables.values() if x.chunks))
    log.info(f"Saving model parameters to {args.coverage_model}")
    models.to_zarr(args.coverage_model, mode="w")
    return


def run_coverage_infer(args: argparse.Namespace) -> None:
    """Create corrected LSV coverage file"""
    log = get_logger()
    log.info(f"Opening coverage from {len(args.original_psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.original_psicov)
    log.info("Setting up model matrix of all factors")
    factors = _get_factors(psicov.prefixes, args).load().chunk({"prefix": 1})
    log.info(f"Opening up model parameters from {args.coverage_model}")
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
    log.info(f"Updating {psicov}")
    psicov = psicov.updated(
        adj_bootstrap_psi,
        adj_raw_psi,
        model_path=f"{args.coverage_model}",
    )
    log.info(f"Saving corrected coverage to {args.corrected_psicov}")
    psicov.to_zarr(args.corrected_psicov, show_progress=args.show_progress)
    return


def args_pipeline(parser: argparse.ArgumentParser) -> None:
    """arguments for pipeline through model/inference steps of factors/coverage"""
    _args_factors(parser)
    _args_ruv_parameters(parser, default_max_new_factors=0)
    parser.add_argument(
        "output_dir",
        type=NewResolvedPath,
        help="Path for new directory for output files",
    )
    parser.add_argument(
        "psicov",
        nargs="+",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        help="Paths for input psi coverage files to perform batch correction on",
    )
    factor_prefixes = parser.add_argument_group("factor modeling experiments arguments")
    factor_prefixes_ex = factor_prefixes.add_mutually_exclusive_group()
    factor_prefixes_ex.add_argument(
        "--factors-prefixes-include",
        nargs="+",
        type=str,
        default=None,
        action=StoreRequiredUniqueActionFactory(),
        help="Prefixes to explicitly include in factors-model step"
        " (default: include all)",
    )
    factor_prefixes_ex.add_argument(
        "--factors-prefixes-exclude",
        nargs="+",
        type=str,
        default=None,
        action=StoreRequiredUniqueActionFactory(),
        help="Prefixes to explicitly exclude in factors-model step,"
        " i.e. use all others (default: include all)",
    )
    coverage_prefixes = parser.add_argument_group(
        "coverage modeling experiments arguments"
    )
    coverage_prefixes_ex = coverage_prefixes.add_mutually_exclusive_group()
    coverage_prefixes_ex.add_argument(
        "--coverage-prefixes-include",
        nargs="+",
        type=str,
        default=None,
        action=StoreRequiredUniqueActionFactory(),
        help="Prefixes to explicitly include in coverage-model step"
        " (default: include all)",
    )
    coverage_prefixes_ex.add_argument(
        "--coverage-prefixes-exclude",
        nargs="+",
        type=str,
        default=None,
        action=StoreRequiredUniqueActionFactory(),
        help="Prefixes to explicitly exclude in coverage-model step,"
        " i.e. use all others (default: include all)",
    )
    resources_args(parser, use_dask=True)
    return


def run_pipeline(args: argparse.Namespace) -> None:
    """model factors, get updated factors, model coverage, get updated coverage"""
    log = get_logger()
    log.info(f"Opening coverage from {len(args.psicov)} PSI coverage files")
    psicov = nm.PsiCoverage.from_zarr(args.psicov)
    psicov_subset: nm.PsiCoverage  # used for modeling steps
    log.info(f"Correcting coverage for {psicov}")
    log.info("Setting up model matrix of known factors")
    factors = _get_factors(psicov.prefixes, args)

    if args.ruv_max_new_factors > 0:
        log.info(f"Modeling up to {args.ruv_max_new_factors} unknown factors")
        if args.factors_prefixes_include:
            log.info("Modeling factors using specified prefixes")
            psicov_subset = psicov[args.factors_prefixes_include]
        elif args.factors_prefixes_exclude:
            log.info("Modeling factors excluding specified prefixes")
            psicov_subset = psicov[
                [x for x in psicov.prefixes if x not in args.factors_prefixes_exclude]
            ]
        else:
            psicov_subset = psicov
        log.info(f"Learning model for unknown confounding factors from {psicov_subset}")
        factors_model = mc.ModelUnknownConfounders.train(
            psicov_subset.raw_psi,
            psicov_subset.event_passed,
            psicov_subset.lsv_offsets.values,
            factors.sel(prefix=psicov_subset.prefixes),
            max_new_factors=args.ruv_max_new_factors,
            max_events=args.ruv_max_events,
            dim_prefix="prefix",
            dim_factor="factor",
            dim_ecidx="ec_idx",
            log=log,
            show_progress=args.show_progress,
        )
        # report about explained variance
        explained_variance = factors_model.explained_variance.values
        if len(explained_variance):
            log.info(
                f"Learned {factors_model.num_factors} factors explaining"
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
        factors_model_path = args.output_dir / "factors_model.zarr"
        log.info(f"Saving factors model to {factors_model_path}")
        factors_model.to_zarr(factors_model_path)

        # infer factors for all experiments
        log.info(f"Inferring unknown factors for {psicov}")
        unknown_factors = factors_model.predict(
            psicov.raw_psi, psicov.event_passed, factors
        )
        if args.show_progress:
            unknown_factors = unknown_factors.persist()
            progress(unknown_factors.data)
        unknown_factors = unknown_factors.load()
        # update factors
        factors = xr.concat(
            [factors, unknown_factors.rename(new_factor="factor")], dim="factor"
        )
        factors_path = args.output_dir / "factors.zarr"
        log.info(f"Saving combined known/unknown factors to {factors_path}")
        factors.rename("factors").to_dataset().to_zarr(factors_path, mode="w")
    # done modeling/updating factors if modeling unknown factors

    if args.coverage_prefixes_include:
        log.info("Modeling coverage using specified prefixes")
        psicov_subset = psicov[args.coverage_prefixes_include]
    elif args.coverage_prefixes_exclude:
        log.info("Modeling coverage excluding specified prefixes")
        psicov_subset = psicov[
            [x for x in psicov.prefixes if x not in args.coverage_prefixes_exclude]
        ]
    else:
        psicov_subset = psicov
    log.info(f"Learning model for observed coverage from {psicov_subset}")
    bootstrap_model = mc.infer_model_params(
        psicov_subset.bootstrap_psi,
        psicov_subset.event_passed,
        factors.sel(prefix=psicov_subset.prefixes),
        complete=False,
        dim_prefix="prefix",
        dim_factor="factor",
    )
    raw_model = mc.infer_model_params(
        psicov_subset.raw_psi,
        psicov_subset.event_passed,
        factors.sel(prefix=psicov_subset.prefixes),
        complete=False,
        dim_prefix="prefix",
        dim_factor="factor",
    )
    models = xr.Dataset(dict(bootstrap_model=bootstrap_model, raw_model=raw_model))
    if args.show_progress:
        models = models.persist()
        progress(*(x.data for x in models.variables.values() if x.chunks))
    models.load().chunk({})
    models_path = args.output_dir / "coverage_model.zarr"
    log.info(f"Saving coverage model to {models_path}")
    models.to_zarr(models_path, mode="w")

    log.info(f"Updating coverage for {psicov}")
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
    # replace original coverage with adjusted coverage
    psicov = psicov.updated(adj_bootstrap_psi, adj_raw_psi, model_path=f"{models_path}")
    psicov_path = args.output_dir / "corrected.psicov"
    log.info(f"Saving corrected coverage to {psicov_path}")
    psicov.to_zarr(psicov_path, show_progress=args.show_progress)
    return


subcommand_factors_model = GenericSubcommand(
    "Build model of unknown confounding factors using input PsiCoverage",
    args_factors_model,
    run_factors_model,
)
subcommand_factors_infer = GenericSubcommand(
    "Use unknown confounding factors model on inputs to infer known/unknown factors",
    args_factors_infer,
    run_factors_infer,
)
subcommand_coverage_model = GenericSubcommand(
    "Build model of PsiCoverage using input coverage and factors",
    args_coverage_model,
    run_coverage_model,
)
subcommand_coverage_infer = GenericSubcommand(
    "Use PsiCoverage model to infer corrected coverage given inputs",
    args_coverage_infer,
    run_coverage_infer,
)
subcommand_pipeline = GenericSubcommand(
    "majiq-moccasin pipeline for PsiCoverage batch correction with"
    " known/unknown factors",
    args_pipeline,
    run_pipeline,
)


SUBPARSER_SOURCES: Dict[str, GenericSubcommand] = {
    "factors-model": subcommand_factors_model,
    "factors-infer": subcommand_factors_infer,
    "coverage-model": subcommand_coverage_model,
    "coverage-infer": subcommand_coverage_infer,
    "pipeline": subcommand_pipeline,
}


def main() -> None:
    """Entry-point into multiple tools using subcommands"""
    # build parser
    parser = argparse.ArgumentParser(
        description="Tools to model factors and PsiCoverage using MOCCASIN"
    )
    # add subparsers
    subparsers = parser.add_subparsers(required=True, help="")
    for src_name, src_module in SUBPARSER_SOURCES.items():
        src_parser = subparsers.add_parser(
            src_name,
            help=src_module.DESCRIPTION,
            description=src_module.DESCRIPTION,
        )
        src_parser.set_defaults(func=src_module.run)
        src_module.add_args(src_parser)

    # check length of input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # parse arguments now
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
