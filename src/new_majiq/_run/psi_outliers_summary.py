"""
psi_outliers_summary.py

Identify aberrant splicing events of interests in cases vs controls, and get
summary statistics over different thresholds, experiments

Author: Joseph K Aicher
"""

import argparse
from typing import List, Optional, cast

import numpy as np
import xarray as xr
from dask.delayed import Delayed
from dask.distributed import progress

import new_majiq as nm
from new_majiq._offsets import groupmax
from new_majiq._run._majiq_args import (
    ExistingResolvedPath,
    NewResolvedPath,
    StoreRequiredUniqueActionFactory,
    check_nonnegative_factory,
    check_range_factory,
    resources_args,
)
from new_majiq._run._run import GenericSubcommand
from new_majiq._stats import histogram
from new_majiq.core._workarounds import _load_zerodim_variables

DESCRIPTION = (
    "Compute outlier summary statistics from PSI cases (N ~ 1) vs PSI controls (N >> 1)"
)


def add_args(parser: argparse.ArgumentParser) -> None:
    # input information
    parser.add_argument(
        "cases",
        metavar="cases_psicov",
        type=ExistingResolvedPath,
        help="Path to PsiCoverage file with case experiments",
    )
    parser.add_argument(
        "controls",
        metavar="controls_summary",
        type=ExistingResolvedPath,
        help="Path to PsiControlsSummary file for comparison to controls",
    )
    parser.add_argument(
        "splicegraph",
        type=ExistingResolvedPath,
        help="Path to splicegraph over which control, case events were defined",
    )
    # output information
    parser.add_argument(
        "outliers",
        metavar="out_outliers",
        type=NewResolvedPath,
        help="Path for output zarr with summary of outliers given different"
        " thresholds",
    )
    # thresholds for case (posterior) quantiles
    parser.add_argument(
        "--alpha",
        metavar="A",
        # floats on [0, 1]
        type=check_range_factory(float, 0, 1, True, True),
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="Threshold for case quantiles in both directions (two-sided"
        " comparison). Quantiles are 0.5 * alpha and 1.0 - 0.5 * alpha."
        " If not specified, use same values as from input controls"
        " (default: %(default)s)",
    )
    parser.add_argument(
        "--min-experiments",
        metavar="X",
        type=check_nonnegative_factory(float, True),
        action=StoreRequiredUniqueActionFactory(),
        default=[nm.constants.DEFAULT_OUTLIERS_MINEXPERIMENTS],
        nargs="+",
        help="Only consider events for which the fraction (value < 1) or"
        " absolue number (value >= 1) of control experiments that must"
        " independently pass quantification filters in order for outlier status"
        " to be considered (default: %(default)s)",
    )
    parser.add_argument(
        "--min-gap-bins",
        metavar="N",
        type=check_nonnegative_factory(int, True),
        default=80,
        help="The number of bins on (0, 1] to identify the number of events or"
        " genes that were called as outliers",
    )
    parser.add_argument(
        "--events-mask",
        metavar="mask",
        type=ExistingResolvedPath,
        action=StoreRequiredUniqueActionFactory(),
        nargs="+",
        default=None,
        help="If specified, path to zarr xarray dataset from"
        " nm-shared-base-events with boolean mask over events that should be"
        " considered per case (masking for two-pass build) (default: None)",
    )

    # resources
    resources_args(parser, use_dask=True)
    return


def run(args: argparse.Namespace) -> None:
    log = nm.logger.get_logger()
    log.info(f"Loading cases from {args.cases}")
    cases = nm.PsiCoverage.from_zarr(args.cases)
    log.info(f"Loading controls summary from {args.controls}")
    controls = nm.PsiControlsSummary.from_zarr(args.controls)
    log.info(f"Comparing cases {cases} vs {controls}")
    outliers = nm.PsiOutliers(cases, controls)
    outliers_ds = outliers.dataset(cases_alpha=args.alpha)
    dpsi_gap = outliers_ds["dpsi_quantile_gap"]
    log.info(f"Applying controls min-experiments thresholds ({args.min_experiments})")
    passed = controls.passed_min_experiments(args.min_experiments)
    dpsi_gap = dpsi_gap.where(passed)
    if args.events_mask:
        log.info(f"Applying mask from {args.events_mask}")
        mask = xr.open_mfdataset(
            args.events_mask,
            engine="zarr",
            combine="nested",
            concat_dim="derived_name",
            join="override",
            compat="override",
            coords="minimal",
            data_vars="minimal",
        )["ec_mask"].sel(derived_name=outliers_ds["prefix"])
        dpsi_gap = dpsi_gap.where(mask)
    log.info(f"Loading events from {args.cases} and {args.splicegraph}")
    sg = nm.SpliceGraph.from_zarr(args.splicegraph)
    events = cases.get_events(sg.introns, sg.junctions)
    log.info("Summarizing gaps between quantiles to event, gene levels")
    dpsi_gap = dpsi_gap.chunk({"ec_idx": None})
    dpsi_gap_events = xr.apply_ufunc(
        groupmax,
        dpsi_gap,
        cases.lsv_idx.load(),
        xr.DataArray(np.empty(events.num_events), dims=["e_idx"]),
        input_core_dims=[["ec_idx"], ["ec_idx"], ["e_idx"]],
        output_core_dims=[["e_idx"]],
        dask="parallelized",
        output_dtypes=[dpsi_gap.dtype],
    ).chunk({"e_idx": None})
    dpsi_gap_genes = xr.apply_ufunc(
        groupmax,
        dpsi_gap,
        xr.DataArray(events.connection_gene_idx(), dims=["ec_idx"]),
        xr.DataArray(np.empty(len(sg.genes)), dims=["gene_idx"]),
        input_core_dims=[["ec_idx"], ["ec_idx"], ["gene_idx"]],
        output_core_dims=[["gene_idx"]],
        dask="parallelized",
        output_dtypes=[dpsi_gap.dtype],
    ).chunk({"gene_idx": None})
    log.info("Summarizing gaps between quantiles to event, gene levels")
    min_gap_bins = xr.DataArray(
        _ := np.linspace(0, 1, 1 + args.min_gap_bins, dtype=dpsi_gap.dtype)[:-1],
        [("min_gap", _)],
    )
    events_quantified = dpsi_gap_events.count("e_idx")
    genes_quantified = dpsi_gap_genes.count("gene_idx")
    dpsi_gap_events_bins = xr.apply_ufunc(
        histogram,
        dpsi_gap_events,
        # count events on (0, 1] (not [0, 1))
        np.array(np.finfo(dpsi_gap_events.dtype).tiny, dtype=dpsi_gap_events.dtype),
        np.array(1 + np.finfo(dpsi_gap_events.dtype).tiny, dtype=dpsi_gap_events.dtype),
        min_gap_bins,
        input_core_dims=[["e_idx"], [], [], ["min_gap"]],
        output_core_dims=[["min_gap"]],
        dask="parallelized",
        output_dtypes=[np.int64],
    )
    dpsi_gap_genes_bins = xr.apply_ufunc(
        histogram,
        dpsi_gap_genes,
        # count genes on (0, 1] (not [0, 1))
        np.array(np.finfo(dpsi_gap_genes.dtype).tiny, dtype=dpsi_gap_genes.dtype),
        np.array(1 + np.finfo(dpsi_gap_genes.dtype).tiny, dtype=dpsi_gap_genes.dtype),
        min_gap_bins,
        input_core_dims=[["gene_idx"], [], [], ["min_gap"]],
        output_core_dims=[["min_gap"]],
        dask="parallelized",
        output_dtypes=[np.int64],
    )
    log.info("Calculating number of outlier for different thresholds on quantile gap")

    def reverse_cumsum_gufunc(x: np.ndarray) -> np.ndarray:
        """cumsum in reverse direction on last axis"""
        out = np.empty(x.shape, dtype=x.dtype)
        np.cumsum(x[..., ::-1], out=out[..., ::-1], axis=-1)
        return out

    dpsi_gap_events_sf = xr.apply_ufunc(
        reverse_cumsum_gufunc,
        dpsi_gap_events_bins,
        input_core_dims=[["min_gap"]],
        output_core_dims=[["min_gap"]],
        dask="parallelized",
        output_dtypes=[dpsi_gap_events_bins.dtype],
    )
    dpsi_gap_genes_sf = xr.apply_ufunc(
        reverse_cumsum_gufunc,
        dpsi_gap_genes_bins,
        input_core_dims=[["min_gap"]],
        output_core_dims=[["min_gap"]],
        dask="parallelized",
        output_dtypes=[dpsi_gap_genes_bins.dtype],
    )

    log.info(f"Saving comparisons to {args.outliers}")
    save_ds = xr.Dataset(
        {
            "quantified_events": events_quantified,
            "quantified_genes": genes_quantified,
            "outlier_events": dpsi_gap_events_sf,
            "outlier_genes": dpsi_gap_genes_sf,
        }
    )
    save_ds_future = cast(
        Delayed,
        save_ds.pipe(_load_zerodim_variables).to_zarr(
            args.outliers,
            mode="w",
            consolidated=True,
            compute=False,
        ),
    )
    if args.show_progress:
        save_ds_future = save_ds_future.persist()
        progress(save_ds_future)
    else:
        save_ds_future.compute()
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
