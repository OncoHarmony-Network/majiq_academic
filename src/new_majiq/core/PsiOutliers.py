"""
PsiOutliers.py

Identify aberrant splicing events of interests in cases vs controls

Author: Joseph K Aicher
"""

import re
from functools import cached_property
from typing import Dict, Final, List, Optional, Tuple

import numpy as np
import pandas as pd
import xarray as xr
from dask.distributed import progress

import new_majiq.constants as constants
from new_majiq.logger import get_logger

from .PsiControlsSummary import PsiControlsSummary
from .PsiCoverage import PsiCoverage
from .SpliceGraph import SpliceGraph


class PsiOutliers(object):
    """Outliers in PSI between cases and controls

    Parameters
    ----------
    cases: PsiCoverage
        PsiCoverage for prefixes that will be treated as cases
    controls: PsiControlsSummary
        Summary of PsiCoverage for prefixes that are treated as controls
    """

    def __init__(self, cases: PsiCoverage, controls: PsiControlsSummary):
        if not cases.events_df.equals(controls.events_df):
            raise ValueError(
                "cases and controls are not defined with the same events"
                f" ({cases.events_df=}, {controls.events_df=})"
            )
        if prefix_overlap := set(cases.prefixes) & set(controls.prefixes):
            # warning if there is overlap between prefixes
            get_logger().warning(
                "PsiOutliers cases and controls have overlapping prefixes"
                f" ({prefix_overlap})"
            )

        # save members of PsiOutliers
        self.cases: Final[PsiCoverage] = cases
        self.controls: Final[PsiControlsSummary] = controls
        return

    @property
    def cases_q(self) -> xr.DataArray:
        """Which quantiles will be evaluated for cases

        Notes
        -----
        We now use the same values of alpha when comparing cases and controls.
        Previously, these were kept as separate dimensions to evaluate using
        different quantiles in cases vs controls.
        This didn't appear to make a huge difference and optimizing this extra
        parameter would be more trouble than it's probably worth, but you could
        easily change the dimension name here (alpha to cases_alpha, etc.)
        in the future.
        """
        return self.controls.q

    @cached_property
    def cases_psi_quantile(self) -> xr.DataArray:
        """Posterior quantiles of PSI (approximate posterior)"""
        return self.cases.approximate_quantile(self.cases_q)

    @cached_property
    def dpsi_quantile_gap(self) -> xr.DataArray:
        """Gap in PSI between extreme quantiles of cases vs controls.

        Gap in PSI between extreme quantiles of cases vs controls.
        Positive values imply separation between the ranges (for same alpha),
        while negative values quantify amount of overlap.
        """
        gap_case_higher = self.cases_psi_quantile.sel(
            is_lb=True
        ) - self.controls.psi_quantile.sel(is_lb=False)
        gap_case_lower = self.controls.psi_quantile.sel(
            is_lb=True
        ) - self.cases_psi_quantile.sel(is_lb=False)
        return gap_case_higher.clip(min=gap_case_lower)

    @cached_property
    def tail_probability(self) -> xr.DataArray:
        """Probability that cases psi < ub quantile or psi > lb quantile

        Probability that cases psi is less than upper bound quantile or greater
        than lower bound quantile so as to make the probability low when cases
        are far from controls and high when cases are close to controls.

        Low confidence case PSI will be more difficult to have values close to
        0 or 1, while high confidence case PSI will more easily go to 0 or 1
        depending on how far the case is from the controls.
        """
        # evaluate CDF at lower and upper bound
        cdf_probabilities = self.cases.approximate_cdf(self.controls.psi_quantile)
        # probability for upper bound should be the cdf itself
        p_ub = cdf_probabilities.sel(is_lb=False)
        # probability for lower bound is 1 - cdf (P > lower bound)
        p_lb = 1 - cdf_probabilities.sel(is_lb=True)
        # when p_lb < 0.5, the control ub is less than case median
        # when p_ub < 0.5, the control lb is more than case median
        # otherwise, the median is between the control lb and ub.
        # Both values are >= 0.5, and we take the smaller value.
        # return p_lb.where(p_lb < 0.5, p_ub.where(p_ub < 0.5, p_ub.clip(min=p_lb)))
        return p_lb.clip(max=p_ub)  # take the minimum value of these

    @cached_property
    def dpsi_lb(self) -> xr.DataArray:
        """dPSI from case quantiles vs controls median.

        dPSI of case quantiles vs controls median.
        Negative values indicate that the controls median is within range of
        case quantiles.
        """
        gap_case_higher = (
            self.cases_psi_quantile.sel(is_lb=True) - self.controls.psi_median
        )
        gap_case_lower = self.controls.psi_median - self.cases_psi_quantile.sel(
            is_lb=False
        )
        return gap_case_higher.clip(min=gap_case_lower)

    @cached_property
    def dpsi_lb_scaled(self) -> xr.DataArray:
        """dPSI of case quantiles vs controls median scaled by controls interquantile range.

        dPSI of case quantiles vs controls median scaled by controls
        interquantile range.
        Negative values indicate that the controls median is within range of
        case quantiles.
        """
        return self.dpsi_lb / self.controls.psi_range.where(self.controls.psi_range > 0)

    @staticmethod
    def summarize_df_genes(df_events: pd.DataFrame) -> pd.DataFrame:
        """Summarize table from :meth:`PsiOutliers.summarize_df_events` to genes

        Parameters
        ----------
        df_events: pd.DataFrame
            Table from :meth:`PsiOutliers.summarize_df_events`.

        Returns
        -------
        pd.DataFrame
            Summary per gene of the event-level statistics.
        """
        # how to aggregate columns per event
        agg_kwargs: Dict[str, Tuple[str, str]] = {
            "gene_name": ("gene_name", "first"),
            "num_events": ("gene_name", "count"),
            "num_events_denovo_junction": ("has_denovo_junction", "sum"),
            "num_events_denovo_intron": ("has_denovo_intron", "sum"),
            "num_events_denovo_exon": ("has_denovo_exon", "sum"),
        }
        # how to summarize events_pass information, if present
        if "events_pass" in df_events.columns:
            df_events = df_events.assign(is_later_pass=lambda df: df.events_pass > 1)
            agg_kwargs = {
                "events_pass_max": ("events_pass", "max"),
                "num_later_pass": ("is_later_pass", "sum"),
                **agg_kwargs,
            }

        # identify cases to summarize
        pattern = re.compile(r"(.*)cases_raw_total")
        case_prefixes = [
            match.group(1)
            for match in (pattern.fullmatch(x) for x in df_events.columns)
            if match is not None
        ]
        # identify values of alpha to summarize
        pattern = re.compile(rf"{case_prefixes[0]}dpsi_lb(.*)")
        alpha_suffixes = [
            match.group(1)
            for match in (pattern.fullmatch(x) for x in df_events.columns)
            if match is not None
        ]

        # define aggregation for each case
        def get_case_kwargs(case_prefix: str) -> Dict[str, Tuple[str, str]]:
            case_agg_kwargs = {
                f"{case_prefix}cases_raw_total_max": (
                    f"{case_prefix}cases_raw_total_gene_max",
                    "first",
                ),
            }
            if len(case_prefixes) > 1:
                case_agg_kwargs = {
                    f"{case_prefix}num_passed": (
                        f"{case_prefix}cases_raw_total",
                        "count",
                    ),
                    **case_agg_kwargs,
                }
            for x in alpha_suffixes:
                case_agg_kwargs = {
                    **case_agg_kwargs,
                    f"{case_prefix}dpsi_quantile_gap{x}": (
                        f"{case_prefix}dpsi_quantile_gap{x}",
                        "max",
                    ),
                    f"{case_prefix}tail_probability{x}": (
                        f"{case_prefix}tail_probability{x}",
                        "min",
                    ),
                    f"{case_prefix}dpsi_lb{x}": (f"{case_prefix}dpsi_lb{x}", "max"),
                }
                # TODO add count above/below threshold as appropriate
            return case_agg_kwargs

        for case_prefix in case_prefixes:
            agg_kwargs = {**agg_kwargs, **get_case_kwargs(case_prefix)}

        # perform the aggregation
        df_genes = df_events.groupby("gene_id").agg(**agg_kwargs)

        # postprocessing for each case
        if len(case_prefixes) > 1:
            for case_prefix in case_prefixes:
                case_passed = df_genes[f"{case_prefix}num_passed"] > 0
                # mask float values to nan whenever no events passed for prefix
                for col in get_case_kwargs(case_prefix).keys():
                    if np.issubdtype(df_genes[col].dtype, np.floating):
                        df_genes[col] = df_genes[col].where(case_passed)

        # sorting for result
        if len(case_prefixes) == 1 and len(alpha_suffixes) == 1:
            df_genes.sort_values("dpsi_quantile_gap", ascending=False, inplace=True)
        if "events_pass_max" in df_genes.columns:
            df_genes.sort_values(
                # use stable sort to keep ordering from previous sort, if applicable
                "events_pass_max",
                ascending=False,
                kind="stable",
                inplace=True,
            )
        return df_genes

    @staticmethod
    def summarize_df_events(df_ecidx: pd.DataFrame) -> pd.DataFrame:
        """Summarize table in format of :meth:`PsiOutliers.to_dataframe` to events

        Parameters
        ----------
        df_ecidx: pd.DataFrame
            Table defined over event connections as created by
            :meth:`PsiOutliers.to_dataframe`.
            Requires columns as from :meth:`Events.ec_dataframe`.
            If present, also summarizes "events_pass" column as introduced by
            :meth:`Events.merge_dataframes`.

        Returns
        -------
        pd.DataFrame
            Table indexed by unique events defined by
            ["gene_id", "ref_exon_start", "ref_exon_end", "event_type"] with
            most extreme values of each statistic per event (note that most
            extreme values may not come from same connection).
        """
        # add additional columns to df_ecidx (copy)
        df_ecidx = df_ecidx.assign(
            is_denovo_junction=lambda df: df["is_denovo"] & ~df["is_intron"],
            is_denovo_intron=lambda df: df["is_denovo"] & df["is_intron"],
            denovo_exon=lambda df: df["ref_exon_denovo"] | df["other_exon_denovo"],
        )
        # how to aggregate columns per event
        agg_kwargs: Dict[str, Tuple[str, str]] = {
            # splicegraph
            "gene_name": ("gene_name", "first"),
            "event_size": ("gene_name", "count"),
            "has_intron": ("is_intron", "any"),
            "has_denovo_junction": ("is_denovo_junction", "any"),
            "has_denovo_intron": ("is_denovo_intron", "any"),
            "has_denovo_exon": ("denovo_exon", "any"),
            # controls, which always named same
            "num_passed": ("num_passed", "first"),
        }
        # if there are two passes, summarize this for the event, make it first
        if "events_pass" in df_ecidx.columns:
            agg_kwargs = {"events_pass": ("events_pass", "first"), **agg_kwargs}
        # identify cases to summarize
        pattern = re.compile(r"(.*)cases_raw_psi_mean")
        case_prefixes = [
            match.group(1)
            for match in (pattern.fullmatch(x) for x in df_ecidx.columns)
            if match is not None
        ]
        # identify values of alpha we have to work with
        pattern = re.compile(rf"{case_prefixes[0]}dpsi_lb(.*)")
        alpha_suffixes = [
            match.group(1)
            for match in (pattern.fullmatch(x) for x in df_ecidx.columns)
            if match is not None
        ]

        # define aggregation per case
        def get_case_kwargs(case_prefix: str) -> Dict[str, Tuple[str, str]]:
            case_agg_kwargs = {
                f"{case_prefix}passed": (f"{case_prefix}cases_raw_psi_mean", "count"),
                f"{case_prefix}cases_raw_total": (
                    f"{case_prefix}cases_raw_coverage",
                    "sum",
                ),
            }
            for x in alpha_suffixes:
                case_agg_kwargs = {
                    **case_agg_kwargs,
                    f"{case_prefix}dpsi_quantile_gap{x}": (
                        f"{case_prefix}dpsi_quantile_gap{x}",
                        "max",
                    ),
                    f"{case_prefix}tail_probability{x}": (
                        f"{case_prefix}tail_probability{x}",
                        "min",
                    ),
                    f"{case_prefix}dpsi_lb{x}": (f"{case_prefix}dpsi_lb{x}", "max"),
                }
            return case_agg_kwargs

        for case_prefix in case_prefixes:
            agg_kwargs = {**agg_kwargs, **get_case_kwargs(case_prefix)}

        # perform the aggregation
        df_events = df_ecidx.groupby(
            ["gene_id", "ref_exon_start", "ref_exon_end", "event_type"]
        ).agg(**agg_kwargs)

        # postprocessing for each case
        for case_prefix in case_prefixes:
            case_passed = df_events[f"{case_prefix}passed"] > 0
            df_events.drop(columns=f"{case_prefix}passed", inplace=True)
            for col in get_case_kwargs(case_prefix).keys():
                # mask columns if they did not pass for specified prefix
                if col in df_events.columns:
                    df_events[col] = df_events[col].where(case_passed)
            # summarize coverage for events along the gene
            df_events[f"{case_prefix}cases_raw_total_gene_max"] = df_events.groupby(
                "gene_id"
            )[f"{case_prefix}cases_raw_total"].transform("max")
            df_events[f"{case_prefix}cases_raw_total_pct"] = (
                df_events[f"{case_prefix}cases_raw_total"]
                / df_events[f"{case_prefix}cases_raw_total_gene_max"]
            )

        # sorting for result
        if len(case_prefixes) == 1 and len(alpha_suffixes) == 1:
            df_events.sort_values("dpsi_quantile_gap", ascending=False, inplace=True)
        if "events_pass" in df_events.columns:
            df_events.sort_values(
                # use stable sort to keep ordering from previous sort, if applicable
                "events_pass",
                ascending=False,
                kind="stable",
                inplace=True,
            )
        return df_events

    def to_dataframe(
        self,
        sg: Optional[SpliceGraph] = None,
        controls_min_experiments: float = constants.DEFAULT_OUTLIERS_MINEXPERIMENTS,
        alpha: Optional[float] = None,
        show_progress: bool = False,
    ) -> pd.DataFrame:
        """Return table of comparisons between cases and controls

        Parameters
        ----------
        sg: Optional[SpliceGraph]
            If provided, splicegraph with introns/junctions consistent with
            events used to annotate resulting dataframe
        controls_min_experiments: float
            Fraction (value < 1) or absolute number (value >= 1) of prefixes
            that must pass individually in control group for comparison to be
            made
        alpha: Optional[float]
            If provided, only show comparisons for specified value of alpha
        show_progress: bool
            show progress bar in dask if enabled

        Returns
        -------
        pd.DataFrame
        """
        # build tables with columns that will be concatenated together
        concat_df: List[pd.DataFrame] = list()
        # add dataframe with events annotations
        if sg is not None:
            concat_df.append(
                self.controls.get_events(sg.introns, sg.junctions).ec_dataframe()
            )
        ds = xr.Dataset(
            {
                # any of the cases passed
                "cases_passed": self.cases.event_passed.any("prefix"),
                "controls_passed": self.controls.passed_min_experiments(
                    controls_min_experiments
                ).isel(min_experiments_f=0, drop=True),
                "controls_psi_median": self.controls.psi_median,
                "cases_raw_psi_mean": self.cases.raw_psi_mean,
                "cases_raw_psi_std": self.cases.raw_psi_std,
                "cases_bootstrap_psi_std": self.cases.bootstrap_psi_std,
                "cases_raw_coverage": self.cases.raw_coverage,
                "controls_psi_quantile": self.controls.psi_quantile,
                "cases_psi_quantile": self.cases_psi_quantile,
                "dpsi_quantile_gap": self.dpsi_quantile_gap,
                "tail_probability": self.tail_probability,
                "dpsi_lb": self.dpsi_lb,
            }
        ).reset_coords("num_passed")
        if alpha is not None:
            try:
                ds = ds.sel(controls_alpha=[alpha])
            except KeyError:
                raise KeyError(
                    f"{alpha = } must be None"
                    f" or one of {ds.controls_alpha.values.tolist()}"
                )
        # load/compute into memory
        if show_progress:
            ds = ds.persist()
            progress(*(x.data for x in ds.variables.values() if x.chunks))
        ds = ds.load()
        # extract comparison information
        ds_comparison = ds[
            [
                k
                for k, v in ds.data_vars.items()
                if "controls_alpha" in v.dims and "is_lb" not in v.dims
            ]
        ]
        for i, prefix in enumerate(ds_comparison["prefix"].values):
            prefix_prefix = f"{prefix}-" if ds_comparison.sizes["prefix"] > 1 else ""
            for j, cur_alpha in enumerate(ds_comparison["controls_alpha"].values):
                alpha_suffix = (
                    f"-alpha_{cur_alpha:0.3f}"
                    if ds_comparison.sizes["controls_alpha"] > 1
                    else ""
                )
                concat_df.append(
                    ds_comparison.isel(prefix=i, controls_alpha=j, drop=True)
                    .to_dataframe()
                    .add_prefix(prefix_prefix)
                    .add_suffix(alpha_suffix)
                )
        # get ec_idx only variables (passed, controls psi_median)
        concat_df.append(
            ds[
                [k for k, v in ds.data_vars.items() if v.dims == ("ec_idx",)]
            ].to_dataframe()
        )
        # get ec_idx, prefix variables (psi_mean/psi_std, etc)
        ds_cases = ds[
            [k for k, v in ds.data_vars.items() if {"prefix", "ec_idx"} == set(v.dims)]
        ]
        for i, prefix in enumerate(ds_cases["prefix"].values):
            prefix_prefix = f"{prefix}-" if ds_comparison.sizes["prefix"] > 1 else ""
            concat_df.append(
                ds_cases.isel(prefix=i, drop=True)
                .to_dataframe()
                .add_prefix(prefix_prefix)
            )
        # extract quantile information for controls (i.e. no prefix)
        df_quantile = (
            ds[
                [
                    k
                    for k, v in ds.data_vars.items()
                    if "is_lb" in v.dims and "prefix" not in v.dims
                ]
            ]
            # index by quantiles, not alpha/is_lb
            .stack(_controls_q=ds["controls_q"].dims)
            .swap_dims(_controls_q="controls_q")
            .drop_vars("_controls_q")
            # to dataframe, put quantiles on columns
            .to_dataframe()
            .unstack("controls_q")
            .sort_index(axis=1)
        )
        df_quantile.columns = [f"{col}_{q:0.3f}" for col, q in df_quantile.columns]
        concat_df.append(df_quantile)
        # extract posterior quantile information for cases (i.e. has prefix)
        df_quantile = (
            ds[
                [
                    k
                    for k, v in ds.data_vars.items()
                    if "is_lb" in v.dims and "prefix" in v.dims
                ]
            ]
            # index by quantiles, not alpha/is_lb
            .stack(_controls_q=ds["controls_q"].dims)
            .swap_dims(_controls_q="controls_q")
            .drop_vars("_controls_q")
            # to dataframe, put quantiles on columns
            .to_dataframe()
            .unstack(["prefix", "controls_q"])
            .sort_index(axis=1)
        )
        if ds.sizes["prefix"] > 1:
            df_quantile.columns = [
                f"{prefix}-{col}_{q:0.3f}" for col, prefix, q in df_quantile.columns
            ]
        else:
            df_quantile.columns = [
                f"{col}_{q:0.3f}" for col, _, q in df_quantile.columns
            ]
        concat_df.append(df_quantile)
        del df_quantile
        return (
            pd.concat(concat_df, axis=1, join="inner")
            .loc[lambda df: df["cases_passed"] & df["controls_passed"]]
            .drop(columns=["cases_passed", "controls_passed"])
        )
