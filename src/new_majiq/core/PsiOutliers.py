"""
PsiOutliers.py

Identify aberrant splicing events of interests in cases vs controls

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import Final, List, Optional

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
        if not cases.events.equals(controls.events):
            raise ValueError(
                "cases and controls are not defined with the same events"
                f" ({cases.events=}, {controls.events=})"
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
    def cases_tail_probability(self) -> xr.DataArray:
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
    def cases_dpsi_lb(self) -> xr.DataArray:
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
    def cases_dpsi_lb_scaled(self) -> xr.DataArray:
        """dPSI of case quantiles vs controls median scaled by interquantile range.

        dPSI of case quantiles vs controls median scaled by interquantile range.
        Negative values indicate that the controls median is within range of
        case quantiles.
        """
        return self.cases_dpsi_lb / self.controls.psi_range.where(
            self.controls.psi_range > 0
        )

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
                self.controls.get_events(sg.introns, sg.junctions).ec_dataframe
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
                "cases_tail_probability": self.cases_tail_probability,
                "cases_dpsi_lb": self.cases_dpsi_lb,
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
