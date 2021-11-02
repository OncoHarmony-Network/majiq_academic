"""
PsiOutliers.py

Identify aberrant splicing events of interests in cases vs controls

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import Final, Optional, Sequence, Union

import xarray as xr

import new_majiq.constants as constants
from new_majiq.logger import get_logger
from new_majiq.PsiCoverage import PsiCoverage


class PsiOutliers(object):
    """Outliers in PSI between cases and controls

    Parameters
    ----------
    cases: PsiCoverage
        PsiCoverage for prefixes that will be treated as cases
    controls: PsiCoverage
        PsiCoverage for prefixes that will be treated as controls
    controls_min_experiments_f: float
        min_experiments used for masking if there are enough controls to
        perform comparison against with cases. If less than 1, a fraction of
        the total number of prefixes, otherwise the absolute number that must
        be passed. Generally set to high percentage so that quantiles are
        reasonable.
    alpha: Union[float, Sequence[float]]
        Threshold for case and control quantiles in both directions (two-sided
        comparison). Quantiles are 0.5 * alpha and 1 - 0.5 * alpha.
        Superseded by cases_alpha and controls_alpha if both specified to allow
        for different quantiles for cases vs controls.
    cases_alpha, controls_alpha: Optional[Union[float, Sequence[float]]]
        If specified, separate thresholds for quantiles for cases, controls.
        Both must be specified, otherwise will raise error.
    """

    def __init__(
        self,
        cases: PsiCoverage,
        controls: PsiCoverage,
        controls_min_experiments_f: float = constants.DEFAULT_OUTLIERS_MINEXPERIMENTS,
        alpha: Union[float, Sequence[float]] = constants.DEFAULT_OUTLIERS_ALPHA,
        cases_alpha: Optional[Union[float, Sequence[float]]] = None,
        controls_alpha: Optional[Union[float, Sequence[float]]] = None,
    ):
        # normalize inputs, check for errors
        if cases.num_connections != controls.num_connections:
            raise ValueError(
                "cases and controls must have the same number of connections"
                f" ({cases.num_connections=}, {controls.num_connections=})"
            )
        if prefix_overlap := set(cases.prefixes) & set(controls.prefixes):
            # warning if there is overlap between prefixes
            get_logger().warning(
                f"Heterogen input groups have overlapping prefixes ({prefix_overlap})"
            )
        # get events that passed in enough events in controls
        controls_passed = controls.passed_min_experiments(controls_min_experiments_f)
        # handle alpha
        if not (cases_alpha and controls_alpha):
            if cases_alpha or controls_alpha:
                raise ValueError(
                    "cases_alpha and controls_alpha must either neither be set"
                    f" or both be set ({cases_alpha = }, {controls_alpha = })"
                )
            elif not alpha:
                raise ValueError(
                    "cases_alpha and controls_alpha are not set, so alpha must"
                    f" have a nonzero value ({alpha = })"
                )
            else:  # alpha and not cases_alpha and not controls_alpha
                cases_alpha = alpha
                controls_alpha = alpha
        # make sure cases_alpha, controls_alpha are Sequence[float]
        if isinstance(cases_alpha, float):
            cases_alpha = [cases_alpha]
        if isinstance(controls_alpha, float):
            controls_alpha = [controls_alpha]
        # make sure all values are all reasonable
        if any(x for x in cases_alpha if not (0 < x < 1)):
            raise ValueError(f"cases_alpha must all be in (0, 1) ({cases_alpha = })")
        if any(x for x in controls_alpha if not (0 < x < 1)):
            raise ValueError(
                f"controls_alpha must all be in (0, 1) ({controls_alpha = })"
            )

        # save members of PsiOutliers
        self.cases: Final[PsiCoverage] = cases
        self.controls: Final[PsiCoverage] = controls.mask_events(controls_passed)
        self.controls_passed: Final[xr.DataArray] = controls_passed
        # quantiles to compute
        self.quantiles: Final[xr.Dataset] = xr.Dataset(
            {},
            dict(
                is_lb=[False, True],
                cases_alpha=cases_alpha,
                controls_alpha=controls_alpha,
            ),
        ).assign_coords(
            cases_q=lambda ds: (q := 0.5 * ds.cases_alpha).where(ds.is_lb, 1 - q),
            controls_q=lambda ds: (q := 0.5 * ds.controls_alpha).where(ds.is_lb, 1 - q),
        )
        return

    @cached_property
    def cases_psi_quantile(self) -> xr.DataArray:
        """Cases posterior quantiles with beta approximation of posterior mixture"""
        return self.cases.approximate_quantile(self.quantiles.cases_q)

    @cached_property
    def controls_raw_psi_quantile(self) -> xr.DataArray:
        """Controls population quantiles using raw posterior means"""
        return self.controls.raw_psi_mean_population_quantile(self.quantiles.controls_q)

    @cached_property
    def controls_bootstrap_psi_quantile(self) -> xr.DataArray:
        """Controls population quantiles using bootstrap posterior means"""
        return self.controls.bootstrap_psi_mean_population_quantile(
            self.quantiles.controls_q
        )

    @property
    def controls_raw_psi_median(self) -> xr.DataArray:
        return self.controls.raw_psi_mean_population_median

    @property
    def controls_bootstrap_psi_median(self) -> xr.DataArray:
        return self.controls.bootstrap_psi_mean_population_median
