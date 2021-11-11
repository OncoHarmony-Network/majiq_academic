"""
PsiOutliers.py

Identify aberrant splicing events of interests in cases vs controls

Author: Joseph K Aicher
"""

from typing import Final, Optional, Sequence, Union

import numpy as np
import xarray as xr

import new_majiq.constants as constants
from new_majiq.logger import get_logger
from new_majiq.PsiControlsSummary import (
    PsiControlsSummary,
    _psirange_from_psiquantiles,
    _q_from_alpha,
)
from new_majiq.PsiCoverage import PsiCoverage


class PsiOutliers(object):
    """Outliers in PSI between cases and controls

    Parameters
    ----------
    cases: PsiCoverage
        PsiCoverage for prefixes that will be treated as cases
    controls: PsiControlsSummary
        Summary of PsiCoverage for prefixes that are treated as controls

    Notes
    -----
    Use PsiOutliers.from_psicov() if you haven't already summarized controls
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

    @classmethod
    def from_psicov(
        cls,
        cases: PsiCoverage,
        controls: PsiCoverage,
        controls_alpha: Union[
            float, Sequence[float]
        ] = constants.DEFAULT_OUTLIERS_ALPHA,
    ) -> "PsiOutliers":
        """Outliers in PSI between cases and controls

        Parameters
        ----------
        cases: PsiCoverage
            PsiCoverage for prefixes that will be treated as cases
        controls: PsiCoverage
            PsiCoverage for prefixes that will be treated as controls
        controls_alpha: Union[float, Sequence[float]]
            threshold used for controls quantiles in both directions (motivated
            by two-sided tests interpretation on CDFs for null-hypothesis
            testing). Quantiles are 0.5 * alpha, 1 - 0.5 * alpha
        """
        controls_summary = PsiControlsSummary.from_psicov(
            controls, alpha=controls_alpha
        )
        return PsiOutliers(cases, controls_summary)

    @property
    def num_connections(self) -> int:
        return self.cases.num_connections

    @property
    def cases_psi_mean(self) -> xr.DataArray:
        """Cases raw posterior means"""
        return self.cases.raw_psi_mean

    def evaluate_cases_alpha(
        self,
        cases_alpha: Optional[Union[float, Sequence[float]]] = None,
        min_experiments_f: Union[
            float, Sequence[float]
        ] = constants.DEFAULT_OUTLIERS_MINEXPERIMENTS,
    ) -> xr.Dataset:
        """Compare ranges from cases/controls

        Parameters
        ----------
        cases_alpha: Optional[Union[float, Sequence[float]]]
            Specify values of alpha (determines two-sided quantiles for each
            case) to use. If None, use the same values of alpha used for
            controls.
        min_experiments_f: Union[float, Sequence[float]]
            Identify if controls has at least this many experiments that passed
            quantification thresholds contributing to psi quantiles. If less
            than 1, proportion of total number of control experiments.
        """
        if cases_alpha is None:
            cases_alpha = self.controls.alpha.values.tolist()
        elif isinstance(cases_alpha, float):
            cases_alpha = [cases_alpha]
        else:
            # make sure sorted and has unique elements
            cases_alpha = sorted(set(cases_alpha))
        ds = (
            xr.Dataset(
                dict(
                    # controls data_vars
                    controls_psi_median=self.controls.psi_median,
                    controls_psi_quantile=self.controls.psi_quantile,
                    controls_psi_range=self.controls.psi_range,
                    controls_passed=self.controls.passed_min_experiments(
                        min_experiments_f
                    ),
                    # cases data_vars
                    cases_psi_mean=self.cases_psi_mean,
                ),
                # controls coordinates implicitly passed through above data_vars
                dict(
                    # add cases_q, cases_alpha
                    cases_q=_q_from_alpha(
                        xr.DataArray(
                            cases_alpha,
                            [("cases_alpha", cases_alpha)],
                            name="cases_alpha",
                        )
                    )
                ),
                dict(controls_prefixes=self.controls.prefixes),
            )
            .assign(
                # compute psi quantiles
                cases_psi_quantile=lambda ds: self.cases.approximate_quantile(
                    ds["cases_q"]
                )
            )
            .assign(
                dpsi_quantile_gap=lambda ds: self._dpsi_quantile_gap(
                    ds["cases_psi_quantile"], ds["controls_psi_quantile"]
                ),
                cases_psi_range=lambda ds: _psirange_from_psiquantiles(
                    ds["cases_psi_quantile"]
                ),
            )
            .assign(
                combined_psi_range=lambda ds: self._combined_psi_range(
                    ds["cases_psi_range"], ds["controls_psi_range"]
                )
            )
            .assign(
                dpsi_quantile_gap_rescaled=lambda ds: (
                    ds["dpsi_quantile_gap"] / ds["combined_psi_range"]
                )
            )
        )
        return ds

    @staticmethod
    def _dpsi_quantile_gap(
        cases_psi_quantile: xr.DataArray, controls_psi_quantile: xr.DataArray
    ) -> xr.DataArray:
        # gap in quantiles when case is higher/lower than controls
        gap_case_higher = cases_psi_quantile.sel(
            is_lb=True
        ) - controls_psi_quantile.sel(is_lb=False)
        gap_case_lower = controls_psi_quantile.sel(is_lb=True) - cases_psi_quantile.sel(
            is_lb=False
        )
        # get the one value which is positive (or zero, if ranges overlap)
        return (
            gap_case_higher.clip(min=gap_case_lower)
            .clip(min=0)
            .rename("dpsi_quantile_gap")
        )

    @staticmethod
    def _combined_psi_range(
        cases_psi_range: xr.DataArray, controls_psi_range: xr.DataArray
    ) -> xr.DataArray:
        """treat ranges as independent/orthogonal length scales"""
        return np.sqrt(np.square(cases_psi_range) + np.square(controls_psi_range))
