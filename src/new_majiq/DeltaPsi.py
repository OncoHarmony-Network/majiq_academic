"""
DeltaPsi.py

Quantify differences between two groups of replicate experiments

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import Final, cast

import numpy as np
import xarray as xr

import new_majiq.beta_mixture as bm
import new_majiq.constants as constants
from new_majiq.DPsiPrior import DPsiPrior
from new_majiq.PMFSummaries import PMFSummaries
from new_majiq.PsiCoverage import PsiCoverage


class DeltaPsi(object):
    """DeltaPsi between two groups of replicate experiments

    Parameters
    ----------
    psi1, psi2: PsiCoverage
        PsiCoverage files to be treated as replicates
    prior: DPsiPrior
        Prior on deltapsi for inference
    min_experiments_f: float
        Only pass events that pass this many experiments in both groups
    name1, name2: str
        Names to use for the aggregated groups coverage/PSI
    """

    def __init__(
        self,
        psi1: PsiCoverage,
        psi2: PsiCoverage,
        prior: DPsiPrior,
        psibins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        name1: str = "psi1",
        name2: str = "psi2",
    ):
        # save aggregate psi coverage, original prior
        self.psi1: Final[PsiCoverage] = psi1.sum(
            name1, min_experiments_f=min_experiments_f
        )
        self.psi2: Final[PsiCoverage] = psi2.sum(
            name2, min_experiments_f=min_experiments_f
        )
        self.prior: Final[DPsiPrior] = prior
        self.psibins: Final[int] = psibins
        return

    def rebin(self, psibins: int) -> "DeltaPsi":
        """Get DeltaPsi with different bins"""
        return DeltaPsi(
            self.psi1,
            self.psi2,
            self.prior,
            psibins=psibins,
            name1=self.psi1.prefixes[0],
            name2=self.psi2.prefixes[0],
        )

    @cached_property
    def passed(self) -> xr.DataArray:
        return self.psi1.event_passed.squeeze(
            "prefix", drop=True
        ) & self.psi2.event_passed.squeeze("prefix", drop=True)

    @cached_property
    def discrete_logprior(self) -> xr.DataArray:
        return self.prior.discretized_logpmf(psibins=self.psibins)

    @cached_property
    def discrete_logposterior(self) -> xr.DataArray:
        a1 = self.psi1.approximate_alpha.squeeze("prefix", drop=True).expand_dims(
            mix1=1
        )
        b1 = self.psi1.approximate_beta.squeeze("prefix", drop=True).expand_dims(mix1=1)
        a2 = self.psi2.approximate_alpha.squeeze("prefix", drop=True).expand_dims(
            mix2=1
        )
        b2 = self.psi2.approximate_beta.squeeze("prefix", drop=True).expand_dims(mix2=1)
        return xr.apply_ufunc(
            bm.dpsi_discrete,
            a1,
            b1,
            a2,
            b2,
            self.discrete_logprior,
            input_core_dims=[["mix1"], ["mix1"], ["mix2"], ["mix2"], ["pmf_bin"]],
            output_core_dims=[["pmf_bin"]],
            dask="parallelized",
        )

    @cached_property
    def discrete_posterior(self) -> PMFSummaries:
        return PMFSummaries(cast(xr.DataArray, np.exp(self.discrete_logposterior)))

    @property
    def discrete_posterior_mean(self) -> xr.DataArray:
        return self.discrete_posterior.mean

    @property
    def discrete_posterior_variance(self) -> xr.DataArray:
        return self.discrete_posterior.variance

    @property
    def discrete_posterior_std(self) -> xr.DataArray:
        return cast(xr.DataArray, np.sqrt(self.discrete_posterior_variance))

    def probability_changing(
        self, changing_threshold: float = constants.DEFAULT_DPSI_CHANGING_THRESHOLD
    ):
        return 1 - self.discrete_posterior.interval_probability(
            -changing_threshold, changing_threshold
        )

    def probability_nonchanging(
        self,
        nonchanging_threshold: float = constants.DEFAULT_DPSI_NONCHANGING_THRESHOLD,
    ):
        return self.discrete_posterior.interval_probability(
            -nonchanging_threshold, nonchanging_threshold
        )

    def dataset(
        self,
        changing_threshold: float = constants.DEFAULT_DPSI_CHANGING_THRESHOLD,
        nonchanging_threshold: float = constants.DEFAULT_DPSI_NONCHANGING_THRESHOLD,
    ) -> xr.Dataset:
        # get PSI datasets
        psi1_ds = (
            self.psi1.dataset()
            .drop_vars("any_passed")
            .squeeze("prefix", drop=True)
            .pipe(
                lambda ds: ds.rename(
                    {k: f"{self.psi1.prefixes[0]}_{k}" for k in ds.variables.keys()}
                )
            )
        )
        psi2_ds = (
            self.psi2.dataset()
            .drop_vars("any_passed")
            .squeeze("prefix", drop=True)
            .pipe(
                lambda ds: ds.rename(
                    {k: f"{self.psi2.prefixes[0]}_{k}" for k in ds.variables.keys()}
                )
            )
        )
        # get deltapsi dataset
        dpsi_ds = xr.Dataset(
            {
                "passed": self.passed,
                "dpsi_mean": self.discrete_posterior_mean,
                "dpsi_std": self.discrete_posterior_std,
                "probability_changing": self.probability_changing(changing_threshold),
                "probability_nonchanging": self.probability_nonchanging(
                    nonchanging_threshold
                ),
            }
        ).reset_coords(drop=True)
        return xr.merge([dpsi_ds, psi1_ds, psi2_ds], compat="override", join="exact")
