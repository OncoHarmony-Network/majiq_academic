"""
DeltaPsi.py

Quantify differences between two groups of replicate experiments

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import TYPE_CHECKING, Final, cast

import numpy as np
import xarray as xr
from scipy.special import logsumexp

import new_majiq.beta_mixture as bm
import new_majiq.constants as constants

from .DPsiPrior import DPsiPrior
from .PMFSummaries import PMFSummaries
from .PsiCoverage import PsiCoverage

if TYPE_CHECKING:
    from ..voila.DeltaPsiDataset import DeltaPsiDataset


class DeltaPsi(object):
    """Compute DeltaPsi between two groups of PsiCoverage (replicate assumption)

    Compute DeltaPsi between two groups of PsiCoverage under assumption that
    the experiments in each group are replicates with the given DPsiPrior.

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
        if psi1.num_connections != psi2.num_connections:
            raise ValueError(
                "psi1 and psi2 must have the same number of connections"
                f" ({psi1.num_connections=}, {psi2.num_connections=})"
            )
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

    @property
    def num_connections(self) -> int:
        return self.psi1.num_connections

    @property
    def name1(self) -> str:
        return self.psi1.prefixes[0]

    @property
    def name2(self) -> str:
        return self.psi2.prefixes[0]

    def __repr__(self) -> str:
        return (
            f"DeltaPsi[{self.num_connections}] for {self.name1} vs {self.name2}"
            f" computed with {self.psibins} psi bins and {self.prior}"
        )

    def rebin(self, psibins: int) -> "DeltaPsi":
        """Get DeltaPsi with different bins"""
        return DeltaPsi(
            self.psi1,
            self.psi2,
            self.prior,
            psibins=psibins,
            name1=self.name1,
            name2=self.name2,
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
        a1 = (
            self.psi1.approximate_alpha.squeeze("prefix", drop=True)
            .expand_dims(mix1=1)
            .where(self.passed)
        )
        b1 = (
            self.psi1.approximate_beta.squeeze("prefix", drop=True)
            .expand_dims(mix1=1)
            .where(self.passed)
        )
        a2 = (
            self.psi2.approximate_alpha.squeeze("prefix", drop=True)
            .expand_dims(mix2=1)
            .where(self.passed)
        )
        b2 = (
            self.psi2.approximate_beta.squeeze("prefix", drop=True)
            .expand_dims(mix2=1)
            .where(self.passed)
        )
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
    def discrete_bootstrap_logposterior(self) -> xr.DataArray:
        """Average bootstrap replicates after inference of deltapsi"""
        a1 = self.psi1.bootstrap_alpha.rename(prefix="mix1").where(self.passed)
        b1 = self.psi1.bootstrap_beta.rename(prefix="mix1").where(self.passed)
        a2 = self.psi2.bootstrap_alpha.rename(prefix="mix2").where(self.passed)
        b2 = self.psi2.bootstrap_beta.rename(prefix="mix2").where(self.passed)
        per_bootstrap = xr.apply_ufunc(
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
        # take sum over bootstrap replicates, which isn't normalized
        unnormalized = xr.apply_ufunc(
            logsumexp,
            per_bootstrap,
            input_core_dims=[["bootstrap_replicate"]],
            kwargs={"axis": -1},
            dask="parallelized",
        )
        # normalize each resulting distribution
        logZ = xr.apply_ufunc(
            logsumexp,
            unnormalized,
            input_core_dims=[["pmf_bin"]],
            kwargs={"axis": -1},
            dask="parallelized",
        )
        return unnormalized / logZ

    @cached_property
    def discrete_posterior(self) -> PMFSummaries:
        """:class:`PMFSummaries` for dpsi posterior from approximated PSI posteriors"""
        return PMFSummaries(cast(xr.DataArray, np.exp(self.discrete_logposterior)))

    @cached_property
    def discrete_bootstrap_posterior(self) -> PMFSummaries:
        """:class:`PMFSummaries` for average bootstrapped dpsi posteriors"""
        return PMFSummaries(
            cast(xr.DataArray, np.exp(self.discrete_bootstrap_logposterior))
        )

    @property
    def discrete_posterior_mean(self) -> xr.DataArray:
        return self.discrete_posterior.mean

    @property
    def discrete_posterior_variance(self) -> xr.DataArray:
        return self.discrete_posterior.variance

    @property
    def discrete_posterior_std(self) -> xr.DataArray:
        return cast(xr.DataArray, np.sqrt(self.discrete_posterior_variance))

    @property
    def discrete_bootstrap_posterior_mean(self) -> xr.DataArray:
        """Posterior mean on deltapsi from discrete_bootstrap_posterior"""
        return self.discrete_bootstrap_posterior.mean

    @property
    def discrete_bootstrap_posterior_variance(self) -> xr.DataArray:
        """Posterior variance on deltapsi from discrete_bootstrap_posterior"""
        return self.discrete_bootstrap_posterior.variance

    @property
    def discrete_bootstrap_posterior_std(self) -> xr.DataArray:
        """Posterior standard deviaion on deltapsi from discrete_bootstrap_posterior"""
        return cast(xr.DataArray, np.sqrt(self.discrete_bootstrap_posterior_variance))

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

    def bootstrap_probability_changing(
        self, changing_threshold: float = constants.DEFAULT_DPSI_CHANGING_THRESHOLD
    ):
        return 1 - self.discrete_bootstrap_posterior.interval_probability(
            -changing_threshold, changing_threshold
        )

    def bootstrap_probability_nonchanging(
        self,
        nonchanging_threshold: float = constants.DEFAULT_DPSI_NONCHANGING_THRESHOLD,
    ):
        return self.discrete_bootstrap_posterior.interval_probability(
            -nonchanging_threshold, nonchanging_threshold
        )

    @property
    def dataset(self) -> "DeltaPsiDataset":
        """Reduce to :class:`DeltaPsiDataset` used for VOILA visualization"""
        from ..voila.DeltaPsiDataset import DeltaPsiDataset

        return DeltaPsiDataset.from_deltapsi(self)
