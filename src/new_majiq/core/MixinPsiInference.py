"""
MixinPsiInference.py

Abstract mixin classes to define how different PSI quantities are computed
from alpha, beta, etc.

Author: Joseph K Aicher
"""

from abc import ABC, abstractmethod
from functools import cached_property
from typing import Sequence, Tuple, Union, cast

import numpy as np
import xarray as xr

import new_majiq.beta_mixture as bm
import new_majiq.constants as constants
from new_majiq._stats import nanmedian, nanquantile


class MixinRawPsi(ABC):
    """Methods for PSI inference when raw_alpha, raw_beta are defined"""

    @property
    @abstractmethod
    def raw_alpha(self) -> xr.DataArray:
        """array(...) of alpha for raw posterior"""
        ...

    @property
    @abstractmethod
    def raw_beta(self) -> xr.DataArray:
        """array(...) of beta for raw posterior"""
        ...

    @cached_property
    def raw_alpha_plus_beta(self) -> xr.DataArray:
        return self.raw_alpha + self.raw_beta

    @cached_property
    def raw_posterior_mean(self) -> xr.DataArray:
        """array(...) means of raw posterior distribution on PSI"""
        return self.raw_alpha / (self.raw_alpha_plus_beta)

    @cached_property
    def raw_posterior_variance(self) -> xr.DataArray:
        """array(...) variances of raw posterior distribution on PSI"""
        mean = self.raw_posterior_mean
        return mean * np.subtract(1, mean) / np.add(1, self.raw_alpha_plus_beta)

    @cached_property
    def raw_posterior_std(self) -> xr.DataArray:
        """array(...) standard deviations of raw posterior distribution"""
        return cast(xr.DataArray, np.sqrt(self.raw_posterior_variance))

    @property
    def raw_psi_mean(self) -> xr.DataArray:
        """array(...) means of raw posterior distribution on PSI (alias)"""
        return self.raw_posterior_mean

    @property
    def raw_psi_variance(self) -> xr.DataArray:
        """array(...) variances of raw posterior distribution on PSI (alias)"""
        return self.raw_posterior_variance

    @property
    def raw_psi_std(self) -> xr.DataArray:
        """array(...) standard deviations of raw posterior distribution (alias)"""
        return self.raw_posterior_std


class MixinApproximatePsi(ABC):
    """Methods for PSI inference when approximate_alpha, approximate_beta defined"""

    @property
    @abstractmethod
    def approximate_alpha(self) -> xr.DataArray:
        """array(...) of alpha for smooth posterior (approximation of bootstrap mixture)"""
        ...

    @property
    @abstractmethod
    def approximate_beta(self) -> xr.DataArray:
        """array(...) of beta for smooth posterior (approximation of bootstrap mixture)"""
        ...

    @cached_property
    def approximate_alpha_plus_beta(self) -> xr.DataArray:
        return self.approximate_alpha + self.approximate_beta

    @cached_property
    def bootstrap_posterior_mean(self) -> xr.DataArray:
        """array(...) means of bootstrap posterior distribution on PSI"""
        return self.approximate_alpha / self.approximate_alpha_plus_beta

    @cached_property
    def bootstrap_posterior_variance(self) -> xr.DataArray:
        """array(...) variances of bootstrap posterior distribution on PSI"""
        mean = self.bootstrap_posterior_mean
        return mean * np.subtract(1, mean) / np.add(1, self.approximate_alpha_plus_beta)

    @cached_property
    def bootstrap_posterior_std(self) -> xr.DataArray:
        """array(...) standard deviations of bootstrap posterior distribution"""
        return cast(xr.DataArray, np.sqrt(self.bootstrap_posterior_variance))

    @property
    def bootstrap_psi_mean(self) -> xr.DataArray:
        """array(...) means of bootstrap posterior distribution on PSI (alias)"""
        return self.bootstrap_posterior_mean

    @property
    def bootstrap_psi_variance(self) -> xr.DataArray:
        """array(...) variances of bootstrap posterior distribution on PSI (alias)"""
        return self.bootstrap_posterior_variance

    @property
    def bootstrap_psi_std(self) -> xr.DataArray:
        """array(...) standard deviations of bootstrap posterior distribution (alias)"""
        return self.bootstrap_posterior_std

    def approximate_cdf(
        self,
        x: Union[xr.DataArray, Sequence[float]],
    ) -> xr.DataArray:
        """Compute cdf of approximate/smoothed bootstrapped posterior

        Parameters
        ----------
        x: Union[xr.DataArray, Sequence[float]]
            potential realizations of distribution to evaluate CDF at

        Returns
        -------
        xr.DataArray
            Return array(...) of cdf probabilities per connection. If
            `x` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "x"
        """
        return _compute_posterior_cdf(self.approximate_alpha, self.approximate_beta, x)

    def approximate_quantile(
        self,
        quantiles: Union[xr.DataArray, Sequence[float]] = [0.1, 0.9],
    ) -> xr.DataArray:
        """Compute quantiles of approximate/smoothed bootstrapped posterior

        Parameters
        ----------
        quantiles: Union[xr.DataArray, Sequence[float]]
            quantiles of distribution to compute

        Returns
        -------
        xr.DataArray
            Return array(...) of quantiles per connection. If
            `quantiles` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "quantiles"
        """
        return _compute_posterior_quantile(
            self.approximate_alpha, self.approximate_beta, quantiles=quantiles
        )

    def approximate_discretized_pmf(
        self, nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS
    ) -> xr.DataArray:
        """Compute discretized PMF of approximate/smoothed bootstrap posterior

        Parameters
        ----------
        nbins: int
            Number of uniform bins on [0, 1] on which probability mass will be
            computed
        """
        return _compute_posterior_discretized_pmf(
            self.approximate_alpha, self.approximate_beta, nbins=nbins
        )

    def approximate_updf(
        self,
        nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
    ) -> xr.DataArray:
        """Compute unnormalized PDF of approximate/smoothed bootstrap posterior

        Parameters
        ----------
        nbins: int
            Compute PDF over endpoints of uniformly spaced bins on [0, 1].
            (the first and last values are computed at the midpoints in order
            to handle singularities at {0, 1} when either of the beta
            distribution parameters are less than 1).

        Notes
        -----
        This is appropriate for plotting because

        - it is much faster than computing the PMF
        - the usual plots are qualitative with arbitrary scale, so the
          normalization constant is irrelevant
        """
        return _compute_beta_updf(self.approximate_alpha, self.approximate_beta)


class MixinBootstrapPsi(MixinApproximatePsi, ABC):
    """Methods for PSI inference when bootstrap_alpha, bootstrap_beta defined"""

    @property
    @abstractmethod
    def bootstrap_alpha(self) -> xr.DataArray:
        """array(..., bootstrap_replicate) of alpha for posterior bootstrap replicates"""
        ...

    @property
    @abstractmethod
    def bootstrap_beta(self) -> xr.DataArray:
        """array(..., bootstrap_replicate) of beta for posterior bootstrap replicates"""
        ...

    @cached_property
    def _approximate_params(self) -> Tuple[xr.DataArray, xr.DataArray]:
        """Beta distribution parameters matching mean, variance of bootstrap mixture

        In many cases, we operate on the bootstrapped distributions as a single
        distribution by treating it as a uniform mixture over bootstrap
        replicates.
        This mixture is an poorly-behaved model for fixed number of bootstrap replicates as the total coverage increases (the bootstrap replicates behave as atoms).
        This motivates making a smooth approximation by a single beta
        distribution.
        Here, we approximate the beta mixture by matching mean and variance,
        which we prefer in most cases.
        """
        _params = xr.apply_ufunc(
            bm.approximation,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            output_core_dims=[["_param_idx"]],
            dask="parallelized",
            dask_gufunc_kwargs=dict(output_sizes={"_param_idx": 2}),
        )
        alpha = _params.isel(_param_idx=0)
        beta = _params.isel(_param_idx=1)
        return (alpha, beta)

    @property
    def approximate_alpha(self) -> xr.DataArray:
        """array(prefix, ec_idx) alpha parameter of approximated bootstrap posterior

        In many cases, we operate on the bootstrapped distributions as a single
        distribution by treating it as a uniform mixture over bootstrap
        replicates.
        This mixture is an poorly-behaved model for fixed number of bootstrap replicates as the total coverage increases (the bootstrap replicates behave as atoms).
        This motivates making a smooth approximation by a single beta
        distribution.
        Here, we approximate the beta mixture by matching mean and variance,
        which we prefer in most cases.
        """
        return self._approximate_params[0]

    @property
    def approximate_beta(self) -> xr.DataArray:
        """array(prefix, ec_idx) beta parameter of approximated bootstrap posterior

        In many cases, we operate on the bootstrapped distributions as a single
        distribution by treating it as a uniform mixture over bootstrap
        replicates.
        This mixture is an poorly-behaved model for fixed number of bootstrap replicates as the total coverage increases (the bootstrap replicates behave as atoms).
        This motivates making a smooth approximation by a single beta
        distribution.
        Here, we approximate the beta mixture by matching mean and variance,
        which we prefer in most cases.
        """
        return self._approximate_params[1]

    @cached_property
    def _bootstrap_moments(self) -> Tuple[xr.DataArray, xr.DataArray]:
        """Compute mean, variance of bootstrap posterior mixture"""
        _moments = xr.apply_ufunc(
            bm.moments,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            output_core_dims=[["_moment_idx"]],
            dask="parallelized",
            dask_gufunc_kwargs=dict(output_sizes={"_moment_idx": 2}),
        )
        agg_mean = _moments.isel(_moment_idx=0)
        agg_variance = _moments.isel(_moment_idx=1)
        return (agg_mean, agg_variance)

    @property
    def bootstrap_posterior_mean(self) -> xr.DataArray:
        """array(...) means of mixtures of bootstrapped posteriors

        array(...) means of mixture of bootstrapped posterior distribution on PSI
        """
        return self._bootstrap_moments[0]

    @property
    def bootstrap_posterior_variance(self) -> xr.DataArray:
        """array(...) variances of mixtures of bootstrapped posteriors

        array(...) variances of mixtures of bootstrapped posterior
        distributions on PSI
        """
        return self._bootstrap_moments[1]

    @cached_property
    def bootstrap_psi_mean_legacy(self) -> xr.DataArray:
        """array(...) median of means of bootstrapped posteriors

        array(...) median of means of bootstrapped posterior distributions on PSI.

        Notes
        -----
        This is what was reported in MAJIQ v1 and v2.
        We have observed that if we increase the number of bootstrap replicates,
        both estimates tend close (but not exactly) to the raw posterior mean,
        which we now prefer.
        """
        return xr.apply_ufunc(
            bm.means_median,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            dask="parallelized",
        )

    def bootstrap_cdf(
        self,
        x: Union[xr.DataArray, Sequence[float]],
    ) -> xr.DataArray:
        """Compute cdf of mixture of bootstrapped posterior distribution

        Parameters
        ----------
        x: Union[xr.DataArray, Sequence[float]]
            potential realizations of distribution to evaluate CDF at

        Returns
        -------
        xr.DataArray
            Return array(...) of CDF probabilities per connection. If
            `x` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "x"
        """
        return _compute_posterior_cdf(self.bootstrap_alpha, self.bootstrap_beta, x)

    def bootstrap_quantile(
        self,
        quantiles: Union[xr.DataArray, Sequence[float]] = [0.1, 0.9],
    ) -> xr.DataArray:
        """Compute quantiles of mixture of bootstrapped posterior distributions

        Parameters
        ----------
        quantiles: Union[xr.DataArray, Sequence[float]]
            quantiles of distribution to compute

        Returns
        -------
        xr.DataArray
            Return array(...) of quantiles per connection. If
            `quantiles` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "quantiles"

        Notes
        -----
        Please use `approximate_quantile` instead, which is faster, and what we
        think is a better representation of PSI
        """
        return _compute_posterior_quantile(
            self.bootstrap_alpha, self.bootstrap_beta, quantiles=quantiles
        )

    def bootstrap_discretized_pmf(
        self, nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS
    ) -> xr.DataArray:
        """Compute discretized PMF of bootstrap posterior mixture

        Parameters
        ----------
        nbins: int
            Number of uniform bins on [0, 1] on which probability mass will be
            computed
        """
        return _compute_posterior_discretized_pmf(
            self.bootstrap_alpha, self.bootstrap_beta, nbins=nbins
        )


class MixinRawPsiMeanPopulation(ABC):
    """Methods for summarizing raw_psi_mean over a population (dimension prefix)"""

    @property
    @abstractmethod
    def raw_psi_mean(self):
        """array(...) means of raw posterior distribution on PSI"""
        ...

    @cached_property
    def _raw_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        result = self.raw_psi_mean
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": None})
        return result

    @cached_property
    def raw_psi_mean_population_median(self) -> xr.DataArray:
        """array(...) median over prefixes of `raw_psi_mean`"""
        return xr.apply_ufunc(
            nanmedian,
            self._raw_psi_mean_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    def raw_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
    ) -> xr.DataArray:
        """empirical quantiles over prefixes of `raw_psi_mean`

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles over quantified population to compute
        quantiles_dim_name: str
            Name of dimension in output array matching `quantiles`

        Returns
        -------
        xr.DataArray
            array(..., `quantiles_dim_name`) of quantiles per connection
            over quantified prefixes
        """
        return _compute_population_quantile(
            self._raw_psi_mean_core_prefix,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )


class MixinBootstrapPsiMeanPopulation(ABC):
    """Methods for summarizing bootstrap_psi_mean over a population (dimension prefix)"""

    @property
    @abstractmethod
    def bootstrap_psi_mean(self):
        """array(...) means of bootstrap posterior distribution on PSI"""
        ...

    @cached_property
    def _bootstrap_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        result = self.bootstrap_psi_mean
        if "prefix" not in result.dims:
            raise ValueError(
                "Missing required dimension 'prefix' for population summaries"
            )
        if result.chunks:
            result = result.chunk({"prefix": None})
        return result

    @cached_property
    def bootstrap_psi_mean_population_median(self) -> xr.DataArray:
        """array(...) median over prefixes of `bootstrap_psi_mean`"""
        return xr.apply_ufunc(
            nanmedian,
            self._bootstrap_psi_mean_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    def bootstrap_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
    ) -> xr.DataArray:
        """empirical quantiles over prefixes of `bootstrap_psi_mean`

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles over quantified population to compute
        quantiles_dim_name: str
            Name of dimension in output array matching `quantiles`

        Returns
        -------
        xr.DataArray
            array(..., `quantiles_dim_name`) of quantiles per connection
            over quantified prefixes
        """
        return _compute_population_quantile(
            self._bootstrap_psi_mean_core_prefix,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )


def _compute_population_quantile(
    x: xr.DataArray,
    quantiles: Sequence[float],
    quantile_dim_name: str = "population_quantile",
) -> xr.DataArray:
    """Compute quantiles (ignoring nan) over dimension prefix"""
    quantiles_xr = xr.DataArray(quantiles, [(quantile_dim_name, quantiles)])
    return xr.apply_ufunc(
        nanquantile,
        x,
        quantiles_xr,
        input_core_dims=[["prefix"], [quantile_dim_name]],
        output_core_dims=[[quantile_dim_name]],
        dask="allowed",
    )


def _compute_posterior_cdf(
    a: xr.DataArray,
    b: xr.DataArray,
    x: Union[xr.DataArray, Sequence[float]],
    mix_dim: str = "bootstrap_replicate",
) -> xr.DataArray:
    """Caclulate cdf over posterior distribution

    Parameters
    ----------
    a, b: xr.DataArray
        Parameters of posterior distributions
    x: Union[xr.DataArray, Sequence[float]]
        potential realizations of the distribution at which to evaluate the CDF
    mix_dim: str
        Dimension of `a` and `b` over which distribution is a mixture. If not
        found in a or b, treat as a single beta distribution
    """
    if not isinstance(x, xr.DataArray):
        x_arr = np.array(x, dtype=a.dtype)
        if x_arr.ndim > 1:
            raise ValueError("Unable to handle non-xarray multi-dimensional x")
        elif x_arr.ndim == 0:
            x_arr = x_arr[np.newaxis]
        x = xr.DataArray(x_arr, [("x", x_arr)])
    # if mixture dimension is not present, treat as one-component mixture
    if mix_dim not in a.dims:
        a = a.expand_dims(**{mix_dim: 1})
    if mix_dim not in b.dims:
        b = b.expand_dims(**{mix_dim: 1})
    if a.sizes[mix_dim] != b.sizes[mix_dim]:
        raise ValueError(
            f"a and b must share same size for dimension {mix_dim}"
            f" ({a.sizes = }, {b.sizes = })"
        )
    return xr.apply_ufunc(
        bm.cdf,
        x,
        a,
        b,
        input_core_dims=[[], [mix_dim], [mix_dim]],
        dask="parallelized",
    )


def _compute_posterior_quantile(
    a: xr.DataArray,
    b: xr.DataArray,
    quantiles: Union[xr.DataArray, Sequence[float]] = [0.1, 0.9],
    mix_dim: str = "bootstrap_replicate",
) -> xr.DataArray:
    """Calculate quantiles over posterior distribution

    Parameters
    ----------
    a, b: xr.DataArray
        Parameters of posterior distributions
    quantiles: Union[xr.DataArray, Sequence[float]]
        quantiles to compute over posterior distribution. If Sequence[float],
        add "quantiles" dimension to result for input quantiles
    mix_dim: str
        Dimension of `a` and `b` over which distribution is a mixture. If not
        found in a or b, treat as a single beta distribution
    """
    if not isinstance(quantiles, xr.DataArray):
        quantiles_arr = np.array(quantiles, dtype=a.dtype)
        if quantiles_arr.ndim > 1:
            raise ValueError("Unable to handle non-xarray multi-dimensional quantiles")
        elif quantiles_arr.ndim == 0:
            quantiles_arr = quantiles_arr[np.newaxis]
        quantiles = xr.DataArray(quantiles_arr, [("quantiles", quantiles_arr)])
    # if mixture dimension is not present, treat as one-component mixture
    if mix_dim not in a.dims:
        a = a.expand_dims(**{mix_dim: 1})
    if mix_dim not in b.dims:
        b = b.expand_dims(**{mix_dim: 1})
    if a.sizes[mix_dim] != b.sizes[mix_dim]:
        raise ValueError(
            f"a and b must share same size for dimension {mix_dim}"
            f" ({a.sizes = }, {b.sizes = })"
        )
    return xr.apply_ufunc(
        bm.quantile,
        quantiles,
        a,
        b,
        input_core_dims=[[], [mix_dim], [mix_dim]],
        dask="allowed",
    )


def _compute_posterior_discretized_pmf(
    a: xr.DataArray,
    b: xr.DataArray,
    nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
    mix_dim: str = "bootstrap_replicate",
) -> xr.DataArray:
    """Compute discretized PMF of posterior (mixture)

    Parameters
    ----------
    a, b: xr.DataArray
        Parameters of posterior distributions
        Dimensions: [..., ?mix_dim]
    nbins: int
        Number of uniform bins on [0, 1] on which probability mass will be
        computed
    mix_dim: str
        Dimension of `a` and `b` over which distribution is a mixture. If not
        found in a or b, treat as a single beta distribution

    Returns
    -------
    xr.DataArray
        discretized PMF for posterior (dimensions: [..., pmf_bin])
    """
    endpoints = np.linspace(0, 1, 1 + nbins, dtype=a.dtype)
    dummy_bins = xr.DataArray(
        np.empty(nbins, dtype=a.dtype),
        {
            "pmf_bin_start": ("pmf_bin", endpoints[:-1]),
            "pmf_bin_end": ("pmf_bin", endpoints[1:]),
        },
        dims=["pmf_bin"],
    )
    # if mixture dimension is not present, treat as one-component mixture
    if mix_dim not in a.dims:
        a = a.expand_dims(**{mix_dim: 1})
    if mix_dim not in b.dims:
        b = b.expand_dims(**{mix_dim: 1})
    return xr.apply_ufunc(
        bm.pmf,
        a,
        b,
        dummy_bins,
        input_core_dims=[
            [mix_dim],
            [mix_dim],
            ["pmf_bin"],
        ],
        output_core_dims=[["pmf_bin"]],
        dask="allowed",
    )


def _compute_beta_updf(
    a: xr.DataArray,
    b: xr.DataArray,
    nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
) -> xr.DataArray:
    """Compute unnormalized PDF of beta distribution

    Parameters
    ----------
    nbins: int
        Compute PDF over endpoints of uniformly spaced bins on [0, 1].
        (the first and last values are computed at the midpoints in order
        to handle singularities at {0, 1} when either of the beta
        distribution parameters are less than 1).

    Notes
    -----
    This is appropriate for plotting because

    - it is much faster than computing the PMF
    - the usual plots are qualitative with arbitrary scale, so the
      normalization constant is irrelevant
    """
    if nbins < 2:
        raise ValueError("approximate_updf requires at least 2 bins")
    psi_arr = np.linspace(0, 1, 1 + nbins, dtype=a.dtype)
    psi_calc = psi_arr.copy()
    psi_calc[0] = 0.5 * psi_calc[1]
    psi_calc[-1] = 0.5 + 0.5 * psi_calc[-2]
    psi = xr.DataArray(psi_calc, [("psi", psi_arr)], name="psi")
    logpsi = np.log(psi)
    log1mpsi = np.log1p(-psi)
    pdf = np.exp((a - 1) * logpsi + (b - 1) * log1mpsi)
    return pdf
