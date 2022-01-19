"""
PsiCoverage.py

PSI and total coverage (raw and bootstrapped). Converted to/from EventsCoverage.
This allows simplification of the workflow with MOCCASIN bootstrap correction
and more readily parallelized analysis of arbitrarily many files by handling
dependences between junctions.

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import (
    Any,
    Dict,
    Final,
    Hashable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

import dask.array as da
import numpy as np
import numpy.typing as npt
import xarray as xr
from dask.delayed import Delayed
from dask.distributed import progress

import new_majiq._offsets as _offsets
import new_majiq.beta_mixture as bm
import new_majiq.constants as constants
from new_majiq._stats import nanmedian, nanquantile
from new_majiq.experiments import bam_experiment_name
from new_majiq.logger import get_logger

from .Events import Events, _Events
from .EventsCoverage import EventsCoverage
from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions
from .SJExperiment import SJExperiment


def min_experiments(min_experiments_f: float, num_experiments: int) -> float:
    if min_experiments_f < 1:
        min_experiments_f *= num_experiments
    return max(1, min(min_experiments_f, num_experiments))


class PsiCoverage(object):
    """Summarized raw and bootstrap coverage over LSVs for one or more experiments.

    Summarized raw and bootstrap coverage over LSVs for one or more
    experiments as input for quantification.
    Coverage is a total readrate over all bins, excluding stacks and after any
    preceding batch correction steps, ready for quantification.
    Per-experiment coverage stored independently over "prefix" dimension, where
    prefixes originate as the prefix from BAM file names
    (i.e. foo/experiment1.bam -> experiment1).
    Coverage is accompanied by boolean array indicating whether an event is
    "passed for quantification" for each experiment
    (:attr:`PsiCoverage.passed`).

    Provides functionality for combining and summarizing over experiments and
    multiple :class:`PsiCoverage` objects.
    Functions and attributes enable computation of PSI posterior statistics
    under MAJIQ models for splicing quantification.
    Computations are performed over xarray objects.
    When loading :class:`PsiCoverage` from Zarr files, data/computations
    will be loaded/performed lazily using Dask.
    Testing of these computations have been performed over local clusters using
    threads rather than processes (expensive computations generally release the
    GIL).

    Generally, for point estimates of location, quantification with raw coverage
    should be preferred, as bootstrap estimates converge very closely to raw
    estimates as the number of bootstrap replicates becomes large.
    For estimates of variability, quantification with bootstrap coverage should
    be used to account for additional per-bin readrate variability that isn't
    fully captured by the Bayesian model on its own.

    Underlying coverage is stored as the total number of reads over the event
    and the proportion of reads per intron/junction.
    This requires twice the uncompressed memory vs the number of reads per
    intron/junction, but permits easier lazy computation with Dask over large
    datasets.

    Parameters
    ----------
    df: xr.Dataset
        Data variables:
            - event_passed[prefix, ec_idx]
            - raw_total[prefix, ec_idx]
            - raw_psi[prefix, ec_idx]
            - bootstrap_total[prefix, ec_idx, bootstrap_replicate]
            - bootstrap_psi[prefix, ec_idx, bootstrap_replicate]
        Coordinates:
            - lsv_offsets[offset_idx]
            - prefix[prefix]
        Derived (from _offsets):
            - event_size[ec_idx]
            - lsv_idx[ec_idx]
    events: xr.Dataset
        dataset that can be loaded along with matching introns/junctions as
        Events

    See Also
    --------
    PsiCoverage.from_sj_lsvs : Create :class:`PsiCoverage` from :class:`SJExperiment` and :class:`Events`
    PsiCoverage.from_events_coverage : Create :class:`PsiCoverage` from :class:`EventsCoverage`
    PsiCoverage.from_zarr : Load :class:`PsiCoverage` from one or more Zarr files
    PsiCoverage.updated : Create updated :class:`PsiCoverage` with updated arrays
    PsiCoverage.sum : Summed :class:`PsiCoverage` over current prefixes
    PsiCoverage.mask_events : Create updated :class:`PsiCoverage` passing only specified events
    PsiCoverage.__getitem__ : Get :class:`PsiCoverage` for subset of prefixes
    """

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
        """Initialize :class:`PsiCoverage` with specified xarray datasets

        Parameters
        ----------
        df: xr.Dataset
            Data variables:
                - event_passed[ec_idx, prefix]
                - raw_total[ec_idx, prefix]
                - raw_psi[prefix, ec_idx]
                - bootstrap_total[ec_idx, prefix, bootstrap_replicate]
                - bootstrap_psi[ec_idx, prefix, bootstrap_replicate]
            Coordinates:
                - lsv_offsets[offset_idx]
                - prefix[prefix]
            Derived (from _offsets):
                - event_size[ec_idx]
                - lsv_idx[ec_idx]
        events: xr.Dataset
            dataset that can be loaded along with matching introns/junctions as
            Events
        """
        offsets = df["lsv_offsets"].load().values
        if offsets[0] != 0:
            raise ValueError("offsets[0] must be zero")
        if offsets[-1] != df.sizes["ec_idx"]:
            raise ValueError("offsets[-1] must equal number of event connections")
        event_size = np.diff(offsets)
        if "event_size" not in df.variables:
            df = df.assign_coords(
                event_size=("ec_idx", np.repeat(event_size, event_size))
            )
        if "lsv_idx" not in df.variables:
            df = df.assign_coords(
                lsv_idx=("ec_idx", np.repeat(np.arange(len(event_size)), event_size)),
            )
        if df["event_passed"].dtype != bool:
            # for some reason, this is sometimes saved as int8, ensure consistent type
            df["event_passed"] = df["event_passed"].astype(bool)
        self.df: Final[xr.Dataset] = df.transpose(
            "ec_idx", "prefix", ..., "bootstrap_replicate"
        )
        self.events: Final[xr.Dataset] = events
        return

    @property
    def num_connections(self) -> int:
        """Total number of connections over all events where coverage defined"""
        return self.df.sizes["ec_idx"]

    @property
    def num_bootstraps(self) -> int:
        """Number of bootstrap replicates used for bootstraped coverage estimates"""
        return self.df.sizes["bootstrap_replicate"]

    @property
    def num_prefixes(self) -> int:
        """Number of independent units of analysis"""
        return self.df.sizes["prefix"]

    @property
    def prefixes(self) -> List[str]:
        """Prefixes: Names of independent units of analysis

        Names of independent units of analysis (e.g. experiments, aggregations
        over experiments). Name derived from prefixes of BAM files used as
        experiment names.
        """
        return self.df["prefix"].values.tolist()

    def __repr__(self) -> str:
        """String representation of PsiCoverage"""
        MAX_PREFIXES_END = 1  # how many prefixes on either end to display
        if self.num_prefixes > 2 * MAX_PREFIXES_END:
            print_prefixes_list = [
                *self.prefixes[:MAX_PREFIXES_END],
                *([] if self.num_prefixes <= 2 * MAX_PREFIXES_END else ["..."]),
                *self.prefixes[-MAX_PREFIXES_END:],
            ]
        else:
            print_prefixes_list = self.prefixes
        print_prefixes = ", ".join(print_prefixes_list)
        return (
            f"PsiCoverage[{self.num_connections}]"
            f" for {self.num_prefixes} experiments [{print_prefixes}]"
        )

    def __getitem__(self, prefixes) -> "PsiCoverage":
        """Subset :py:class:`PsiCoverage` to selected prefixes"""
        if isinstance(prefixes, str):
            # make sure that prefixes is a sequence
            prefixes = [prefixes]
        return PsiCoverage(self.df.sel(prefix=prefixes), self.events)

    @property
    def event_passed(self) -> xr.DataArray:
        """array(prefix, ec_idx) indicating if event passed"""
        return cast(xr.DataArray, self.df["event_passed"].reset_coords(drop=True))

    @property
    def event_size(self) -> xr.DataArray:
        """array(ec_idx) total number of connections from same event

        Number of connections belonging to event for each event connection,
        which is used to set the parameters for the prior distribution
        """
        return cast(xr.DataArray, self.df["event_size"].reset_coords(drop=True))

    @property
    def lsv_offsets(self) -> xr.DataArray:
        """array(offset_idx) offsets for events into ec_idx"""
        return cast(xr.DataArray, self.df["lsv_offsets"].reset_coords(drop=True))

    @property
    def lsv_idx(self) -> xr.DataArray:
        """array(ec_idx) index identifying event it belongs to"""
        return cast(xr.DataArray, self.df["lsv_idx"].reset_coords(drop=True))

    @property
    def raw_total(self) -> xr.DataArray:
        """array(prefix, ec_idx) raw total reads over event

        py:class:`xr.DataArray` (prefix, ec_idx) raw total reads over event
        (i.e. sum over all event connections per event)
        """
        return cast(xr.DataArray, self.df["raw_total"].reset_coords(drop=True))

    @property
    def raw_psi(self) -> xr.DataArray:
        """array(prefix, ec_idx) percentage of raw_total for connection

        py:class:`xr.DataArray` (prefix, ec_idx) percentage of raw_total for
        connection (maximum likelihood estimate of PSI over raw reads)
        """
        return cast(xr.DataArray, self.df["raw_psi"].reset_coords(drop=True))

    @property
    def bootstrap_total(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) bootstrapped raw_total

        py:class:`xr.DataArray` (prefix, ec_idx, bootstrap_replicate)
        bootstrapped total reads over event (i.e. sum over all event
        connections per event)
        """
        return cast(xr.DataArray, self.df["bootstrap_total"].reset_coords(drop=True))

    @property
    def bootstrap_psi(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) bootstrapped raw_psi

        py:class:`xr.DataArray` (prefix, ec_idx, bootstrap_replicate)
        percentage of bootstrap_total for connection (maximum likelihood
        estimate of PSI over raw reads)
        """
        return cast(xr.DataArray, self.df["bootstrap_psi"].reset_coords(drop=True))

    @cached_property
    def raw_coverage(self) -> xr.DataArray:
        """array(prefix, ec_idx) coverage for individual connection (psi * total)"""
        return self.raw_psi * self.raw_total

    @cached_property
    def bootstrap_coverage(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) bootstrapped raw_coverage"""
        return self.bootstrap_psi * self.bootstrap_total

    @cached_property
    def alpha_prior(self) -> xr.DataArray:
        """array(ec_idx) alpha parameter of prior distribution on PSI for connection"""
        return 1 / self.event_size.astype(self.raw_psi.dtype)

    @cached_property
    def beta_prior(self) -> xr.DataArray:
        """array(ec_idx) beta parameter of prior distribution on PSI for connection"""
        return 1 - self.alpha_prior

    @cached_property
    def raw_alpha(self) -> xr.DataArray:
        """array(prefix, ec_idx) alpha parameter of raw posterior"""
        return (self.raw_coverage + self.alpha_prior).where(self.event_passed)

    @cached_property
    def bootstrap_alpha(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) alpha parameter of bootstrapped posterior"""
        return (self.bootstrap_coverage + self.alpha_prior).where(self.event_passed)

    @cached_property
    def raw_beta(self) -> xr.DataArray:
        """array(prefix, ec_idx) beta parameter of raw posterior"""
        return 1 + self.raw_total - self.raw_alpha

    @cached_property
    def bootstrap_beta(self) -> xr.DataArray:
        """array(prefix, ec_idx, bootstrap_replicate) beta parameter of bootstrapped posterior"""
        return 1 + self.bootstrap_total - self.bootstrap_alpha

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
    def raw_posterior_mean(self) -> xr.DataArray:
        """array(prefix, ec_idx) means of raw posterior distribution on PSI"""
        return self.raw_alpha / np.add(1, self.raw_total)

    @cached_property
    def raw_posterior_variance(self) -> xr.DataArray:
        """array(prefix, ec_idx) variances of raw posterior distribution on PSI"""
        mean = self.raw_posterior_mean
        return mean * np.subtract(1, mean) / np.add(1, self.raw_total)

    @cached_property
    def raw_posterior_std(self) -> xr.DataArray:
        """array(prefix, ec_idx) standard deviations of raw posterior distribution"""
        return cast(xr.DataArray, np.sqrt(self.raw_posterior_variance))

    @property
    def raw_psi_mean(self) -> xr.DataArray:
        """array(prefix, ec_idx) means of raw posterior distribution on PSI

        array(prefix, ec_idx) means of raw posterior distribution on PSI.
        Alias for :py:meth:`PsiCoverage.raw_posterior_mean`
        """
        return self.raw_posterior_mean

    @property
    def raw_psi_variance(self) -> xr.DataArray:
        """array(prefix, ec_idx) variances of raw posterior distribution on PSI

        array(prefix, ec_idx) variances of raw posterior distribution on PSI.
        Alias for :py:meth:`PsiCoverage.raw_posterior_variance`
        """
        return self.raw_posterior_variance

    @property
    def raw_psi_std(self) -> xr.DataArray:
        """array(prefix, ec_idx) standard deviations of raw posterior distribution

        array(prefix, ec_idx) standard deviations of raw posterior distribution
        on PSI. Alias for :py:meth:`PsiCoverage.raw_posterior_std`
        """
        return self.raw_posterior_std

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
        """array(prefix, ec_idx) means of mixtures of bootstrapped posteriors

        array(prefix, ec_idx) means of mixture of bootstrapped posterior
        distribution on PSI
        """
        return self._bootstrap_moments[0]

    @property
    def bootstrap_posterior_variance(self) -> xr.DataArray:
        """array(prefix, ec_idx) variances of mixtures of bootstrapped posteriors

        array(prefix, ec_idx) variances of mixtures of bootstrapped posterior
        distributions on PSI
        """
        return self._bootstrap_moments[1]

    @cached_property
    def bootstrap_posterior_std(self) -> xr.DataArray:
        """array(prefix, ec_idx) standard deviations of mixtures of bootstrapped posteriors

        array(prefix, ec_idx) standard deviations of mixtures of bootstrapped
        posterior distributions on PSI
        """
        return cast(xr.DataArray, np.sqrt(self.bootstrap_posterior_variance))

    @property
    def bootstrap_psi_mean(self) -> xr.DataArray:
        """array(prefix, ec_idx) means of mixtures of bootstrapped posteriors

        array(prefix, ec_idx) means of mixture of bootstrapped posterior
        distribution on PSI.
        Alias for `:py:meth:`PsiCoverage.bootstrap_posterior_mean`
        """
        return self.bootstrap_posterior_mean

    @cached_property
    def bootstrap_psi_mean_legacy(self) -> xr.DataArray:
        """array(prefix, ec_idx) median of means of bootstrapped posteriors

        array(prefix, ec_idx) median of means of bootstrapped posterior
        distributions on PSI.

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

    @property
    def bootstrap_psi_variance(self) -> xr.DataArray:
        """array(prefix, ec_idx) variances of mixtures of bootstrapped posteriors

        array(prefix, ec_idx) variances of mixtures of bootstrapped posterior
        distributions on PSI.
        Alias for `:py:meth:`PsiCoverage.bootstrap_posterior_variance`
        """
        return self.bootstrap_posterior_variance

    @property
    def bootstrap_psi_std(self) -> xr.DataArray:
        """array(prefix, ec_idx) standard deviations of mixtures of bootstrapped posteriors

        array(prefix, ec_idx) standard deviations of mixtures of bootstrapped
        posterior distributions on PSI.
        Alias for `:py:meth:`PsiCoverage.bootstrap_posterior_std`
        """
        return self.bootstrap_posterior_std

    @staticmethod
    def _compute_posterior_quantile(
        a: xr.DataArray,
        b: xr.DataArray,
        quantiles: Union[xr.DataArray, Sequence[float]] = [0.1, 0.9],
        mix_dim: str = "bootstrap_replicate",
    ) -> xr.DataArray:
        if not isinstance(quantiles, xr.DataArray):
            quantiles_arr = np.array(quantiles, dtype=a.dtype)
            if quantiles_arr.ndim > 1:
                raise ValueError(
                    "Unable to handle non-xarray multi-dimensional quantiles"
                )
            elif quantiles_arr.ndim == 0:
                quantiles_arr = quantiles_arr[np.newaxis]
            quantiles = xr.DataArray(quantiles_arr, [("quantiles", quantiles_arr)])
        # if mixture dimension is not present, treat as one-component mixture
        if mix_dim not in a.dims:
            a = a.expand_dims(**{mix_dim: 1})
        if mix_dim not in b.dims:
            b = b.expand_dims(**{mix_dim: 1})
        return xr.apply_ufunc(
            bm.quantile,
            quantiles,
            a,
            b,
            input_core_dims=[[], [mix_dim], [mix_dim]],
            dask="allowed",
        )

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
            Return array(ec_idx, ...) of quantiles per connection. If
            `quantiles` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "quantiles"

        Notes
        -----
        Please use :py:meth:`PsiCoverage.approximate_quantile` instead, which
        is faster, and what we think is a better representation of PSI
        """
        return self._compute_posterior_quantile(
            self.bootstrap_alpha, self.bootstrap_beta, quantiles=quantiles
        )

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
            Return array(ec_idx, ...) of quantiles per connection. If
            `quantiles` is not :py:class:`xr.DataArray`, dimension over
            quantiles will be "quantiles"

        See Also
        --------
        :py:meth:`PsiCoverage.bootstrap_quantile`
        """
        return self._compute_posterior_quantile(
            self.approximate_alpha, self.approximate_beta, quantiles=quantiles
        )

    @staticmethod
    def _compute_posterior_discretized_pmf(
        a: xr.DataArray,
        b: xr.DataArray,
        nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
        mix_dim: str = "bootstrap_replicate",
    ) -> xr.DataArray:
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
        return self._compute_posterior_discretized_pmf(
            self.bootstrap_alpha, self.bootstrap_beta, nbins=nbins
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
        return self._compute_posterior_discretized_pmf(
            self.approximate_alpha, self.approximate_beta, nbins=nbins
        )

    @cached_property
    def _raw_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        result = self.raw_psi_mean
        if result.chunks:
            result = result.chunk({"prefix": None})
        return result

    @cached_property
    def _bootstrap_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        result = self.bootstrap_psi_mean
        if result.chunks:
            result = result.chunk({"prefix": None})
        return result

    @cached_property
    def raw_psi_mean_population_median(self) -> xr.DataArray:
        """array(ec_idx) median over prefixes of :py:meth:`PsiCoverage.raw_psi_mean`"""
        return xr.apply_ufunc(
            nanmedian,
            self._raw_psi_mean_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    @cached_property
    def bootstrap_psi_mean_population_median(self) -> xr.DataArray:
        """array(ec_idx) median over prefixes of :py:meth:`PsiCoverage.bootstrap_psi_mean`"""
        return xr.apply_ufunc(
            nanmedian,
            self._bootstrap_psi_mean_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    @staticmethod
    def _compute_population_quantile(
        x: xr.DataArray,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
    ) -> xr.DataArray:
        quantiles_xr = xr.DataArray(quantiles, [(quantile_dim_name, quantiles)])
        return xr.apply_ufunc(
            nanquantile,
            x,
            quantiles_xr,
            input_core_dims=[["prefix"], [quantile_dim_name]],
            output_core_dims=[[quantile_dim_name]],
            dask="allowed",
        )

    def raw_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
    ) -> xr.DataArray:
        """empirical quantiles over prefixes of :py:meth:`PsiCoverage.raw_psi_mean`

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles over quantified population to compute
        quantiles_dim_name: str
            Name of dimension in output array matching `quantiles`

        Returns
        -------
        xr.DataArray
            array(ec_idx, `quantiles_dim_name`) of quantiles per connection
            over quantified prefixes
        """
        return self._compute_population_quantile(
            self._raw_psi_mean_core_prefix,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )

    def bootstrap_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
    ) -> xr.DataArray:
        """empirical quantiles over prefixes of :py:meth:`PsiCoverage.bootstrap_psi_mean`

        Parameters
        ----------
        quantiles: Sequence[float]
            quantiles over quantified population to compute
        quantiles_dim_name: str
            Name of dimension in output array matching `quantiles`

        Returns
        -------
        xr.DataArray
            array(ec_idx, `quantiles_dim_name`) of quantiles per connection
            over quantified prefixes
        """
        return self._compute_population_quantile(
            self._bootstrap_psi_mean_core_prefix,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )

    @classmethod
    def from_events_coverage(
        cls,
        events_coverage: EventsCoverage,
        minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
        minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
    ) -> "PsiCoverage":
        """Create :py:class:`PsiCoverage` from :py:class:`EventsCoverage`

        Parameters
        ----------
        minreads, minbins: float
            Quantifiability thresholds

        Returns
        -------
        PsiCoverage
        """
        # get offsets as int (not uint)
        offsets: npt.NDArray[np.int64] = np.array(
            events_coverage.events._offsets, dtype=np.int64
        )
        # get whether individual connection passes thresholds
        passed = (events_coverage.numreads >= minreads) & (
            events_coverage.numbins >= minbins
        )
        # get whether any connection in event passed, per connection
        event_passed = _offsets.offset_logical_or(passed, offsets)
        # get total coverage per event, per connection
        raw_total = _offsets.offsetsum(
            events_coverage.numreads, offsets, axes=[0, -1, 0]
        )
        bootstrap_total = _offsets.offsetsum(
            events_coverage.bootstraps, offsets, axes=[0, -1, 0]
        )
        # get psi per connection
        with np.errstate(divide="ignore", invalid="ignore"):
            raw_psi = np.where(raw_total > 0, events_coverage.numreads / raw_total, 0)
            bootstrap_psi = np.where(
                bootstrap_total > 0, events_coverage.bootstraps / bootstrap_total, 0
            )
        # return dataset with matched values
        return cls(
            xr.Dataset(
                data_vars=dict(
                    event_passed=("ec_idx", event_passed),
                    raw_total=("ec_idx", raw_total),
                    raw_psi=("ec_idx", raw_psi),
                    bootstrap_total=(
                        ("ec_idx", "bootstrap_replicate"),
                        bootstrap_total,
                    ),
                    bootstrap_psi=(("ec_idx", "bootstrap_replicate"), bootstrap_psi),
                ),
                coords=dict(
                    lsv_offsets=("offset_idx", offsets),
                ),
                attrs=dict(
                    minreads=minreads,
                    minbins=minbins,
                    bam_path=events_coverage.bam_path,
                    bam_version=events_coverage.bam_version,
                ),
            ).expand_dims(prefix=[bam_experiment_name(events_coverage.bam_path)]),
            events_coverage.events.save_df,
        )

    @classmethod
    def from_sj_lsvs(
        cls,
        sj: SJExperiment,
        lsvs: Events,
        minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
        minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
        num_bootstraps: int = constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        pvalue_threshold: float = constants.DEFAULT_COVERAGE_STACK_PVALUE,
    ) -> "PsiCoverage":
        """Create :class:`PsiCoverage` from :class:`SJExperiment` and :class:`Events`.

        Parameters
        ----------
        sj: SJExperiment
            Intron and junction coverage from an experiment
        lsvs: Events
            Events over which PsiCoverage will be defined
        minreads, minbins: float
            Quantifiability thresholds
        num_bootstraps: int
            The number of bootstrap replicates for bootstrapped estimates
        pvalue_threshold: float
            P-value threshold for removing stacks under leave-one-out Poisson
            model of per-bin read coverage, for both raw and bootstrapped
            coverage (Set to nonpositive value to skip stack detection)

        Returns
        -------
        PsiCoverage

        Notes
        -----
        The pvalue_threshold for stack removal is applied to both raw and
        bootstrapped coverage. This differs from the behavior in MAJIQ v2,
        where stack removal was only applied to bootstrapped coverage.
        In this sense "raw" coverage is only after stack detection.
        """
        lsv_coverage = EventsCoverage.from_events_and_sj(
            lsvs, sj, num_bootstraps=num_bootstraps, pvalue_threshold=pvalue_threshold
        )
        return PsiCoverage.from_events_coverage(lsv_coverage, minreads, minbins)

    @classmethod
    def from_zarr(cls, path: Union[str, Path, List[Union[str, Path]]]) -> "PsiCoverage":
        """Load :py:class:`PsiCoverage` from one or more specified paths

        Load :py:class:`PsiCoverage` from one or more specified paths.
        Prefixes will be concatenated (overlapping prefixes will use values
        from the first file it is found in).

        Parameters
        ----------
        path: Union[str, Path, List[Union[str, Path]]]
            path or paths with PsiCoverage saved in zarr format

        Returns
        -------
        PsiCoverage
            PsiCoverage for prefixes found in all specified files

        Notes
        -----
        Does not check that events are same in each input file. It will fail if
        they are not the same size, which should catch most cases, but be wary
        that events information is derived from the first file alone.
        """
        if not isinstance(path, list):
            path = [path]
        df = xr.open_mfdataset(
            path,
            engine="zarr",
            group=constants.NC_PSICOVERAGE,
            combine="nested",
            concat_dim="prefix",
            join="override",
            compat="override",
            coords="minimal",
            data_vars="minimal",
        )
        if len(path) > 1:
            # attributes are defined by path[0]. We'd rather just have none
            df.attrs.clear()
        events_df = xr.open_zarr(path[0], group=constants.NC_EVENTS)
        return cls(df, events_df)

    def updated(
        self,
        bootstrap_psi: Optional[xr.DataArray],
        raw_psi: Optional[xr.DataArray],
        **update_attrs,
    ) -> "PsiCoverage":
        """Create updated :py:class:`PsiCoverage` with new values of psi

        Parameters
        ----------
        bootstrap_psi, raw_psi: Optional[xr.DataArray]
            If specified, new values of bootstrap_psi, raw_psi to use in new
            PsiCoverage (with all other variables equal)
        update_attrs:
            Additional kwargs are set as attributes to the dataset used to
            construct the resulting PsiCoverage

        Returns
        -------
        PsiCoverage
            Updated :py:class:`PsiCoverage` with new values of psi
        """
        df = self.df
        # update psi arrays
        if bootstrap_psi is not None:
            if set(self.bootstrap_psi.dims) != set(bootstrap_psi.dims):
                raise ValueError("bootstrap_psi doesn't have same named axes")
            df = df.assign(bootstrap_psi=bootstrap_psi)
        if raw_psi is not None:
            if set(self.raw_psi.dims) != set(raw_psi.dims):
                raise ValueError("raw_psi doesn't have same named axes")
            df = df.assign(raw_psi=raw_psi)
        # update/add attributes
        df = df.assign_attrs(**update_attrs)
        # return resulting PsiCoverage object
        return PsiCoverage(df, self.events)

    def _save_df(
        self,
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        remove_bam_attrs: bool = False,
    ) -> xr.Dataset:
        """Prepare dataset of psicoverage that will be saved

        Prepare dataset of psicoverage that will be saved. This sets the
        chunking, clears any encodings from before, removes attrs that have to
        do with BAM, etc.
        """
        USE_CHUNKS = dict(
            ec_idx=ec_chunksize, bootstrap_replicate=None, offset_idx=None
        )
        save_df = self.df.drop_vars(["event_size", "lsv_idx"])
        if save_df.sizes["ec_idx"] > 0:
            save_df = save_df.chunk(USE_CHUNKS)  # type: ignore[arg-type]
        # clear any previous encodings from before, before saving
        for v in save_df.variables.values():
            v.encoding.clear()
        # remove attributes referring to bam?
        if remove_bam_attrs:
            for k in [x for x in save_df.attrs.keys() if str(x).startswith("bam")]:
                save_df.attrs.pop(k, None)
        return save_df

    def to_zarr(
        self,
        path: Union[str, Path],
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        consolidated: bool = True,
        show_progress: bool = False,
    ) -> None:
        """Save :py:class:`PsiCoverage` to specified path

        Parameters
        ----------
        path: Union[str, Path]
            Path for output file in zarr format
        ec_chunksize: int
            How to chunk event connections to prevent memory from getting to
            large when loading many samples simultaneously
        consolidated: bool
            When saving the file make sure that it is consolidated. In general,
            if you are appending a bunch of files together, it can make sense
            to set consolidated=False, and consolidate on the last write (only
            consolidate once). But, don't forget to consolidate at the end.
        show_progress: bool
            Attempt to show progress on distributed cluster for Dask
        """
        save_df = self._save_df(ec_chunksize=ec_chunksize)
        save_df_future = cast(
            Delayed,
            save_df.to_zarr(
                path,
                mode="w",
                group=constants.NC_PSICOVERAGE,
                consolidated=False,
                compute=False,
            ),
        )
        if show_progress:
            save_df_future = save_df_future.persist()
            progress(save_df_future)
        else:
            save_df_future.compute()
        self.events.chunk(self.events.sizes).to_zarr(
            path, mode="a", group=constants.NC_EVENTS, consolidated=consolidated
        )
        return

    def to_zarr_slice(
        self,
        path: Union[str, Path],
        prefix_slice: slice,
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
    ) -> None:
        """Save :py:class:`PsiCoverage` to specified path for specified slice on prefix

        Save :py:class:`PsiCoverage` to specified path for specified slice.
        Typically run after :py:meth:`PsiCoverage.to_zarr_slice_init`

        Parameters
        ----------
        path: Union[str, Path]
            Path with output file in zarr format with metadata initialized by
            :py:meth:`PsiCoverage.to_zarr_slice_init`
        prefix_slice: slice
            Slice of prefix dimension in output zarr store to save current
            PsiCoverage
        ec_chunksize: int
            How to chunk event connections to prevent memory from getting to
            large when loading many samples simultaneously

        See Also
        --------
        PsiCoverage.to_zarr_slice_init
        """
        self._save_df(ec_chunksize=ec_chunksize).drop_vars("prefix").pipe(
            lambda x: x.drop_vars(
                [k for k, v in x.variables.items() if "prefix" not in v.dims]
            )
        ).to_zarr(
            path, group=constants.NC_PSICOVERAGE, region=dict(prefix=prefix_slice)
        )
        return

    @classmethod
    def to_zarr_slice_init(
        cls,
        path: Union[str, Path],
        events_df: xr.Dataset,
        prefixes: List[str],
        num_bootstraps: int,
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        cov_dtype: type = np.float32,
        psicov_attrs: Dict[Hashable, Any] = dict(),
    ) -> None:
        """Initialize zarr store for saving :py:class:`PsiCoverage` over many writes

        Initialize zarr for :py:class:`PsiCoverage` over many prefixes.
        Saves all information except dimensions that are prefix-specific.
        This enables multithreaded (or multiprocess) write with
        :py:meth:`PsiCoverage.to_zarr_slice`

        Parameters
        ----------
        path: Union[str, Path]
            Path for output Zarr for psicoverage output
        events_df: xr.Dataset
            Dataset encoding psicoverage events (Events.save_df,
            PsiCoverage.events, etc.)
        prefixes: List[str]
            Values for the prefix dimension coordinate
        num_bootstraps: int
            Number of bootstrap replicates that will be used
        ec_chunksize: int
            How to chunk event connections to prevent memory from getting to
            large when loading many samples simultaneously
        cov_dtype: type
            What type to use for psi/total_coverage arrays
        psicov_attrs: Dict[Hashable, Any]
            Attributes to include

        See Also
        --------
        PsiCoverage.to_zarr_slice
        """
        # force events to be saved as single chunk (no benefit for chunking here)
        events_df = events_df.chunk(events_df.sizes)
        # save events
        events_df.to_zarr(path, mode="w", group=constants.NC_EVENTS, consolidated=False)
        # dims for skeleton
        raw_dims = ("ec_idx", "prefix")
        bootstrap_dims = (*raw_dims, "bootstrap_replicate")
        # shapes for skeleton
        raw_shape = (events_df.sizes["ec_idx"], len(prefixes))
        bootstrap_shape = (*raw_shape, num_bootstraps)
        # chunksizes for skeleton
        raw_chunks = (ec_chunksize, 1)
        bootstrap_chunks = (*raw_chunks, None)
        # arrays for skeleton
        raw_arr = da.empty(raw_shape, dtype=cov_dtype, chunks=raw_chunks)
        passed_arr = da.empty(raw_shape, dtype=bool, chunks=raw_chunks)
        bootstrap_arr = da.empty(
            bootstrap_shape, dtype=cov_dtype, chunks=bootstrap_chunks
        )
        # save metadata for skeleton
        xr.Dataset(
            dict(
                bootstrap_psi=(bootstrap_dims, bootstrap_arr),
                bootstrap_total=(bootstrap_dims, bootstrap_arr),
                event_passed=(raw_dims, passed_arr),
                raw_psi=(raw_dims, raw_arr),
                raw_total=(raw_dims, raw_arr),
            ),
        ).to_zarr(
            path,
            mode="a",
            compute=False,
            group=constants.NC_PSICOVERAGE,
            consolidated=False,
        )
        # save offsets, prefixes, and attributes
        add_offsets = (
            events_df[["_offsets"]]
            .reset_coords()
            .astype(int)
            .set_coords("_offsets")
            .rename_dims(e_offsets_idx="offset_idx")
            .rename_vars(_offsets="lsv_offsets")
            .assign_coords(prefix=("prefix", prefixes))
        )
        add_offsets.attrs = psicov_attrs  # overwrite attrs with what we want
        add_offsets.to_zarr(
            path, mode="a", group=constants.NC_PSICOVERAGE, consolidated=True
        )
        return

    @classmethod
    def convert_sj_batch(
        cls,
        sjs: Sequence[Path],
        lsvs: Events,
        path: Path,
        minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
        minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
        num_bootstraps: int = constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        pvalue_threshold: float = constants.DEFAULT_COVERAGE_STACK_PVALUE,
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        imap_unordered_fn=map,
    ) -> None:
        """Load PsiCoverage from sj paths, save to single output path

        Parameters
        ----------
        sjs: Sequence[Path]
            Paths to input SJ files that will have PsiCoverage evaluated for
        lsvs: Events
            Events over which PsiCoverage will be calculated
        path: Path
            Output path for PsiCoverage zarr file
        minreads, minbins: float
            Quantifiability thresholds
        num_bootstraps: int
            The number of bootstrap replicates for bootstrapped estimates
        pvalue_threshold: float
            P-value threshold for removing stacks under leave-one-out Poisson
            model of per-bin read coverage, for both raw and bootstrapped
            coverage (Set to nonpositive value to skip stack detection)
        imap_unordered_fn
            Loading/saving of input SJ files will be passed through this
            function, which can enable concurrency if appropriate
        """
        log = get_logger()

        # how to go from sj path to psicoverage file
        def sj_to_psicov(sj_path: Path) -> PsiCoverage:
            return PsiCoverage.from_sj_lsvs(
                SJExperiment.from_zarr(sj_path),
                lsvs,
                minreads=minreads,
                minbins=minbins,
                num_bootstraps=num_bootstraps,
                pvalue_threshold=pvalue_threshold,
            )

        if len(sjs) == 0:
            raise ValueError("At least one SJ file must be processed")
        elif len(sjs) == 1:
            # if there is only one file, don't bother with imap_unordered_fn
            log.info("Inferring PsiCoverage from %s", sjs[0])
            psi_coverage = sj_to_psicov(sjs[0])
            log.info("Saving %s to %s", psi_coverage, path)
            psi_coverage.to_zarr(path, ec_chunksize=ec_chunksize)
        else:
            # precompute prefixes to use
            log.debug("Precomputing prefixes corresponding to input SJ files")
            prefixes = [
                bam_experiment_name(SJExperiment.original_path_from_zarr(x))
                for x in sjs
            ]
            # we have more than one input file
            log.info("Saving event information and metadata to %s", path)
            PsiCoverage.to_zarr_slice_init(
                path,
                lsvs.save_df,
                prefixes,
                num_bootstraps,
                psicov_attrs=dict(
                    sj=[str(x) for x in sjs],
                    minreads=minreads,
                    minbins=minbins,
                ),
                ec_chunksize=ec_chunksize,
            )

            def job_fn(sj_idx: int, sj: Path) -> Path:
                sj_to_psicov(sj).to_zarr_slice(
                    path,
                    slice(sj_idx, 1 + sj_idx),
                    ec_chunksize=ec_chunksize,
                )
                return sj

            jobs = imap_unordered_fn(lambda x: job_fn(x[0], x[1]), list(enumerate(sjs)))
            for ndx, sj in enumerate(jobs, 1):
                log.info("Saved coverage from %s (%d / %d)", sj, ndx, len(sjs))
        return

    @cached_property
    def num_passed(self) -> xr.DataArray:
        """array(ec_idx) number of experiments passed for each connection"""
        return self.event_passed.sum("prefix")

    def passed_min_experiments(
        self,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
    ) -> xr.DataArray:
        """Return array(ec_idx) boolean mask for events passing min_experiments

        Parameters
        ----------
        min_experiments_f: float
            Threshold for group filters. This specifies the fraction (value <
            1) or absolute number (value >= 1) of prefixes that must pass
            individually for the event to be considered as passed for the group

        Returns
        -------
        xr.DataArray
            array(ec_idx) indicate whether the event (for each event connection)
            passed min_experiments
        """
        return self.num_passed >= min_experiments(min_experiments_f, self.num_prefixes)

    def sum(
        self,
        new_prefix: str,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
    ) -> "PsiCoverage":
        """Create aggregated :py:class:`PsiCoverage` with sum coverage over prefixes

        Parameters
        ----------
        new_prefix: str
            Prefix for summarized :py:class:`PsiCoverage`
        min_experiments_f: float
            Threshold for group filters. This specifies the fraction (value <
            1) or absolute number (value >= 1) of prefixes that must pass
            individually for the event to be considered as passed for the group

        Returns
        -------
        PsiCoverage
            Sum coverage over prefixes, with passed being defined over group
            filters
        """
        if self.num_prefixes > 1:
            event_passed = self.passed_min_experiments(min_experiments_f)
            raw_total = self.raw_total.sum("prefix")
            raw_coverage = (self.raw_total * self.raw_psi).sum("prefix")
            raw_psi = (raw_coverage / raw_total.where(raw_total > 0)).fillna(0)
            bootstrap_total = self.bootstrap_total.sum("prefix")
            bootstrap_coverage = (self.bootstrap_total * self.bootstrap_psi).sum(
                "prefix"
            )
            bootstrap_psi = (
                bootstrap_coverage / bootstrap_total.where(bootstrap_total > 0)
            ).fillna(0)
            df = xr.Dataset(
                data_vars=dict(
                    event_passed=event_passed,
                    raw_total=raw_total,
                    raw_psi=raw_psi,
                    bootstrap_total=bootstrap_total,
                    bootstrap_psi=bootstrap_psi,
                ),
                coords=dict(
                    lsv_offsets=self.lsv_offsets,
                    event_size=self.event_size,
                    lsv_idx=self.lsv_idx,
                ),
                attrs=dict(original_prefix=self.prefixes),
            ).expand_dims(prefix=[new_prefix])
        else:
            df = self.df.assign_coords(prefix=[new_prefix])
        return PsiCoverage(df, self.events)

    def mask_events(self, passed: xr.DataArray) -> "PsiCoverage":
        """Return :py:class:`PsiCoverage` passing only events that are passed in input

        Parameters
        ----------
        passed: xr.DataArray
            boolean array(ec_idx) where connections marked as not passed
            (False) will be marked as not passed (False) in resulting
            :py:class:`PsiCoverage`
        Returns
        -------
        PsiCoverage
            Same coverage but passing only events that are passed in input (and
            in the original object)
        """
        return PsiCoverage(
            self.df.assign(event_passed=self.event_passed & passed), self.events
        )

    def get_events(
        self,
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> Events:
        """Construct :py:class:`Events` using saved dataset and introns, junctions

        Parameters
        ----------
        introns: GeneIntrons
        junctions: GeneJunctions

        Returns
        -------
        Events
        """
        if self.events.intron_hash != introns.checksum():
            raise ValueError("GeneIntrons checksums do not match")
        if self.events.junction_hash != junctions.checksum():
            raise ValueError("GeneJunctions checksums do not match")
        return Events(
            _Events(
                introns._gene_introns,
                junctions._gene_junctions,
                self.events.ref_exon_idx,
                self.events.event_type,
                self.events._offsets,
                self.events.is_intron,
                self.events.connection_idx,
            )
        )

    def dataset(
        self,
        properties: Sequence[str] = constants.DEFAULT_PSI_PROPERTIES,
        quantiles: Sequence[float] = list(),
        psibins: Optional[int] = None,
    ) -> xr.Dataset:
        """Extract selected properties into single :py:class:`xr.Dataset`

        Parameters
        ----------
        properties: Sequence[str]
            PsiCoverage properties to request.
        quantiles: Sequence[float]
            If non-empty, calculate quantiles of posterior distribution
        psibins: Optional[int]
            If specified, calculate discretized approximation to posterior
            distribution with this many bins

        Notes
        -----
        quantiles and psibins computations use the approximate posterior
        distribution Beta(approximate_alpha, approximate_beta) because:

        - by default, it's 30 times faster than using the bootstrap mixture,
        - for high coverage events, the bootstrap distribution is discrete (for
          each bootstrap replicate), so the approximate distribution is a
          better representation of our desired model of variability in PSI.
        """
        # initialize variables to return with noting if any experiment passed
        quantify_vars: Dict[str, xr.DataArray] = {
            "any_passed": self.event_passed.any("prefix")
        }
        # add properties
        for x in properties:
            quantify_vars[x] = getattr(self, x)
        if len(quantiles):
            quantify_vars["psi_quantile"] = self.approximate_quantile(quantiles)
        if psibins:
            quantify_vars["psi_pmf"] = self.approximate_discretized_pmf(psibins)
        return xr.Dataset(quantify_vars).reset_coords(drop=True)  # type: ignore[arg-type]
