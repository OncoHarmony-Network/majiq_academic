"""
DPsiPrior.py

Define prior on deltapsi using input PsiCoverage

Author: Joseph K Aicher
"""

from typing import List, Optional, Union

import numpy as np
import xarray as xr
from scipy.special import logsumexp
from scipy.stats import beta as beta_dist

import new_majiq.constants as constants
from new_majiq.logger import get_logger
from new_majiq.PsiCoverage import PsiCoverage, min_experiments


def get_empirical_dpsi(
    psi1: PsiCoverage,
    psi2: PsiCoverage,
    minreads: float = constants.DEFAULT_DPSI_PRIOR_MINREADS,
    min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
) -> xr.DataArray:
    """Get high confidence empirical deltapsi from input groups of experiments

    Get high confidence empirical deltapsi from input groups of experiments
    meeting criteria:
    + one junction only per event
    + binary events only
    + must have passed at least min_experiments events
    + must also pass additional minreads criteria (higher confidence, etc.)

    Parameters
    ----------
    psi1, psi2: PsiCoverage
        psi coverage for two groups of experiments. Read evidence will be
        combined per group.
    minreads: float
        Additional criteria (beyond having previously passed) for being
        considered
    min_experiments_f: float
        Proportion (if < 1) or number of experiments that must pass all
        criteria in order to be considered

    Returns
    -------
    xr.DataArray
        difference in raw_psi between sum of two groups for selected events
        Dimensions: lsv_idx
    """
    reads1 = psi1.raw_total * psi1.raw_psi
    passed1 = (psi1.event_passed & (reads1 >= minreads)).sum(
        "prefix"
    ) >= min_experiments(min_experiments_f, psi1.df.sizes["prefix"])
    reads2 = psi2.raw_total * psi2.raw_psi
    passed2 = (psi2.event_passed & (reads2 >= minreads)).sum(
        "prefix"
    ) >= min_experiments(min_experiments_f, psi2.df.sizes["prefix"])
    # mask includes potentially duplicate events
    passed = passed1 & passed2 & (psi1.event_size == 2)
    # get dpsi that passed, keeping only one per lsv
    dpsi = (
        (
            (reads2.sum("prefix") / psi2.raw_total.sum("prefix"))
            - (reads1.sum("prefix") / psi1.raw_total.sum("prefix"))
        )
        .where(passed)
        .load()
        .dropna("ec_idx")
        .groupby("lsv_idx")
        .first()
    )
    return dpsi


def infer_pmix_given_dpsi(
    dpsi: xr.DataArray,
    a: xr.DataArray,
    pmix: xr.DataArray,
) -> xr.DataArray:
    """Get probability of membership to mixture given observaion (E-step)

    Parameters
    ----------
    dpsi: xr.DataArray (dimension: lsv_idx)
    a, pmix: xr.DataArray (dimension: mixture_component)

    Returns
    -------
    pmix_given_dpsi: xr.DataArray (dimensions: lsv_idx, mixture_component)
    """
    if (a < 1).any():
        raise ValueError(
            "Prior may not have distinct modes at endpoints (reuire a >= 1)"
        )
    if (pmix <= 0).any():
        raise ValueError("prior miture probabilities must be positive")
    # get result in logspace
    likelihood = xr.apply_ufunc(beta_dist.logpdf, dpsi, a, a, -1, 2) + np.log(pmix)
    # normalize per lsv_idx
    logZ = xr.apply_ufunc(
        logsumexp,
        likelihood,
        input_core_dims=[["mixture_component"]],
        kwargs=dict(axis=-1),
    )
    return np.exp(likelihood - logZ)


def fit_pmix(pmix_given_dpsi: xr.DataArray) -> xr.DataArray:
    """Fit pmix using pmix_given_dpsi (M-step)"""
    x = pmix_given_dpsi.sum("lsv_idx")
    return x / x.sum()  # make sure result is normalized


def fit_a(
    dpsi: xr.DataArray, pmix_given_dpsi: xr.DataArray, force_slab: bool = True
) -> xr.DataArray:
    """Fit a using dpsi, pmix_given_dpsi (M-step) by method of moments"""
    variance = xr.dot(
        dpsi, dpsi, pmix_given_dpsi, dims="lsv_idx"
    ) / pmix_given_dpsi.sum("lsv_idx")
    # get beta parameter that has this variance (when scaled on [-1, 1])
    a = 0.5 * (1 / variance - 1)
    if force_slab:
        # force the first value of a to always be uniform distribution
        a.loc[dict(mixture_component=0)] = 1.0
    # return result
    return a


def discrete_log_prior(
    pmix: xr.DataArray, a: xr.DataArray, psibins: int
) -> xr.DataArray:
    """Get discretized logprior for deltapsi with 2 * psibins deltapsi bins"""
    endpoints = xr.DataArray(
        np.linspace(-1, 1, 1 + 2 * psibins),
        dims="pmf_bin",
    )
    cdf = xr.dot(xr.apply_ufunc(beta_dist.cdf, endpoints, a, a, -1, 2), pmix)
    pmf = (
        cdf.isel(pmf_bin=slice(1, None)) - cdf.isel(pmf_bin=slice(None, -1))
    ).assign_coords(
        pmf_bin_start=endpoints.isel(pmf_bin=slice(None, -1)),
        pmf_bin_end=endpoints.isel(pmf_bin=slice(1, None)),
    )
    PSEUDO = 1e-20
    return np.log(PSEUDO + pmf)


def fit_discrete_dpsi_prior(
    dpsi: xr.DataArray,
    a: Union[List[float], xr.DataArray] = constants.DEFAULT_DPSI_PRIOR_A,
    pmix: Union[List[float], xr.DataArray] = constants.DEFAULT_DPSI_PRIOR_PMIX,
    min_lsvs: int = constants.DEFAULT_DPSI_PRIOR_MINLSV,
    n_update_a: int = constants.DEFAULT_DPSI_PRIOR_MAXITER,
    n_update_pmix: Optional[int] = None,  # None --> 1 + n_update_a
    psibins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
) -> xr.DataArray:
    """Given observed dpsi, fit mixture parameters, compute discretized bins

    Parameters
    ----------
    dpsi: xr.DataArray (dimension: lsv_idx)
    a, pmix: Union[List[float], xr.DataArray (dimension: mixture_component)]
        Default parameters for prior. Must have same length.
    min_lsvs: int
        If less than this many LSVs in dpsi, don't fit, use default/initial
        prior
    n_update_a: int
        Number of iterations to update a during M step
    n_update_pmix: Optional[int]
        Number of iterations to update pmix. If none, use 1 + n_update_a.
    psibins: int
        Number of bins for PSI (2 * psibins for deltapsi)
    """
    # normalize inputs to be usable
    if not isinstance(a, xr.DataArray):
        a = xr.DataArray(a, dims=["mixture_component"])
    if not isinstance(pmix, xr.DataArray):
        pmix = xr.DataArray(pmix, dims=["mixture_component"])
    if n_update_pmix is None:
        n_update_pmix = 1 + n_update_a
    # get logger
    log = get_logger()
    # are we doing EM iteration to update prior?
    n_update = max(0, n_update_a, n_update_pmix)
    if n_update and dpsi.sizes["lsv_idx"] < min_lsvs:
        log.info(
            f"Only {dpsi.sizes['lsv_idx']} dpsi LSVs provided to fit prior"
            f" (min_lsvs={min_lsvs}). Using default prior"
        )
        n_update = 0
    if n_update:
        log.info(
            "Updating initial deltapsi prior with parameters"
            f" (a={[f'{x:.2e}' for x in a.values]},"
            f" pmix={[f'{x:.2e}' for x in pmix.values]})"
        )
    for i in range(n_update):
        # E step
        pmix_given_dpsi = infer_pmix_given_dpsi(dpsi, a, pmix)
        # M step
        if i < n_update_pmix:
            pmix = fit_pmix(pmix_given_dpsi)
        if i < n_update_a:
            a = fit_a(dpsi, pmix_given_dpsi)
    log.info(
        f"Using deltapsi prior with parameters (a={[f'{x:.2e}' for x in a.values]},"
        f" pmix={[f'{x:.2e}' for x in pmix.values]})"
    )
    return discrete_log_prior(pmix, a, psibins)
