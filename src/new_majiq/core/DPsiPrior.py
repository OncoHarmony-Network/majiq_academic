"""
DPsiPrior.py

Define prior on deltapsi using input PsiCoverage

Author: Joseph K Aicher
"""

from typing import Final, List, Optional, Union, cast

import numpy as np
import xarray as xr
from dask.distributed import progress
from scipy.special import logsumexp
from scipy.stats import beta as beta_dist

import new_majiq.constants as constants
from new_majiq.logger import get_logger

from .PsiCoverage import PsiCoverage, min_experiments


class DPsiPrior(object):
    """Prior on DeltaPsi as weighted mixture of beta distributions (over [-1, 1])

    Parameters
    ----------
    a, pmix: Union[List[float], xr.DataArray]
        Default parameters for prior. Must have same length.
        If xr.DataArray, must have dimension mixture_component.
        a is the parameters for each component beta.
        pmix is the probability of each component.
    """

    def __init__(
        self,
        a: Union[List[float], xr.DataArray] = constants.DEFAULT_DPSI_PRIOR_A,
        pmix: Union[List[float], xr.DataArray] = constants.DEFAULT_DPSI_PRIOR_PMIX,
    ):
        # normalize inputs to be usable
        if not isinstance(a, xr.DataArray):
            a = xr.DataArray(a, dims=["mixture_component"])
        if not isinstance(pmix, xr.DataArray):
            pmix = xr.DataArray(pmix, dims=["mixture_component"])
        if a.sizes["mixture_component"] != pmix.sizes["mixture_component"]:
            raise ValueError(
                f"a and pmix must have the same size ({a.sizes = }, {pmix.sizes = })"
            )
        self.a: Final[xr.DataArray] = a
        self.pmix: Final[xr.DataArray] = pmix
        return

    def __repr__(self) -> str:
        return (
            f"DPsiPrior(a={[f'{x:.1e}' for x in self.a.values]},"
            f" pmix={[f'{x:.2e}' for x in self.pmix.values]})"
        )

    def discretized_logpmf(
        self, psibins: int = constants.DEFAULT_QUANTIFY_PSIBINS, PSEUDO: float = 1e-20
    ) -> xr.DataArray:
        """Get discretized logprior for deltapsi with 2 * psibins bins

        Parameters
        ----------
        psibins: int
            How many bins to discretize psi with (twice as many for deltapsi)
        PSEUDO: float
            Add a small pseudocount to the probability of each bin for logprior
        """
        endpoints = xr.DataArray(
            np.linspace(-1, 1, 1 + 2 * psibins),
            dims="pmf_bin",
        )
        cdf = xr.dot(
            xr.apply_ufunc(beta_dist.cdf, endpoints, self.a, self.a, -1, 2),
            self.pmix,
            dims=["mixture_component"],
        )
        pmf = (
            cdf.isel(pmf_bin=slice(1, None)) - cdf.isel(pmf_bin=slice(None, -1))
        ).assign_coords(
            pmf_bin_start=endpoints.isel(pmf_bin=slice(None, -1)),
            pmf_bin_end=endpoints.isel(pmf_bin=slice(1, None)),
        )
        PSEUDO = 1e-20
        return np.log(PSEUDO + pmf)

    def empirical_update(
        self,
        psi1: PsiCoverage,
        psi2: PsiCoverage,
        minreads: float = constants.DEFAULT_DPSI_PRIOR_MINREADS,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        min_lsvs: int = constants.DEFAULT_DPSI_PRIOR_MINLSV,
        n_update_a: int = constants.DEFAULT_DPSI_PRIOR_MAXITER,
        n_update_pmix: Optional[int] = None,
        legacy: bool = False,
        show_progress: bool = False,
    ) -> "DPsiPrior":
        """Use reliable binary events from psi1,2 to return updated prior

        Use high confidence empirical deltapsi from input groups of experiments
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
        min_lsvs: int
            If less than this many binary events meeting criteria, don't
            attempt an update to prior
        n_update_a: int
            Number of iterations to update `a` during M step
        n_update_pmix: Optional[int]
            Optional number of iterations to update `pmix` during M step.
            If not specified, use 1 + n_update_a.
        legacy: bool
            If True, use old implementation in v2 that does hard binning of
            observed dpsi into hard-coded bins at differences of 5% and 30%.
        show_progress: bool
            Attempt to show progress on distributed cluster for Dask
        """
        log = get_logger()
        # how many updates do we want to make?
        if n_update_pmix is None:
            n_update_pmix = 1 + n_update_a
        if max(n_update_a, n_update_pmix) < 1:
            # if we don't want to update anything, don't
            return self
        # don't bother getting empirical dpsi if futile
        if psi1.num_events < max(1, min_lsvs):
            log.info(
                "It is impossible for enough events (%d) to be identified."
                " Will not adjust prior.",
                min_lsvs,
            )
            return self

        # otherwise, get empirical dpsi to make adjustment
        dpsi = self.get_empirical_dpsi(
            psi1,
            psi2,
            minreads=minreads,
            min_experiments_f=min_experiments_f,
            show_progress=show_progress,
        )
        # do we have enough observations to do adjustment?
        if dpsi.sizes["lsv_idx"] < min_lsvs:
            log.info(
                f"Only {dpsi.sizes['lsv_idx']} reliable binary events identified"
                " to update deltapsi prior."
                f" Will not adjust prior since less than threshold of {min_lsvs}."
            )
            return self
        # otherwise
        log.info(f"Adjusting deltapsi prior using {dpsi.sizes['lsv_idx']} events")
        # update current parameters
        if legacy:
            log.info("Requested deprecated legacy approach for setting prior")
            return self.legacy_empirical_replace(dpsi)
        else:
            return self.empirical_update_EM(
                dpsi, self.a, self.pmix, n_update_a, n_update_pmix
            )

    @classmethod
    def empirical_update_EM(
        cls,
        dpsi: xr.DataArray,
        a: xr.DataArray,
        pmix: xr.DataArray,
        n_update_a: int,
        n_update_pmix: int,
    ) -> "DPsiPrior":
        for i in range(max(n_update_a, n_update_pmix)):
            # E step: which mixture component given observed dpsi
            pmix_given_dpsi = cls.infer_pmix_given_dpsi(dpsi, a, pmix)
            # M step (method of moments for beta parameters)
            if i < n_update_pmix:
                pmix = cls.fit_pmix(pmix_given_dpsi)
            if i < n_update_a:
                a = cls.fit_a(dpsi, pmix_given_dpsi)
        return DPsiPrior(a, pmix)

    @classmethod
    def legacy_empirical_replace(cls, dpsi: xr.DataArray) -> "DPsiPrior":
        """This isn't an update so much as replacement"""
        mask_slab = cast(xr.DataArray, np.abs(dpsi) > 0.30)
        mask_spike = cast(xr.DataArray, np.abs(dpsi) <= 0.05)
        mask_center = cast(xr.DataArray, ~(mask_slab | mask_spike))
        pmix_given_dpsi = xr.concat(
            [1.0 * mask_slab, 1.0 * mask_center, 1.0 * mask_spike],
            dim="mixture_component",
        )
        pmix = cls.fit_pmix(pmix_given_dpsi)
        a = cls.fit_a(dpsi, pmix_given_dpsi)
        return DPsiPrior(a, pmix)

    @staticmethod
    def get_empirical_dpsi(
        psi1: PsiCoverage,
        psi2: PsiCoverage,
        minreads: float = constants.DEFAULT_DPSI_PRIOR_MINREADS,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        show_progress: bool = False,
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
        show_progress: bool
            Attempt to show dask progress bar for distributed cluster

        Returns
        -------
        xr.DataArray
            difference in raw_psi between sum of two groups for selected events
            Dimensions: lsv_idx
        """
        reads1 = psi1.raw_total * psi1.raw_psi
        passed1 = (psi1.event_passed & (reads1 >= minreads)).sum(
            "prefix"
        ) >= min_experiments(min_experiments_f, psi1.num_prefixes)
        reads2 = psi2.raw_total * psi2.raw_psi
        passed2 = (psi2.event_passed & (reads2 >= minreads)).sum(
            "prefix"
        ) >= min_experiments(min_experiments_f, psi2.num_prefixes)
        # mask includes potentially duplicate events
        passed = passed1 & passed2 & (psi1.event_size == 2)
        # summarize over prefixes
        reads1 = reads1.sum("prefix")
        reads2 = reads2.sum("prefix")
        raw_total1 = psi1.raw_total.sum("prefix")
        raw_total2 = psi2.raw_total.sum("prefix")
        # get dpsi that passed, keeping only one per lsv
        dpsi = (
            (reads2 / raw_total2.where(raw_total2 > 0))
            - (reads1 / raw_total1.where(raw_total1 > 0))
        ).where(passed)
        if show_progress and dpsi.chunks:
            dpsi = dpsi.persist()
            progress(dpsi.data)
        dpsi = (
            dpsi.load()
            .assign_coords(lsv_idx=psi1.lsv_idx)
            .dropna("ec_idx")
            .groupby("lsv_idx")
            .first()
        )
        return dpsi

    @staticmethod
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

    @staticmethod
    def fit_pmix(pmix_given_dpsi: xr.DataArray) -> xr.DataArray:
        """Fit pmix using pmix_given_dpsi (M-step)"""
        x = pmix_given_dpsi.sum("lsv_idx")
        return x / x.sum()  # make sure result is normalized

    @staticmethod
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
        return a.where(a.notnull(), 1.0)  # if invalid prior, set to uniform
