"""
moccasin.py

Adaptation of MOCCASIN to new-majiq

Author: Joseph K Aicher
"""

import itertools
import numpy as np
import xarray as xr
import dask.array as da
import new_majiq.constants as constants

from new_majiq.Quantifier import (
    QuantifiableEvents,
    QuantifierThresholds,
)

from pathlib import Path
from typing import (
    Final,
    Union,
)


def tmp_psi_for_moccasin(
    lsv_coverage_file: Union[str, Path],
    thresholds: QuantifierThresholds = QuantifierThresholds(),
) -> xr.Dataset:
    """reparameterize bootstrap coverage to psi and total_coverage

    Batch correction is performed on scaled psi (by MLE) per bootstrap
    replicate over many samples. Except for trivial cases, this is too large to
    fit and memory. To ease chunking automatically with dask, we precompute
    these values once and keep them as temporary files that can easily be read
    in. (otherwise computing psi directly from majiq files in chunks requires
    operating over entire ec_idx dimension unless doing custom chunking that
    standard libraries have difficulty with)

    Parameters
    ----------
    lsv_coverage_file: Union[str, Path]
        Path to LSVCoverage dataset (i.e. majiq file)
    thresholds: QuantifierThresholds
        thresholds used to determine quantifiability

    Returns
    -------
    xr.Dataset
        Variables:
            psi (ec_idx, bootstrap_replicate)
            total_coverage (ec_idx, bootstrap_replicate)
        Coordinates:
            _offsets (_eidx_offset)
    """
    qe = QuantifiableEvents.from_quantifier_group([lsv_coverage_file], thresholds)
    with xr.open_zarr(
        lsv_coverage_file, group="events_coverage"
    ) as df_ec, xr.open_zarr(lsv_coverage_file, group="events") as df_e:
        offsets = df_e._offsets.values.astype(np.int64)
        bootstraps = df_ec.bootstraps.transpose("bootstrap_replicate", "ec_idx").values
    # get coverage over ec_idx for events defined by offsets
    total_coverage = np.repeat(
        np.add.reduceat(bootstraps, offsets[:-1], axis=-1), np.diff(offsets), axis=-1
    )
    with np.errstate(invalid="ignore"):  # ignore invalid when total_coverage = 0
        psi = bootstraps / total_coverage
    # mask events with no coverage or that didn't pass quantifiability thresholds
    psi = np.where(qe.event_connection_passed_mask & (total_coverage > 0), psi, np.nan)
    return xr.Dataset(
        {
            "psi": (("bootstrap_replicate", "ec_idx"), psi),
            "total_coverage": (("bootstrap_replicate", "ec_idx"), total_coverage),
        },
        {
            "_offsets": ("_eidx_offset", offsets),
        },
    )


def bootstraps_from_tmp_psi(
    psi: xr.DataArray,
    total_coverage: xr.DataArray,
    offsets: xr.DataArray,
    original_bootstraps: xr.DataArray,
) -> xr.DataArray:
    """invert (potentially corrected) psi values back to lsv_coverage bootstraps

    Parameters
    ----------
    psi: xr.DataArray
        values of psi, potentially shifted by correction
    total_coverage: xr.DataArray
        total coverage per event (over ec_idx) expected after correction
    offsets: xr.DataArray
        offsets for events
    original_bootstraps: xr.DataArray
        retain old values if unable to perform correction

    Notes
    -----
    Clips psi to non-negative values, renormalizes each event to total_coverage
    and maps back to original_bootstraps
    """
    psi_arr = psi.clip(min=0).transpose("bootstrap_replicate", "ec_idx").values
    total_coverage_arr = total_coverage.transpose(
        "bootstrap_replicate", "ec_idx"
    ).values
    offsets_arr = offsets.values
    new_total_coverage_arr = np.repeat(
        np.add.reduceat(psi_arr, offsets_arr[:-1], axis=-1),
        np.diff(offsets_arr),
        axis=-1,
    )
    # corrected values
    # TODO probably need to add with np.errstate(invalid="ignore") here
    bootstraps = xr.DataArray(
        np.where(
            new_total_coverage_arr > 0,
            # sum over event will now add up to total_coverage_arr
            psi_arr * total_coverage_arr / new_total_coverage_arr,
            np.nan,
        ),
        dims=["bootstrap_replicate", "ec_idx"],
    )
    # we don't allow one bootstrap replicate to be corrected but not others
    corrected = bootstraps.notnull().all("bootstrap_replicate")
    return (
        original_bootstraps
        # replace with corrected values if was able to do correction
        .where(~corrected, bootstraps)
        # add boolean mask indicating if correction was performed
        .assign_coords(corrected=corrected)
    )


def _silent_linalg_solve(*args, **kwargs) -> np.ndarray:
    """wraps internals of np.linalg.solve

    np.linalg.solve solves Ax = b in vectorized manner, but throws
    exception if it is able to detect any of the matrices as singular
    (doesn't always work due to floating point rounding)
    It uses np.linalg._umath_linalg.solve1 (or 2) to do the underlying
    calculation, which returns nan when it detects singular matrices.
    But it raises warnings in these cases, so we wrap it with appropriate
    context manager to silence the warning (which we expect and are
    depending on).
    """
    with np.errstate(invalid="ignore"):
        return np.linalg._umath_linalg.solve1(*args, **kwargs)


def infer_model_params(
    uncorrected: xr.DataArray, factors: xr.DataArray, complete: bool = False
) -> xr.DataArray:
    """get coefficients for OLS on uncorrected data given factors

    Parameters
    ----------
    uncorrected: xr.DataArray
        Observed values per prefix. NaN values indicate missing data.
        Dimensions: (prefix, ...)
    factors: xr.DataArray
        Confounding and non-confounding factors which are used to predict the
        values in uncorrected. Must have all prefixes found in uncorrected
        Dimensions: (prefix, factor)
    complete: bool
        Indicates if we should handle missing data in uncorrected. If the data
        are complete, can be much faster doing this

    Returns
    -------
    xr.DataArray
        Coefficients such that xr.dot(factors, coeff, dim="factor") is the OLS
        estimator of uncorrected given observed data
        Dimensions: (..., factor)
    """
    # get factors only for prefixes in uncorrected, make same dtype as uncorrected
    factors = factors.sel(prefix=uncorrected.prefix).astype(uncorrected.dtype)
    gramian: xr.DataArray  # X^T X for factors (potentially given missing data)
    projection: xr.DataArray  # X^T Y (observed uncorrected data onto factors)
    if not complete:
        missing = uncorrected.notnull()
        gramian = xr.dot(
            missing, factors.rename(factor="fsolve"), factors, dims="prefix"
        )
        gramian = gramian.where(
            # mask singular gramians as nan
            xr.apply_ufunc(
                np.linalg.matrix_rank,
                gramian,
                input_core_dims=[["fsolve", "factor"]],
                dask="parallelized",
            )
        )
        projection = xr.dot(
            missing, uncorrected, factors.rename(factor="fsolve"), dims="prefix"
        )
    else:
        gramian = xr.dot(factors.rename(factor="fsolve"), factors, dims="prefix")
        projection = xr.dot(uncorrected, factors.rename(factor="fsolve"), dims="prefix")
    # solve for OLS parameters
    return xr.apply_ufunc(
        _silent_linalg_solve,
        gramian,
        projection,
        input_core_dims=[["fsolve", "factor"], ["fsolve"]],
        output_core_dims=[["factor"]],
        dask="parallelized",
        output_dtypes=[gramian.dtype],
    )


def correct_with_model(
    uncorrected: xr.DataArray, params: xr.DataArray, factors: xr.DataArray
) -> xr.DataArray:
    """correct uncorrected given factors and params

    Parameters
    ----------
    uncorrected: xr.DataArray
        Observed values per prefix. NaN values indicate missing data.
        Dimensions: (?prefix, ...)
    params: xr.DataArray
        Model parameters from infer_model_params()
    factors: xr.DataArray
        Confounding and non-confounding factors which are used to predict the
        values in uncorrected. Must have all prefixes found in uncorrected if
        uncorrected has prefixes (if uncorrected does not have prefix, then
        factors must not either).
        Coordinate "confounding" is boolean mask over factors indicating if
        factor is confounding vs unconfounding
        Dimensions: (?prefix, factor)
        Coordinates: confounding (dim: factor)

    Returns
    -------
    xr.DataArray
        Corrected values. Residuals from predictions with all factors (vs true
        values) added to predictions with nonconfounding factors only
    """
    if "prefix" in uncorrected.dims:
        factors = factors.sel(prefix=uncorrected.prefix)
    else:
        if "prefix" in factors.dims:
            raise ValueError(
                "factors cannot have prefix dimension when uncorrected does not"
            )
    factors = factors.astype(uncorrected.dtype)  # same type as uncorrected
    predictions = xr.dot(factors, params, dims="factor")  # all factors
    residuals = uncorrected - predictions  # residuals from all factors
    unconfounded_predictions = xr.dot(
        # zeroing out confounding factors
        factors.where(~factors["confounding"], 0),  # type: ignore
        params,
        dims="factor",
    )
    return unconfounded_predictions + residuals


class ModelUnknownConfounders(object):
    """Model for how to use observations and factors to find unknown confounders"""

    def __init__(
        self,
        original_ecidx: xr.DataArray,
        model_params: xr.DataArray,
        singular_values: xr.DataArray,
        ec_vectors: xr.DataArray,
        total_variance: float,
    ):
        """Model of unknown confounders with specified parameters

        Parameters
        ----------
        original_ecidx: xr.DataArray
            indexes into ec_idx from original event connections
            (dimension: top_ec_idx)
        model_params: xr.DataArray
            coefficients for linear model to get residuals on which unknown
            confounders are computed
            (dimension: top_ec_idx, factor)
        singular_values: xr.DataArray
            singular values for unknown confounders found from training data
            (dimension: new_factor)
        ec_vectors: xr.DataArray
            right-singular vectors over top_ec_idx used to project to new
            factors (scaled appropriately by singular values)
            (dimension: new_factor, top_ec_idx)
        total_variance: float
            sum of variances scaled by the number of experiments used to train.
            This allows us to calculate the variance explained by the singular
            values
        """
        self.original_ecidx: Final[xr.DataArray] = original_ecidx
        self.model_params: Final[xr.DataArray] = model_params
        self.singular_values: Final[xr.DataArray] = singular_values
        self.ec_vectors: Final[xr.DataArray] = ec_vectors
        self.total_variance: Final[float] = total_variance
        return

    @classmethod
    def train(
        cls,
        uncorrected: xr.DataArray,
        offsets: np.ndarray,
        factors: xr.DataArray,
        max_new_factors: int,
        max_events: int = constants.DEFAULT_MOCCASIN_RUV_MAX_EVENTS,
    ) -> "ModelUnknownConfounders":
        """Learn model for unknown confounders using residuals from OLS model

        Parameters
        ----------
        uncorrected: xr.DataArray (dims: prefix, ..., ec_idx)
            data from which variations are used to identify confounding factors
        offsets: np.ndarray
            offsets into ec_idx describing different events (only use 0 or 1
            connection per event in model)
        factors: xr.DataArray (dims: prefix, factor)
            confounding/nonconfounding factors used to predict uncorrected.
            Residuals from OLS model used for SVD
        max_new_factors: int
            maximum number of new factors the model will try to produce
        max_events: int
            maximum number of events used to create the model

        Returns
        -------
        ModelUnknownConfounders
            Model picks top events that are fully quantified and have highest
            variance.
            Model uses factors to build model predicting uncorrected.
            Residuals are used as input into SVD, and appropriate parameters
            saved to be able to produce new unknown confounders for any future
            sample
        """
        factors = (
            factors.sel(prefix=uncorrected.prefix).astype(uncorrected.dtype).load()
        )
        if (factors_rank := np.linalg.matrix_rank(factors.values)) < factors.sizes[
            "factor"
        ]:
            raise ValueError(
                f"factors is not full rank ({factors_rank} < {factors.sizes['factor']}"
            )
        # get indexes for connections for top LSVs
        uncorrected_medians = cls._median_extra_dims(uncorrected)
        original_ecidx = (
            # get variance across samples of psi (median over bootstraps, etc.)
            xr.Dataset(
                {
                    "variance": uncorrected_medians.var("prefix"),
                    "complete": uncorrected_medians.count("prefix")
                    == uncorrected.sizes["prefix"],
                },
            )
            # compute into in-memory pandas dataframe
            .to_dataframe()
            # get indexes with maximum variance, only 1 per LSV
            .assign(lsv_idx=np.repeat(np.arange(len(offsets) - 1), np.diff(offsets)))
            # keep only complete indexes, get top variances, 1 per LSV
            .loc[lambda df: df["complete"]]
            .sort_values("variance", ascending=False)
            .drop_duplicates("lsv_idx", keep="first")
            .head(max_events)
            # make appropriate DataArray of indexes for this
            .pipe(
                lambda x: xr.DataArray(
                    x.index.sort_values().values, dims=["top_ec_idx"]
                )
            )
        )
        # get data with which we will train the model
        # TODO maybe want to persist this in memory if it isn't too big
        x = cls._median_extra_dims(
            uncorrected.isel(ec_idx=original_ecidx), ec_idx="top_ec_idx"
        )
        model_params = infer_model_params(x, factors, complete=True)
        x_residuals = x - xr.dot(factors, model_params, dims="factor")
        # perform sparse SVD on x_residuals
        prefix_vectors, singular_values, ec_vectors = da.linalg.svd_compressed(
            x_residuals.transpose("prefix", "top_ec_idx").data, max_new_factors
        )
        # load model parameters into memory
        local = xr.Dataset(
            {
                "model_params": model_params,
                "prefix_vectors": xr.DataArray(
                    prefix_vectors,
                    dims=["prefix", "new_factor"],
                ),
                "singular_values": xr.DataArray(singular_values, dims=["new_factor"]),
                "ec_vectors": xr.DataArray(
                    ec_vectors,
                    dims=["new_factor", "top_ec_idx"],
                ),
                "total_variance": xr.dot(x_residuals, x_residuals),
            },
            {
                "prefix": factors["prefix"],
                "new_factor": list(
                    itertools.islice(
                        # make sure there are no conflicts in names of unknown factors
                        (
                            name
                            for name in (f"_U_{n}" for n in itertools.count(1))
                            if name not in factors["factor"]
                        ),
                        max_new_factors,
                    )
                ),
                "confounding": ("new_factor", [True for _ in range(max_new_factors)]),
            },
        ).compute()
        # determine if new factors from training data (prefix_vectors) are
        # linearly independent from existing factors (we only want the linearly
        # independent ones). We can do this by looking at residuals from trying
        # to reconstruct prefix_vectors using factors
        _, old_new_residuals, _, _ = np.linalg.lstsq(
            factors.transpose("prefix", "factor").values,
            local.prefix_vectors.transpose("prefix", "new_factor").values,
            rcond=None,
        )
        # residuals > 0, accounting for potential accumulation of error
        keep_new_factors = old_new_residuals > (
            factors.sizes["prefix"] * np.finfo(old_new_residuals.dtype).eps
        )
        return ModelUnknownConfounders(
            original_ecidx=original_ecidx,
            model_params=local.model_params,
            singular_values=local.singular_values.isel(new_factor=keep_new_factors),
            ec_vectors=local.ec_vectors.isel(new_factor=keep_new_factors),
            total_variance=local.total_variance.values[()],
        )

    @property
    def explained_variance(self) -> xr.DataArray:
        """proportion of training data residuals variance explained by each new factor"""
        return (self.singular_values * self.singular_values) / self.total_variance

    @staticmethod
    def _median_extra_dims(
        x: xr.DataArray, prefix: str = "prefix", ec_idx: str = "ec_idx"
    ) -> xr.DataArray:
        """collapse any dimensions that aren't prefix or ec_idx by median"""
        extra_dims = [x for x in x.dims if x not in (prefix, ec_idx)]
        return x.median(extra_dims)

    def predict(self, uncorrected: xr.DataArray, factors: xr.DataArray) -> xr.DataArray:
        """get unknown confounders using observed data/factors

        Missing values in uncorrected have their residuals imputed to 0 (by
        definition experiments used in training must be present, but not true
        for held out data)
        """
        if "prefix" in uncorrected.dims:
            factors = factors.sel(prefix=uncorrected.prefix)
        else:
            if "prefix" in factors.dims:
                raise ValueError(
                    "factors cannot have prefix dimension when uncorrected does not"
                )
        factors = factors.astype(uncorrected.dtype)
        # get subset of uncorrected to use
        x = uncorrected.isel(ec_idx=self.original_ecidx)  # select subset
        x = self._median_extra_dims(x, ec_idx="top_ec_idx")  # summarize over extra dims
        # get residuals, imputing 0 for missing values possible in held-out data
        residuals = x - xr.dot(factors, self.model_params, dims="factor")
        residuals = residuals.where(residuals.notnull(), 0)
        # project residuals to new factors
        return xr.dot(
            residuals,
            1 / self.singular_values,  # type: ignore
            self.ec_vectors,
            dims="top_ec_idx",
        )
