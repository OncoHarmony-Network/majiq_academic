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

from new_majiq.logger import get_logger
from pathlib import Path
from typing import (
    Final,
    Hashable,
    List,
    Union,
)


# TODO we want to have new PSICoverage from old PSICoverage. Will be removed,
# but note how to renormalize PSI by clipping/dividing by new total coverage
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
    total_coverage_fn = _get_total_coverage_function(offsets.values)
    psi = psi.clip(min=0).chunk({"ec_idx": None})
    new_total_coverage = xr.apply_ufunc(
        total_coverage_fn,
        psi,
        input_core_dims=[["ec_idx"]],
        output_core_dims=[["ec_idx"]],
        output_dtypes=(psi.dtype,),
        dask="parallelized",
    )
    bootstraps = psi * total_coverage / new_total_coverage.where(new_total_coverage > 0)
    # we don't allow one bootstrap replicate to be corrected but not others
    corrected = bootstraps.notnull().all("bootstrap_replicate")
    return (
        original_bootstraps
        # replace with corrected values if was able to do correction
        .where(~corrected, bootstraps)
        # add boolean mask indicating if correction was performed
        .assign_coords(corrected=corrected)
    )


def infer_model_params(
    uncorrected: xr.DataArray,
    passed: xr.DataArray,
    factors: xr.DataArray,
    complete: bool = False,
) -> xr.DataArray:
    """get coefficients for OLS on uncorrected data given factors

    Parameters
    ----------
    uncorrected: xr.DataArray
        Observed values per prefix. NaN values are propagated regardless of
        passed vs not.
        Dimensions: (prefix, ...)
    passed: xr.DataArray
        Boolean indicator if observation passed and used for modeling or
        ignored in modeling.
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
        gramian = xr.dot(
            passed, factors.rename(factor="fsolve"), factors, dims="prefix"
        )
        gramian = gramian.where(
            # mask singular gramians as nan
            xr.apply_ufunc(
                np.linalg.matrix_rank,
                gramian,
                input_core_dims=[["fsolve", "factor"]],
                dask="parallelized",
            )
            == gramian.sizes["factor"]
        )
        projection = xr.dot(
            passed, uncorrected, factors.rename(factor="fsolve"), dims="prefix"
        )
    else:
        gramian = xr.dot(factors.rename(factor="fsolve"), factors, dims="prefix")
        projection = xr.dot(factors.rename(factor="fsolve"), uncorrected, dims="prefix")
    # we can stack extra core dimensions for efficient linear solve (solve2)
    extra_core_dims = list(set(projection.dims) - (set(gramian.dims) | set(passed.dims)))
    # solve for OLS parameters
    params: xr.DataArray
    if extra_core_dims:
        params = xr.apply_ufunc(
            _silent_linalg_solve2,
            gramian,
            projection.stack({"_extra_core_dims": extra_core_dims}).chunk(
                {"_extra_core_dims": None}
            ),
            input_core_dims=[["fsolve", "factor"], ["fsolve", "_extra_core_dims"]],
            output_core_dims=[["factor", "_extra_core_dims"]],
            dask="parallelized",
            output_dtypes=[gramian.dtype],
        )
        params = params.unstack("_extra_core_dims")
        # unstack adds coordinates to dimensions if they weren't unlabeled, so
        # we have to reverse that
        params = params.drop_vars(
            [x for x in params.coords if x not in uncorrected.coords]
        )
    else:
        params = xr.apply_ufunc(
            _silent_linalg_solve1,
            gramian,
            projection,
            input_core_dims=[["fsolve", "factor"], ["fsolve"]],
            output_core_dims=[["factor"]],
            dask="parallelized",
            output_dtypes=[gramian.dtype],
        )
    return params


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

    def to_zarr(self, output: Union[str, Path]) -> None:
        """Save model parameters fo file"""
        xr.Dataset(
            {
                "original_ecidx": self.original_ecidx,
                "model_params": self.model_params,
                "singular_values": self.singular_values,
                "ec_vectors": self.ec_vectors,
            },
            {},
            {
                "total_variance": self.total_variance,
            },
        ).to_zarr(output, mode="w")
        return

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "ModelUnknownConfounders":
        """Load model parameters from file"""
        with xr.open_zarr(path) as df:
            df.load()
            return ModelUnknownConfounders(
                original_ecidx=df.original_ecidx,
                model_params=df.model_params,
                singular_values=df.singular_values,
                ec_vectors=df.ec_vectors,
                total_variance=df.attrs["total_variance"],
            )

    @classmethod
    def train(
        cls,
        uncorrected: xr.DataArray,
        passed: xr.DataArray,
        offsets: np.ndarray,
        factors: xr.DataArray,
        max_new_factors: int = constants.DEFAULT_MOCCASIN_RUV_MAX_FACTORS,
        max_events: int = constants.DEFAULT_MOCCASIN_RUV_MAX_EVENTS,
    ) -> "ModelUnknownConfounders":
        """Learn model for unknown confounders using residuals from OLS model

        Parameters
        ----------
        uncorrected: xr.DataArray (dims: prefix, ..., ec_idx)
            data from which variations are used to identify confounding factors
        passed: xr.DataArray (dims: prefix, ..., ec_idx)
            indicator as to whether values in uncorrected passed and should be
            included in modeling
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
        uncorrected_medians = cls._median_extra_dims(uncorrected.where(passed))
        log = get_logger()
        log.info(
            "We expect likely NaN slice, invalid value in true_divide warnings"
            " in this next step (https://github.com/dask/dask/issues/3245)."
            " They can be safely ignored"
        )
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
        log.info("We do not expect warnings after this point")
        # get data with which we will train the model
        # TODO maybe want to persist this in memory if it isn't too big
        x = cls._median_extra_dims(
            uncorrected.isel(ec_idx=original_ecidx), ec_idx="top_ec_idx"
        )
        model_params = infer_model_params(x, passed, factors, complete=True)
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
        # only keep as many new factors that keep model matrix having full
        # column rank
        combined_factors = np.empty(
            (factors.sizes["prefix"], factors.sizes["factor"] + max_new_factors),
            dtype=factors.dtype,
        )
        combined_factors[:, : factors.sizes["factor"]] = factors.transpose(
            "prefix", "factor"
        ).values
        combined_factors[:, factors.sizes["factor"] :] = local.prefix_vectors.transpose(
            "prefix", "new_factor"
        ).values
        keep_n: int  # we know we can keep at least this many new factors (starting at 0)
        for keep_n in range(max_new_factors):
            # try adding 1 more after the keep_n we know works, still full rank?
            total = factors.sizes["factor"] + keep_n + 1
            if np.linalg.matrix_rank(combined_factors[:, :total]) < total:
                # combined factors now singular so keep_n is maximum new factors
                break
        else:
            # all of combined_factors was still full rank
            keep_n = max_new_factors

        return ModelUnknownConfounders(
            original_ecidx=original_ecidx,
            model_params=local.model_params,
            singular_values=local.singular_values.isel(new_factor=slice(keep_n)),
            ec_vectors=local.ec_vectors.isel(new_factor=slice(keep_n)),
            total_variance=local.total_variance.values[()],
        )

    @property
    def num_factors(self) -> int:
        """Number of new unknown factors from the model"""
        return self.singular_values.sizes["new_factor"]

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

    def predict(
        self, uncorrected: xr.DataArray, passed: xr.DataArray, factors: xr.DataArray
    ) -> xr.DataArray:
        """get unknown confounders using observed data/factors

        Unpassed values have their residuals imputed to 0 (by
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
        x = uncorrected.where(passed)  # mask unpassed values
        x = x.isel(ec_idx=self.original_ecidx)  # select subset
        x = self._median_extra_dims(x, ec_idx="top_ec_idx")  # summarize over extra dims
        # get residuals, imputing 0 for missing values possible in held-out data
        residuals = x - xr.dot(factors, self.model_params, dims="factor")
        residuals = residuals.where(residuals.notnull(), 0)
        # project residuals to new factors
        return xr.dot(
            residuals,
            np.reciprocal(self.singular_values),
            self.ec_vectors,
            dims="top_ec_idx",
        )


def _silent_linalg_solve1(*args, **kwargs) -> np.ndarray:
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


def _silent_linalg_solve2(*args, **kwargs) -> np.ndarray:
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
        return np.linalg._umath_linalg.solve(*args, **kwargs)


def _get_total_coverage_function(offsets: np.ndarray):
    """get function that summarizes totals defined by offsets"""
    offsets = offsets.astype(np.int64)
    event_size = np.diff(offsets)

    def result(x: np.ndarray) -> np.ndarray:
        """get total coverage within each event defined by offsets"""
        if x.shape[-1] != offsets[-1]:
            raise ValueError(
                "input array core dim size does not match offsets"
                f" ({x.shape = }, {offsets[-1] = })"
            )
        return np.repeat(np.add.reduceat(x, offsets[:-1], axis=-1), event_size, axis=-1)

    return result
