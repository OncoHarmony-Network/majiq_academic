"""
stats.py

Majiq implementations of two-sample statistical tests omitting nan values.
Implementations of Welch's t-test, Mann-Whitney U (two-sample Wilcoxon), TNOM,
InfoScore.

This Python module wraps the direct Python C++ bindings to emulate gufunc-like
broadcasting rather than strictly requiring 2D arrays

Author: Joseph K Aicher
"""

from typing import Optional

import numpy as np
import numpy.typing as npt

import new_majiq._stats as _stats


def histogram(
    x: npt._ArrayLikeFloat_co,
    x_min: npt._ArrayLikeFloat_co,
    x_max: npt._ArrayLikeFloat_co,
    dummy_bins: Optional[npt._ArrayLikeFloat_co] = None,
    nbins: Optional[int] = None,
    closed_min: bool = True,
    closed_max: bool = False,
) -> npt.NDArray[np.int64]:
    """Compute histogram over uniform bins ignoring nans.

    gufunc signature: (n),(),(),(b)->(b)

    Compute histogram over observations in `x` core dimension on the range
    `x_min` to `x_max`. Whether `x_min` or `x_max` is included in the range is
    determined by `closed_min` and `closed_max`.
    Number of bins is determined by either `dummy_bins` (size from core
    dimension), or `nbins`, which creates matching `dummy_bins` for
    convenience.
    NaN values and values outside of the range are ignored.

    Notes
    -----
    Internal implementation uses bins that are closed on the minimum value and
    open on the maximum (e.g., a <= x < b).
    Changing `closed_min` and `closed_max` from their defaults is done by
    incrementing `x_min`, `x_max` to the next representable floating point
    number to ensure they are excluded or included in the range used
    internally.

    Parameters
    ----------
    x: array[float]
        Values to be summarized as histogram. NaN and out-of-range values are
        ignored.
    x_min, x_max: float
        Minimum and maximum values for the range of observations to count
    dummy_bins: array[float] OR nbins: int
        Specify the number of bins sized uniformly on the range.
        `dummy_bins` specifies the number of bins in a gufunc style using the
        size of its core dimension (its values are otherwise ignored).
        `nbins` is available for convenience and sets `dummy_bins` to an
        appropriate value.
    closed_min, closed_max: bool
        Indicate whether `x_min`, `x_max` should be included in the range

    Returns
    -------
    array[int]
        Counts of input values on uniform bins between `x_min` and `x_max`,
        where inclusion of endpoints is determined by `closed_min` and
        `closed_max`.
    """
    # get type of x, making it an array if not available directly
    # doing this rather than np.array(x) in case x is a dask array so we don't
    # load it into memory
    try:
        x_dtype = x.dtype  # type: ignore[union-attr]
    except AttributeError:
        x = np.array(x)
        x_dtype = x.dtype
    # set dummy_bins appropriately
    try:
        dummy_nbins = None if dummy_bins is None else dummy_bins.shape[-1]  # type: ignore[union-attr]
    except AttributeError:
        dummy_nbins = np.array(dummy_bins).shape[-1]
    if nbins is not None:
        if dummy_nbins is not None and dummy_nbins != nbins:
            raise ValueError(
                "dummy_bins and nbins were both specified but disagree on"
                f" number of bins ({dummy_nbins = }, {nbins = })"
            )
        dummy_nbins = nbins
    if dummy_nbins is None:
        raise ValueError("At least one of dummy_bins or nbins must be specified")
    dummy_bins = np.empty(dummy_nbins, dtype=np.float64)
    # set x_min, x_max appropriately
    x_min = np.array(x_min, dtype=x_dtype)
    x_max = np.array(x_max, dtype=x_dtype)
    if not closed_min:
        x_min = np.nextafter(x_min, 1 + x_min, dtype=x_dtype)
    if closed_max:
        x_max = np.nextafter(x_max, 1 + x_max, dtype=x_dtype)
    # compute using internal implementation
    return _stats.histogram(x, x_min, x_max, dummy_bins)


def ttest(
    x: npt._ArrayLikeFloat_co, labels: npt._ArrayLikeBool_co
) -> npt.NDArray[np.floating]:
    """Compute p-values for Welch's t-test on input data

    gufunc signature: (n),(n)->()

    Compute p-values for Welch's t-test on input data, using a two-sided
    alternative hypothesis and omitting nan values.

    Parameters
    ----------
    x: array[float]
        test over observations in last axis
    labels: array[bool]
        test over labels in last axis

    Returns
    -------
    array[float]
        broadcast p-values for observations/labels. Invalid tests are nan
    """
    return _stats.ttest(x, labels)


def mannwhitneyu(
    x: npt._ArrayLikeFloat_co,
    labels: npt._ArrayLikeBool_co,
    sortx: Optional[npt._ArrayLikeInt_co] = None,
) -> npt.NDArray[np.floating]:
    """Compute p-values for Mann-Whitney U test on input data

    gufunc signature: (n),(n)(,(n))?->()

    Compute p-values for Mann-Whitney U test on input data, using a two-sided
    alternative hypothesis and omitting nan values.

    Parameters
    ----------
    x: array[float]
        test over observations in last axis
    labels: array[bool]
        test over labels in last axis
    sortx: Optional[array[int]]
        If specified, previously computed values of np.argsort(x, axis=-1)
        that will not be checked

    Returns
    -------
    array[float]
        broadcast p-values for observations/labels. Invalid tests are nan
    """
    if sortx is None:
        sortx = np.argsort(x, axis=-1)
    return _stats.mannwhitneyu(x, sortx, labels)


def infoscore(
    x: npt._ArrayLikeFloat_co,
    labels: npt._ArrayLikeBool_co,
    sortx: Optional[npt._ArrayLikeInt_co] = None,
) -> npt.NDArray[np.floating]:
    """Compute p-values for InfoScore test on input data

    gufunc signature: (n),(n)(,(n))?->()

    Compute p-values for InfoScore test on input data omitting nan values.

    Parameters
    ----------
    x: array[float]
        test over observations in last axis
    labels: array[bool]
        test over labels in last axis
    sortx: Optional[array[int]]
        If specified, previously computed values of np.argsort(x, axis=-1)
        that will not be checked

    Returns
    -------
    array[float]
        broadcast p-values for observations/labels. Invalid tests are nan
    """
    if sortx is None:
        sortx = np.argsort(x, axis=-1)
    return _stats.infoscore(x, sortx, labels)


def tnom(
    x: npt._ArrayLikeFloat_co,
    labels: npt._ArrayLikeBool_co,
    sortx: Optional[npt._ArrayLikeInt_co] = None,
) -> npt.NDArray[np.floating]:
    """Compute p-values for TNOM test on input data

    gufunc signature: (n),(n)(,(n))?->()

    Compute p-values for TNOM test on input data omitting nan values.

    Parameters
    ----------
    x: array[float]
        test over observations in last axis
    labels: array[bool]
        test over labels in last axis
    sortx: Optional[array[int]]
        If specified, previously computed values of np.argsort(x, axis=-1)
        that will not be checked

    Returns
    -------
    array[float]
        broadcast p-values for observations/labels. Invalid tests are nan
    """
    if sortx is None:
        sortx = np.argsort(x, axis=-1)
    return _stats.tnom(x, sortx, labels)
