"""
stats.py

Majiq implementations of two-sample statistical tests omitting nan values.
Implementations of Welch's t-test, Mann-Whitney U (two-sample Wilcoxon), TNOM,
InfoScore.

This Python module wraps the direct Python C++ bindings to emulate gufunc-like
broadcasting rather than strictly requiring 2D arrays

Author: Joseph K Aicher
"""

from typing import List, Optional, Tuple

import numpy as np

import new_majiq._stats as _stats


def ttest(x: np.ndarray, labels: np.ndarray) -> np.ndarray:
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
    broadcast_shape, (x, labels) = _broadcast_to_2d(x, labels)
    pvalues2d = _stats.ttest(x, labels)
    return np.reshape(pvalues2d, broadcast_shape[:-1])


def mannwhitneyu(
    x: np.ndarray, labels: np.ndarray, sortx: Optional[np.ndarray] = None
) -> np.ndarray:
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
    broadcast_shape, (x, sortx, labels) = _broadcast_to_2d(x, sortx, labels)
    pvalues2d = _stats.mannwhitneyu(x, sortx, labels)
    return np.reshape(pvalues2d, broadcast_shape[:-1])


def infoscore(
    x: np.ndarray, labels: np.ndarray, sortx: Optional[np.ndarray] = None
) -> np.ndarray:
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
    broadcast_shape, (x, sortx, labels) = _broadcast_to_2d(x, sortx, labels)
    pvalues2d = _stats.infoscore(x, sortx, labels)
    return np.reshape(pvalues2d, broadcast_shape[:-1])


def tnom(
    x: np.ndarray, labels: np.ndarray, sortx: Optional[np.ndarray] = None
) -> np.ndarray:
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
    broadcast_shape, (x, sortx, labels) = _broadcast_to_2d(x, sortx, labels)
    pvalues2d = _stats.tnom(x, sortx, labels)
    return np.reshape(pvalues2d, broadcast_shape[:-1])


def _broadcast_to_2d(
    *x: np.typing.ArrayLike,
) -> Tuple[Tuple[int, ...], List[np.ndarray]]:
    """broadcast input arrays together, then reshape to 2D

    Returns
    -------
    (shape_original, x2d)
        shape_original: Tuple[int]
            shape after broadcasting, before being reshaped to 2D
        x2d: List[np.ndarray]
            broadcast, then reshaped 2D arrays
    """
    # handle case where no arrays are passed
    if not x:
        return (tuple(), list())
    # broadcast the arrays together
    bx: List[np.ndarray] = np.broadcast_arrays(*x)
    shape_original: Tuple[int, ...] = bx[0].shape
    dim_last = shape_original[-1] if shape_original else 1
    return (shape_original, [y.reshape(-1, dim_last) for y in bx])
