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
