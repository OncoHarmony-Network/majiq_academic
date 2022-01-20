"""
PMFSummaries.py

Helper functions for computing quantities with discretized PMFs.
Generally expect these objects to have distributions on core dimension
"pmf_bin" and coordinates "pmf_bin_start" and "pmf_bin_end" for their endpoints

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import Final, cast

import numpy as np
import xarray as xr


class PMFSummaries(object):
    """Calculate summary statistics for input probability distribution

    Parameters
    ----------
    p: xr.DataArray
        Probability distribution(s) over intervals (dimension: pmf_bin)
        Dimensions: [..., "pmf_bin"]
        Coordinates: pmf_bin_start, pmf_bin_end (..., "pmf_bin")
    """

    def __init__(self, p: xr.DataArray):
        """Calculate summary statistics for input probability distribution

        Parameters
        ----------
        p: xr.DataArray
            Probability distribution(s) over intervals (dimension: pmf_bin)
            Dimensions: [..., "pmf_bin"]
            Coordinates: pmf_bin_start, pmf_bin_end (..., "pmf_bin")
        """
        if any(x not in p.coords for x in ("pmf_bin_start", "pmf_bin_end")):
            raise ValueError(
                "Input PMF does not have required coordinates"
                " `pmf_bin_start` and `pmf_bin_end`"
            )
        if "pmf_bin" not in p.dims:
            raise ValueError("Input PMF does not have required dimension `pmf_bin`")
        Z = p.sum("pmf_bin")
        self.p: Final[xr.DataArray] = p / Z.where(Z > 0)
        return

    @property
    def bin_start(self) -> xr.DataArray:
        """start of each bin"""
        return self.p["pmf_bin_start"]

    @property
    def bin_end(self) -> xr.DataArray:
        """end of each bin"""
        return self.p["pmf_bin_end"]

    @cached_property
    def midpoints(self) -> xr.DataArray:
        """midpoints of PMF bins"""
        return 0.5 * (self.bin_start + self.bin_end)

    @staticmethod
    def interval_width(a: xr.DataArray, b: xr.DataArray) -> xr.DataArray:
        """width of interval [a, b]"""
        return (b - a).clip(min=0.0)

    @cached_property
    def bin_width(self) -> xr.DataArray:
        """extract width of each bin for PMF"""
        return self.interval_width(self.bin_start, self.bin_end)

    @staticmethod
    def f_expectation(p: xr.DataArray, f: xr.DataArray) -> xr.DataArray:
        """expectation of array f over distribution p"""
        return xr.dot(p, f, dims="pmf_bin")

    @cached_property
    def mean(self) -> xr.DataArray:
        """expectation of position in bins"""
        return self.f_expectation(self.p, self.midpoints)

    @cached_property
    def variance_per_bin(self) -> xr.DataArray:
        """variance of position in a selected bin under uniform distribution"""
        # for U ~ Uniform[a, b], Var(U) = (b - a)**2 / 12
        return cast(xr.DataArray, np.square(self.bin_width) / 12)

    @cached_property
    def variance_across_bins(self) -> xr.DataArray:
        """variance choosing between bins (treating as point mass at midpoints)"""
        residual = self.midpoints - self.mean
        return self.f_expectation(self.p, cast(xr.DataArray, np.square(residual)))

    @cached_property
    def variance(self) -> xr.DataArray:
        """variance of position (between bins, and within bins)"""
        mean_variance_per_bin = self.f_expectation(self.p, self.variance_per_bin)
        return mean_variance_per_bin + self.variance_across_bins

    @cached_property
    def standard_deviation(self) -> xr.DataArray:
        """standard deviation (sqrt of variance)"""
        return cast(xr.DataArray, np.sqrt(self.variance))

    def interval_probability(self, a: float, b: float) -> xr.DataArray:
        """compute probability of value on interval [a, b]"""
        if b <= a:
            # there is no interval, so probability should be 0
            return self.f_expectation(self.p, xr.DataArray(0))
        # compute percentage of each interval that overlaps [a, b]
        overlap_start = cast(xr.DataArray, np.maximum(a, self.bin_start))
        overlap_end = cast(xr.DataArray, np.minimum(b, self.bin_end))
        overlap_width = self.interval_width(overlap_start, overlap_end)
        overlap_pct = overlap_width / self.bin_width
        # the probability of being in [a, b] is the expectation of overlap_pct
        return self.f_expectation(self.p, overlap_pct)
