"""
SJBinsReads.py

Wrap shared aspects of junction/intron raw SJ reads

Author: Joseph K Aicher
"""

from typing import Final, Optional

import numpy as np
from numpy.typing import ArrayLike

import new_majiq.constants as constants
from new_majiq.internals import ExperimentStrandness


class SJBinsReads(object):
    def __init__(
        self,
        sj_binsreads,
        strandness: ExperimentStrandness,
        original_path: str,
        original_version: str,
        original_time: str,
    ):
        self._sj_binsreads: Final = sj_binsreads
        self._strandness: Final[ExperimentStrandness] = strandness
        self._original_path: Final[str] = original_path
        self._original_version: Final[str] = original_version
        self._original_time: Final[str] = original_time
        return

    @property
    def strandness(self) -> ExperimentStrandness:
        return self._strandness

    @property
    def original_path(self) -> str:
        return self._original_path

    @property
    def original_version(self) -> str:
        return self._original_version

    @property
    def original_time(self) -> str:
        return self._original_time

    def __eq__(self, other) -> bool:
        try:
            return self._sj_binsreads == self._sj_binsreads
        except AttributeError:
            return False

    def __len__(self) -> int:
        """Number of bins over all regions"""
        return len(self._sj_binsreads)

    @property
    def _sjbin_idx(self) -> np.ndarray:
        """Index over bins/regions"""
        return np.arange(len(self))

    @property
    def _regions(self):
        """Underlying SJRegions the bin reads are defined over"""
        return self._sj_binsreads._regions

    @property
    def _offsets(self) -> np.ndarray:
        """raw offsets for regions (1 more than length of regions) into bins"""
        return self._sj_binsreads._offsets

    @property
    def _region_idx_start(self) -> np.ndarray:
        """offset into start of bins for each region"""
        return self._offsets[:-1]

    @property
    def _region_idx_end(self) -> np.ndarray:
        """offset after end of bins for each region"""
        return self._offsets[1:]

    def _slice_for_region(self, region_idx: int) -> slice:
        """Get slice into bins for speciied region"""
        return slice(self._offsets[region_idx], self._offsets[1 + region_idx])

    @property
    def total_bins(self) -> int:
        """total possible number of bins (defined by positions for junctions)"""
        return self._sj_binsreads.total_bins

    @property
    def bin_reads(self) -> np.ndarray:
        """number of reads for each region/bin"""
        return self._sj_binsreads.bin_reads

    @property
    def bin_idx(self) -> np.ndarray:
        """index of each bin on each region"""
        return self._sj_binsreads.bin_idx

    def numstacks(
        self,
        region_idx: Optional[ArrayLike] = None,
        pvalue_threshold: ArrayLike = constants.DEFAULT_COVERAGE_STACK_PVALUE,
    ) -> ArrayLike:
        """Calculate number of stacks at specified regions/pvalues"""
        if region_idx is None:
            region_idx = np.arange(len(self._regions))
        return self._sj_binsreads.numstacks(region_idx, pvalue_threshold)

    def numbins(
        self,
        region_idx: Optional[ArrayLike] = None,
        minreads: ArrayLike = 0,
    ) -> ArrayLike:
        """Number of bins for regions with at least specified number of reads"""
        if region_idx is None:
            region_idx = np.arange(len(self._regions))
        return self._sj_binsreads.numbins(region_idx, minreads)

    def numreads(
        self,
        region_idx: Optional[ArrayLike] = None,
        numstacks: ArrayLike = 0,
    ) -> ArrayLike:
        """Calculate total number of reads at specified regions with stacks"""
        if region_idx is None:
            region_idx = np.arange(len(self._regions))
        return self._sj_binsreads.numreads(region_idx, numstacks)
