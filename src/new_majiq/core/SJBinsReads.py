"""
SJBinsReads.py

Wrap shared aspects of junction/intron raw SJ reads

Author: Joseph K Aicher
"""

from typing import TYPE_CHECKING, Final, Optional, Union

import numpy as np
import numpy.typing as npt

import new_majiq.constants as constants
from new_majiq.internals import ExperimentStrandness

if TYPE_CHECKING:
    import sparse


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
    def _sjbin_idx(self) -> npt.NDArray[np.int64]:
        """Index over bins/regions"""
        return np.arange(len(self))

    @property
    def _regions(self):
        """Underlying SJRegions the bin reads are defined over"""
        return self._sj_binsreads._regions

    @property
    def _offsets(self) -> npt.NDArray[np.uint64]:
        """raw offsets for regions (1 more than length of regions) into bins"""
        return self._sj_binsreads._offsets

    @property
    def _region_idx_start(self) -> npt.NDArray[np.uint64]:
        """offset into start of bins for each region"""
        return self._offsets[:-1]

    @property
    def _region_idx_end(self) -> npt.NDArray[np.uint64]:
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
    def bin_reads(self) -> npt.NDArray:
        """number of reads for each region/bin"""
        return self._sj_binsreads.bin_reads

    @property
    def bin_idx(self) -> npt.NDArray[np.int32]:
        """index of each bin on each region"""
        return self._sj_binsreads.bin_idx

    @property
    def bin_reads_2d(self) -> Union["sparse.GCXS", npt.NDArray]:
        """2D array representation of reads per region/bin

        Notes
        -----
        If `sparse` module is installed, will use sparse representation over
        internal arrays (no copy). Otherwise, will return dense array, which
        will increase memory utilization, potentially significantly
        """
        try:
            import sparse

            return sparse.GCXS(
                (self.bin_reads, self.bin_idx, self._offsets.view(np.int64)),
                shape=(len(self._regions), self.total_bins),
                compressed_axes=(0,),
            )
        except ModuleNotFoundError:
            # create dense array
            result = np.zeros(
                (len(self._regions), self.total_bins), dtype=self.bin_reads.dtype
            )
            # fill nonzero values
            result[
                np.repeat(np.arange(len(self._regions)), self.numbins()),
                self.bin_idx,
            ] = self.bin_reads
            return result

    def numstacks(
        self,
        region_idx: Optional[npt._ArrayLikeInt_co] = None,
        pvalue_threshold: npt._ArrayLikeFloat_co = constants.DEFAULT_COVERAGE_STACK_PVALUE,
    ) -> npt.NDArray[np.int32]:
        """Calculate number of stacks at specified regions/pvalues"""
        if region_idx is None:
            region_idx = np.arange(len(self._regions))
        return self._sj_binsreads.numstacks(region_idx, pvalue_threshold)

    def numbins(
        self,
        region_idx: Optional[npt._ArrayLikeInt_co] = None,
        minreads: npt._ArrayLikeInt_co = 0,
    ) -> npt.NDArray[np.int32]:
        """Number of bins for regions with at least specified number of reads"""
        if region_idx is None:
            region_idx = np.arange(len(self._regions))
        return self._sj_binsreads.numbins(region_idx, minreads)

    def numreads(
        self,
        region_idx: Optional[npt._ArrayLikeInt_co] = None,
        numstacks: npt._ArrayLikeInt_co = 0,
    ) -> npt.NDArray:
        """Calculate total number of reads at specified regions with stacks"""
        if region_idx is None:
            region_idx = np.arange(len(self._regions))
        return self._sj_binsreads.numreads(region_idx, numstacks)
