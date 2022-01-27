"""
Regions.py

Parent class that wraps new_majiq internals for regions. This is a parent class

Author: Joseph K Aicher
"""

from typing import Final, Optional

import numpy as np
import numpy.typing as npt


class Regions(object):
    def __init__(self, regions):
        self._regions: Final = regions
        return

    def __eq__(self, other):
        try:
            return self._regions == other._regions
        except AttributeError:
            return False

    def __len__(self) -> int:
        """Number of regions"""
        return len(self._regions)

    def overlaps(
        self, other, region_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.bool_]:
        """Get mask over region_idx indicating if they overlap regions in other

        Parameters
        ----------
        other:
            Regions of same type as self
        region_idx: Optional[array_like[int]]
            Indexes for regions in self. If None, compute overlaps for all
            regions in self

        Returns
        -------
        array[bool]
            array matching region_idx. An element is true when there exists a
            feature in other that overlaps self[region_idx].

        Notes
        -----
        For junctions, this just searches for a matching feature.
        Requires that the regions share the same parent.
        """
        if region_idx is None:
            region_idx = np.arange(len(self))
        return self._regions.overlaps(region_idx, other._regions)

    @property
    def _region_idx(self) -> npt.NDArray[np.int64]:
        """Index over regions"""
        return np.arange(len(self))

    @property
    def _parents(self):
        """internals class for parents on which regions defined"""
        return self._regions._parents

    @property
    def _parent_idx_start(self) -> npt.NDArray[np.uint64]:
        """First index in regions for each parent"""
        return self._regions._parent_idx_start

    @property
    def _parent_idx_end(self) -> npt.NDArray[np.uint64]:
        """One after last index in regions for each parent"""
        return self._regions._parent_idx_end

    def _slice_for_parent(self, parent_idx: int) -> slice:
        """Get slice into regions for specified parent"""
        return slice(
            self._parent_idx_start[parent_idx],
            self._parent_idx_end[parent_idx],
        )

    @property
    def start(self) -> npt.NDArray[np.int64]:
        """Start coordinate of each region"""
        return self._regions.start

    @property
    def end(self) -> npt.NDArray[np.int64]:
        """End coordinate of each region"""
        return self._regions.end
