"""
SpliceGraphMask.py

Boolean mask over gene introns and junctions

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Final, Optional, Union

import numpy as np
import numpy.typing as npt
import xarray as xr

import new_majiq.constants as constants
from new_majiq.internals import SpliceGraphMask as _SpliceGraphMask

from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions


class SpliceGraphMask(object):
    """Boolean mask over splicegraph introns and junctions

    Parameters
    ----------
    sg_mask: _SpliceGraphMask
        Internals representation of SpliceGraphMask

    See Also
    --------
    SpliceGraphMask.from_zarr : Load :class:`SpliceGraphMask` from Zarr file
    SpliceGraphMask.from_arrays : Load :class:`SpliceGraphMask` from arrays
    """

    def __init__(self, sg_mask: _SpliceGraphMask):
        self._sg_mask: Final[_SpliceGraphMask] = sg_mask
        return

    @property
    def introns(self) -> GeneIntrons:
        """:class:`GeneIntrons` on which the mask is defined"""
        return GeneIntrons(self._sg_mask._introns)

    @property
    def junctions(self) -> GeneJunctions:
        """:class:`GeneJunctions` on which the mask is defined"""
        return GeneJunctions(self._sg_mask._junctions)

    @property
    def introns_mask(self) -> npt.NDArray[np.bool_]:
        """boolean mask over `self.introns`"""
        return self._sg_mask.introns_values

    @property
    def junctions_mask(self) -> npt.NDArray[np.bool_]:
        """boolean mask over `self.junctions`"""
        return self._sg_mask.junctions_values

    @property
    def df(self) -> xr.Dataset:
        return xr.Dataset(
            dict(
                introns_mask=("gi_idx", self.introns_mask),
                junctions_mask=("gj_idx", self.junctions_mask),
            )
        )

    @property
    def df_introns(self) -> xr.DataArray:
        return self.introns.df.assign(introns_mask=("gi_idx", self.introns_mask))

    @property
    def df_junctions(self) -> xr.DataArray:
        return self.junctions.df.assign(junctions_mask=("gj_idx", self.junctions_mask))

    def to_zarr(self, path: Union[str, Path], mode: str = "w") -> None:
        """Save :class:`SpliceGraphMask` to specified zarr file"""
        save_df = self.df.assign_attrs(
            introns_checksum=self.introns.checksum_nodata(),
            junctions_checksum=self.junctions.checksum_nodata(),
        )
        save_df.chunk(save_df.sizes).to_zarr(path, mode=mode, group=constants.NC_SGMASK)
        return

    @classmethod
    def from_zarr(
        cls, path: Union[str, Path], introns: GeneIntrons, junctions: GeneJunctions
    ) -> "SpliceGraphMask":
        """Load :class:`SpliceGraphMask` from Zarr file to match `introns` and `junctions`"""
        with xr.open_zarr(path, mode="r", group=constants.NC_SGMASK) as ds:
            if ds.introns_checksum != introns.checksum_nodata():
                raise ValueError("Introns checksum does not match")
            if ds.junctions_checksum != junctions.checksum_nodata():
                raise ValueError("Junctions checksum does not match")
            return SpliceGraphMask.from_arrays(
                introns,
                junctions,
                ds["introns_mask"].values,
                ds["junctions_mask"].values,
            )

    @classmethod
    def from_arrays(
        cls,
        introns: GeneIntrons,
        junctions: GeneJunctions,
        introns_mask: Optional[npt._ArrayLikeBool_co] = None,
        junctions_mask: Optional[npt._ArrayLikeBool_co] = None,
    ) -> "SpliceGraphMask":
        """Create :class:`SpliceGraphMask` with specified arrays and regions

        Parameters
        ----------
        introns: GeneIntrons
        junctions: GeneJunctions
        introns_mask, junctions_mask: array_like, bool, optional
            boolean mask arrays with lengths equal to `introns` and
            `junctions`.
            If not specified, use connection flags (passed reliability filters,
            was not simplified).
        """
        if introns_mask is None:
            introns_mask = introns.passed_build & (~introns.simplified)
        if junctions_mask is None:
            junctions_mask = junctions.passed_build & (~junctions.simplified)
        return SpliceGraphMask(
            _SpliceGraphMask(
                introns._gene_introns,
                junctions._gene_junctions,
                introns_mask,
                junctions_mask,
            )
        )
