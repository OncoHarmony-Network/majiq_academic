"""
_workarounds.py

Helper functions to workaround known bugs while waiting for bugfix

Author: Joseph K Aicher
"""

import xarray as xr


def _load_zerodim_variables(x: xr.Dataset) -> xr.Dataset:
    """Load variables in x that have zero-length dimensions

    This addresses exceptions when saving dask arrays with zero-length
    dimension in xarray/zarr (xarray issue #5741 / PR #5742)
    """
    for v in x.variables.values():
        if v.size == 0:
            v.load()
    return x
