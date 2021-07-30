"""
PsiCoverage.py

PSI and total coverage (raw and bootstrapped). Converted to/from EventsCoverage.
This allows simplification of the workflow with MOCCASIN bootstrap correction
and more readily parallelized analysis of arbitrarily many files by handling
dependences between junctions.

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.experiments import bam_experiment_name
from new_majiq.Quantifier import min_experiments
from new_majiq.EventsCoverage import EventsCoverage
from typing import (
    Final,
    List,
    Union,
)
from pathlib import Path


def calculate_psi_coverage(
    events_coverage: EventsCoverage,
    minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
    minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
) -> xr.Dataset:
    """Convert EventsCoverage into PSI/total coverage

    Convert EventsCoverage into PSI/total coverage. This allows for
    junctions/introns to be processed independently.
    """
    # get offsets as int (not uint)
    offsets: np.ndarray = events_coverage.events._offsets.astype(np.int64)
    event_size: np.ndarray = np.diff(offsets)
    # get whether individual connection passes thresholds
    passed = (events_coverage.numreads >= minreads) & (
        events_coverage.numbins >= minbins
    )
    # get whether any connection in event passed, per connection
    event_passed = np.repeat(np.logical_or.reduceat(passed, offsets[:-1]), event_size)
    # get total coverage per event, per connection
    raw_total = np.repeat(
        np.add.reduceat(events_coverage.numreads, offsets[:-1]), event_size
    )
    bootstrap_total = np.repeat(
        np.add.reduceat(events_coverage.bootstraps, offsets[:-1]), event_size, axis=0
    )
    # get psi per connection
    with np.errstate(divide="ignore", invalid="ignore"):
        raw_psi = np.where(raw_total > 0, events_coverage.numreads / raw_total, 0)
        bootstrap_psi = np.where(
            bootstrap_total > 0, events_coverage.bootstraps / bootstrap_total, 0
        )
    # return dataset with matched values
    return xr.Dataset(
        data_vars=dict(
            event_passed=("ec_idx", event_passed),
            raw_total=("ec_idx", raw_total),
            raw_psi=("ec_idx", raw_psi),
            bootstrap_total=(("ec_idx", "bootstrap_replicate"), bootstrap_total),
            bootstrap_psi=(("ec_idx", "bootstrap_replicate"), bootstrap_psi),
        ),
        coords=dict(event_size=("ec_idx", np.repeat(event_size, event_size))),
        attrs=dict(
            minreads=minreads,
            minbins=minbins,
            bam_path=events_coverage.bam_path,
            bam_version=events_coverage.bam_version,
        ),
    )


def save_psi_coverage(psi_coverage: xr.Dataset, path: Union[str, Path]) -> None:
    """Save PSI coverage dataset as zarr"""
    # 8192 junctions * 30 bootstrap replicates * 1000 samples * 8 bytes per ~ 2GB
    USE_CHUNKS = {"ec_idx": 8192, "bootstrap_replicate": None}
    psi_coverage.chunk(USE_CHUNKS).to_zarr(path, mode="w")  # type: ignore
    return


def open_psi_coverage(path: Union[str, Path]) -> xr.Dataset:
    """Load a single psi coverage file"""
    return xr.open_zarr(path)


def open_multi_psi_coverage(paths: List[Path]) -> xr.Dataset:
    """Load multiple psi coverage files together at once"""
    if len(set(bam_experiment_name(x) for x in paths)) < len(paths):
        raise ValueError("paths have non-unique prefixes")
    return xr.open_mfdataset(
        paths,
        engine="zarr",
        group=None,
        concat_dim="prefix",
        preprocess=lambda x: x.expand_dims(
            prefix=[bam_experiment_name(x.encoding["source"])]
        ),
        join="override",
        compat="override",
        coords="minimal",
        data_vars="minimal",
    )


def multi_psi_to_psi_coverage(
    multi_psi_coverage: xr.Dataset,
    min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
) -> xr.Dataset:
    """Aggregate coverage/psi values over all prefixes"""
    event_passed = multi_psi_coverage["event_passed"].sum("prefix") >= min_experiments(
        min_experiments_f, multi_psi_coverage.sizes["prefix"]
    )
    raw_total = multi_psi_coverage["raw_total"].sum("prefix")
    raw_coverage = (
        multi_psi_coverage["raw_total"] * multi_psi_coverage["raw_psi"]
    ).sum("prefix")
    raw_psi = (raw_coverage / raw_total).where(raw_total > 0, 0)
    bootstrap_total = multi_psi_coverage["bootstrap_total"].sum("prefix")
    bootstrap_coverage = (
        multi_psi_coverage["bootstrap_total"] * multi_psi_coverage["bootstrap_psi"]
    ).sum("prefix")
    bootstrap_psi = (bootstrap_coverage / bootstrap_total).where(bootstrap_total > 0, 0)
    return xr.Dataset(
        data_vars=dict(
            event_passed=event_passed,
            raw_total=raw_total,
            raw_psi=raw_psi,
            bootstrap_total=bootstrap_total,
            bootstrap_psi=bootstrap_psi,
        ),
        attrs=dict(original_prefix=multi_psi_coverage["prefix"].values.tolist()),
    )
