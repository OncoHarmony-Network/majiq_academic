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
import new_majiq.beta_mixture as bm

from new_majiq.experiments import bam_experiment_name
from new_majiq.Quantifier import min_experiments
from new_majiq.EventsCoverage import EventsCoverage

from functools import cached_property
from typing import (
    Final,
    List,
    Sequence,
    Tuple,
    Union,
)
from pathlib import Path


def _silent_quantile(*args, **kwargs):
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return bm.quantile(*args, **kwargs)


def _silent_cdf(*args, **kwargs):
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return bm.cdf(*args, **kwargs)


class AbstractPsiCoverage(object):
    def __init__(self, df: xr.Dataset):
        """Compute posterior quantities using information about event total/psi

        Parameters
        ----------
        df: xr.Dataset
            Data variables:
                event_passed[ec_idx]
                raw_total[ec_idx]
                raw_psi[ec_idx]
                bootstrap_total[ec_idx, bootstrap_replicate]
                bootstrap_psi[ec_idx, bootstrap_replicate]
            Coordinates:
                lsv_offsets[offset_idx]
            Derived (from _offsets):
                event_size[ec_idx]
                lsv_idx[ec_idx]
        """
        offsets = df["lsv_offsets"].load().values
        if offsets[0] != 0:
            raise ValueError("offsets[0] must be zero")
        if offsets[-1] != df.sizes["ec_idx"]:
            raise ValueError("offsets[-1] must equal number of event connections")
        event_size = np.diff(offsets)
        self.df: Final[xr.Dataset] = df.assign_coords(
            event_size=("ec_idx", np.repeat(event_size, event_size)),
            lsv_idx=("ec_idx", np.repeat(np.arange(len(event_size)), event_size)),
        )
        return

    @property
    def event_passed(self) -> xr.DataArray:
        """True if event passed per event connection"""
        return self.df["event_passed"]

    @property
    def event_size(self) -> xr.DataArray:
        """Size of event each event connection belongs to"""
        return self.df["event_size"]

    @property
    def lsv_offsets(self) -> xr.DataArray:
        """Offsets into ec_idx for LSVs (start/end)"""
        return self.df["lsv_offsets"]

    @property
    def lsv_idx(self) -> xr.DataArray:
        """Index identifier of the LSV to which each event connection belongs"""
        return self.df["lsv_idx"]

    @property
    def raw_total(self) -> xr.DataArray:
        """Raw total reads across event per event connection"""
        return self.df["raw_total"]

    @property
    def raw_psi(self) -> xr.DataArray:
        """Percentage of self.raw_total attributable to each event connection"""
        return self.df["raw_psi"]

    @property
    def bootstrap_total(self) -> xr.DataArray:
        """Bootstrapped total reads across event per event connection"""
        return self.df["bootstrap_total"]

    @property
    def bootstrap_psi(self) -> xr.DataArray:
        """Percentage of self.bootstrap_total attributable to each event connection"""
        return self.df["bootstrap_psi"]

    @cached_property
    def alpha_prior(self) -> xr.DataArray:
        return np.reciprocal(self.event_size.astype(self.raw_psi.dtype))

    @cached_property
    def beta_prior(self) -> xr.DataArray:
        return np.subtract(1, self.alpha_prior)

    @cached_property
    def raw_alpha(self) -> xr.DataArray:
        return self.alpha_prior + self.raw_psi * self.raw_total

    @cached_property
    def bootstrap_alpha(self) -> xr.DataArray:
        return self.alpha_prior + self.bootstrap_psi * self.bootstrap_total

    @cached_property
    def raw_beta(self) -> xr.DataArray:
        return np.add(1, self.raw_total - self.raw_alpha)

    @cached_property
    def bootstrap_beta(self) -> xr.DataArray:
        return np.add(1, self.bootstrap_total - self.bootstrap_alpha)

    @cached_property
    def raw_posterior_mean(self) -> xr.DataArray:
        return self.raw_alpha / np.add(1, self.raw_total)

    @cached_property
    def raw_posterior_variance(self) -> xr.DataArray:
        mean = self.raw_posterior_mean
        return mean * np.subtract(1, mean) / np.add(1, self.raw_total)

    @cached_property
    def _bootstrap_moments(self) -> Tuple[xr.DataArray, xr.DataArray]:
        # get mean, variance per bootstrap replicate
        mean = self.bootstrap_alpha / np.add(1, self.bootstrap_total)
        variance = mean * np.subtract(1, mean) / np.add(1, self.bootstrap_total)
        # summarize over bootstrap replicates
        agg_mean = mean.mean("bootstrap_replicate")
        agg_variance = (
            # Law of total variance: V(X) = V(E[X|B]) + E[V(X|B)]
            variance.mean("bootstrap_replicate")
            + mean.var("bootstrap_replicate", ddof=0)
        )
        return (agg_mean, agg_variance)

    @property
    def bootstrap_posterior_mean(self) -> xr.DataArray:
        return self._bootstrap_moments[0]

    @property
    def bootstrap_posterior_variance(self) -> xr.DataArray:
        return self._bootstrap_moments[1]

    def bootstrap_quantile(
        self,
        quantiles: Sequence[float] = [0.1, 0.9],
    ) -> xr.DataArray:
        quantiles_arr = xr.DataArray(quantiles, dims="quantiles")
        return xr.apply_ufunc(
            _silent_quantile,
            quantiles_arr,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[[], ["bootstrap_replicate"], ["bootstrap_replicate"]],
            dask="allowed",
        )

    def bootstrap_discretized_pmf(self, nbins: int = 40):
        endpoints = xr.DataArray(
            np.linspace(0, 1, 1 + nbins, dtype=self.bootstrap_psi.dtype), dims="pmf_bin"
        )
        cdf = xr.apply_ufunc(
            _silent_cdf,
            endpoints,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[[], ["bootstrap_replicate"], ["bootstrap_replicate"]],
            dask="allowed",
        )
        pmf = cdf.isel(pmf_bin=slice(1, None)) - cdf.isel(pmf_bin=slice(None, -1))
        return pmf.assign_coords(
            pmf_bin_start=("pmf_bin", endpoints[:-1]),
            pmf_bin_end=("pmf_bin", endpoints[1:]),
        )


class PsiCoverage(AbstractPsiCoverage):
    """PSI coverage for a single experiment at a time"""

    @classmethod
    def from_events_coverage(
        cls,
        events_coverage: EventsCoverage,
        minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
        minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
    ) -> "PsiCoverage":
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
        event_passed = np.repeat(
            np.logical_or.reduceat(passed, offsets[:-1]), event_size
        )
        # get total coverage per event, per connection
        raw_total = np.repeat(
            np.add.reduceat(events_coverage.numreads, offsets[:-1]), event_size
        )
        bootstrap_total = np.repeat(
            np.add.reduceat(events_coverage.bootstraps, offsets[:-1]),
            event_size,
            axis=0,
        )
        # get psi per connection
        with np.errstate(divide="ignore", invalid="ignore"):
            raw_psi = np.where(raw_total > 0, events_coverage.numreads / raw_total, 0)
            bootstrap_psi = np.where(
                bootstrap_total > 0, events_coverage.bootstraps / bootstrap_total, 0
            )
        # return dataset with matched values
        return cls(
            xr.Dataset(
                data_vars=dict(
                    event_passed=("ec_idx", event_passed),
                    raw_total=("ec_idx", raw_total),
                    raw_psi=("ec_idx", raw_psi),
                    bootstrap_total=(
                        ("ec_idx", "bootstrap_replicate"),
                        bootstrap_total,
                    ),
                    bootstrap_psi=(("ec_idx", "bootstrap_replicate"), bootstrap_psi),
                ),
                coords=dict(
                    lsv_offsets=("offset_idx", offsets),
                ),
                attrs=dict(
                    minreads=minreads,
                    minbins=minbins,
                    bam_path=events_coverage.bam_path,
                    bam_version=events_coverage.bam_version,
                ),
            )
        )

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "PsiCoverage":
        """Load a single psi coverage file"""
        # TODO verify correctness of file?
        return cls(xr.open_zarr(path))

    def to_zarr(self, path: Union[str, Path], ec_chunksize: int = 8192) -> None:
        """Save PSI coverage dataset as zarr

        Notes
        -----
        We will typically load chunks across all experiments at once for all
        bootstrap replicates. So file size will be dependent on:
        + ec_chunksize (how many event connections per chunk)
        + num_bootstraps (how many bootstraps are simultaneously there)
        + num_samples (how many prefixes are being simultaneously processed)
        + bytes per record (precision, how many arrays) at once.

        At 8 bytes per record, default num bootstraps, and 1000 samples, this is
        around 2GB with default chunksize of 8192.

        If we know we are processing at most a few samples at a time, this
        should be made higher.
        """
        USE_CHUNKS = dict(
            ec_idx=ec_chunksize, bootstrap_replicate=None, offset_idx=None
        )
        (
            self.df.drop_vars(["event_size", "lsv_idx"])
            .chunk(USE_CHUNKS)  # type: ignore
            .to_zarr(path, mode="w")
        )
        return


class MultiPsiCoverage(AbstractPsiCoverage):
    @classmethod
    def from_mf_zarr(cls, paths: List[Path]) -> "MultiPsiCoverage":
        """Load multiple psi coverage files together at once"""
        if len(set(bam_experiment_name(x) for x in paths)) < len(paths):
            raise ValueError("paths have non-unique prefixes")
        return cls(
            xr.open_mfdataset(
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
        )

    def sum(
        self, min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS
    ) -> PsiCoverage:
        """Aggregate coverage/psi values over all prefixes"""
        event_passed = self.event_passed.sum("prefix") >= min_experiments(
            min_experiments_f, self.df.sizes["prefix"]
        )
        raw_total = self.raw_total.sum("prefix")
        raw_coverage = (self.raw_total * self.raw_psi).sum("prefix")
        raw_psi = (raw_coverage / raw_total).where(raw_total > 0, 0)
        bootstrap_total = self.bootstrap_total.sum("prefix")
        bootstrap_coverage = (self.bootstrap_total * self.bootstrap_psi).sum("prefix")
        bootstrap_psi = (bootstrap_coverage / bootstrap_total).where(
            bootstrap_total > 0, 0
        )
        return PsiCoverage(
            xr.Dataset(
                data_vars=dict(
                    event_passed=event_passed,
                    raw_total=raw_total,
                    raw_psi=raw_psi,
                    bootstrap_total=bootstrap_total,
                    bootstrap_psi=bootstrap_psi,
                ),
                coords=dict(
                    lsv_offsets=("offset_idx", self.lsv_offsets),
                ),
                attrs=dict(original_prefix=self.df["prefix"].values.tolist()),
            )
        )
