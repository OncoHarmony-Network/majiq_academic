"""
PsiCoverage.py

PSI and total coverage (raw and bootstrapped). Converted to/from EventsCoverage.
This allows simplification of the workflow with MOCCASIN bootstrap correction
and more readily parallelized analysis of arbitrarily many files by handling
dependences between junctions.

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import Dict, Final, List, Optional, Sequence, Tuple, Union, cast

import numpy as np
import xarray as xr

import new_majiq._offsets as _offsets
import new_majiq.beta_mixture as bm
import new_majiq.constants as constants
from new_majiq.Events import Events, _Events
from new_majiq.EventsCoverage import EventsCoverage
from new_majiq.experiments import bam_experiment_name
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions


def min_experiments(min_experiments_f: float, num_experiments: int) -> float:
    if min_experiments_f < 1:
        min_experiments_f *= num_experiments
    return max(1, min(min_experiments_f, num_experiments))


class PsiCoverage(object):
    """PSI and total coverage (raw and bootstrapped) for arbitrary number of samples"""

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
        """Compute posterior quantities using information about event total/psi

        Parameters
        ----------
        df: xr.Dataset
            Data variables:
                event_passed[prefix, ec_idx]
                raw_total[prefix, ec_idx]
                raw_psi[prefix, ec_idx]
                bootstrap_total[prefix, ec_idx, bootstrap_replicate]
                bootstrap_psi[prefix, ec_idx, bootstrap_replicate]
            Coordinates:
                lsv_offsets[offset_idx]
                prefix[prefix]
            Derived (from _offsets):
                event_size[ec_idx]
                lsv_idx[ec_idx]
        events: xr.Dataset
            dataset that can be loaded along with matching introns/junctions as
            Events
        """
        offsets = df["lsv_offsets"].load().values
        if offsets[0] != 0:
            raise ValueError("offsets[0] must be zero")
        if offsets[-1] != df.sizes["ec_idx"]:
            raise ValueError("offsets[-1] must equal number of event connections")
        event_size = np.diff(offsets)
        if "event_size" not in df.variables:
            df = df.assign_coords(
                event_size=("ec_idx", np.repeat(event_size, event_size))
            )
        if "lsv_idx" not in df.variables:
            df = df.assign_coords(
                lsv_idx=("ec_idx", np.repeat(np.arange(len(event_size)), event_size)),
            )
        if df["event_passed"].dtype != bool:
            # for some reason, this is sometimes saved as int8, ensure consistent type
            df["event_passed"] = df["event_passed"].astype(bool)
        self.df: Final[xr.Dataset] = df
        self.events: Final[xr.Dataset] = events
        return

    @property
    def num_connections(self) -> int:
        return self.df.sizes["ec_idx"]

    @property
    def num_bootstraps(self) -> int:
        return self.df.sizes["bootstrap_replicate"]

    @property
    def num_prefixes(self) -> int:
        return self.df.sizes["prefix"]

    @property
    def prefixes(self) -> List[str]:
        return self.df["prefix"].values.tolist()

    def __repr__(self) -> str:
        MAX_PREFIXES_END = 1  # how many prefixes on either end to display
        print_prefixes_list = [
            *self.prefixes[:MAX_PREFIXES_END],
            *([] if self.num_prefixes <= 2 * MAX_PREFIXES_END else ["..."]),
            *self.prefixes[-MAX_PREFIXES_END:],
        ]
        print_prefixes = ", ".join(print_prefixes_list)
        return (
            f"PsiCoverage[{self.num_connections}]"
            f" for {self.num_prefixes} experiments [{print_prefixes}]"
        )

    def __getitem__(self, prefixes) -> "PsiCoverage":
        """Get subset of PsiCoverage corresponding to selected prefixes"""
        if isinstance(prefixes, str):
            # make sure that prefixes is a sequence
            prefixes = [prefixes]
        return PsiCoverage(self.df.sel(prefix=prefixes), self.events)

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
    def raw_coverage(self) -> xr.DataArray:
        """Infer LSV coverage at specific event connection from psi/total"""
        return self.raw_psi * self.raw_total

    @cached_property
    def bootstrap_coverage(self) -> xr.DataArray:
        """Infer LSV coverage at specific event connection from psi/total"""
        return self.bootstrap_psi * self.bootstrap_total

    @cached_property
    def alpha_prior(self) -> xr.DataArray:
        return 1 / self.event_size.astype(self.raw_psi.dtype)

    @cached_property
    def beta_prior(self) -> xr.DataArray:
        return 1 - self.alpha_prior

    @cached_property
    def raw_alpha(self) -> xr.DataArray:
        return (self.alpha_prior + self.raw_coverage).where(self.event_passed)

    @cached_property
    def bootstrap_alpha(self) -> xr.DataArray:
        return (self.alpha_prior + self.bootstrap_coverage).where(self.event_passed)

    @cached_property
    def raw_beta(self) -> xr.DataArray:
        return 1 + self.raw_total - self.raw_alpha

    @cached_property
    def bootstrap_beta(self) -> xr.DataArray:
        return 1 + self.bootstrap_total - self.bootstrap_alpha

    @cached_property
    def _approximate_params(self) -> Tuple[xr.DataArray, xr.DataArray]:
        _params = xr.apply_ufunc(
            bm.approximation,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            output_core_dims=[["_param_idx"]],
            dask="parallelized",
            dask_gufunc_kwargs=dict(output_sizes={"_param_idx": 2}),
        )
        alpha = _params.isel(_param_idx=0)
        beta = _params.isel(_param_idx=1)
        return (alpha, beta)

    @property
    def approximate_alpha(self) -> xr.DataArray:
        return self._approximate_params[0]

    @property
    def approximate_beta(self) -> xr.DataArray:
        return self._approximate_params[1]

    @cached_property
    def raw_posterior_mean(self) -> xr.DataArray:
        return self.raw_alpha / np.add(1, self.raw_total)

    @cached_property
    def raw_posterior_variance(self) -> xr.DataArray:
        mean = self.raw_posterior_mean
        return mean * np.subtract(1, mean) / np.add(1, self.raw_total)

    @cached_property
    def raw_posterior_std(self) -> xr.DataArray:
        return cast(xr.DataArray, np.sqrt(self.raw_posterior_variance))

    @property
    def raw_psi_mean(self) -> xr.DataArray:
        """alias for raw_posterior_mean"""
        return self.raw_posterior_mean

    @property
    def raw_psi_variance(self) -> xr.DataArray:
        """alias for raw_posterior_variance"""
        return self.raw_posterior_variance

    @property
    def raw_psi_std(self) -> xr.DataArray:
        """alias for raw_posterior_std"""
        return self.raw_posterior_std

    @cached_property
    def _bootstrap_moments(self) -> Tuple[xr.DataArray, xr.DataArray]:
        _moments = xr.apply_ufunc(
            bm.moments,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            output_core_dims=[["_moment_idx"]],
            dask="parallelized",
            dask_gufunc_kwargs=dict(output_sizes={"_moment_idx": 2}),
        )
        agg_mean = _moments.isel(_moment_idx=0)
        agg_variance = _moments.isel(_moment_idx=1)
        return (agg_mean, agg_variance)

    @property
    def bootstrap_posterior_mean(self) -> xr.DataArray:
        return self._bootstrap_moments[0]

    @property
    def bootstrap_posterior_variance(self) -> xr.DataArray:
        return self._bootstrap_moments[1]

    @cached_property
    def bootstrap_posterior_std(self) -> xr.DataArray:
        return cast(xr.DataArray, np.sqrt(self.bootstrap_posterior_variance))

    @property
    def bootstrap_psi_mean(self) -> xr.DataArray:
        """alias for bootstrap_posterior_mean"""
        return self.bootstrap_posterior_mean

    @property
    def bootstrap_psi_variance(self) -> xr.DataArray:
        """alias for bootstrap_posterior_variance"""
        return self.bootstrap_posterior_variance

    @property
    def bootstrap_psi_std(self) -> xr.DataArray:
        """alias for bootstrap_posterior_std"""
        return self.bootstrap_posterior_std

    @staticmethod
    def _compute_posterior_quantile(
        a: xr.DataArray,
        b: xr.DataArray,
        quantiles: Union[xr.DataArray, Sequence[float]] = [0.1, 0.9],
        mix_dim: str = "bootstrap_replicate",
    ) -> xr.DataArray:
        if not isinstance(quantiles, xr.DataArray):
            quantiles_arr = np.array(quantiles, dtype=a.dtype)
            if quantiles_arr.ndim > 1:
                raise ValueError(
                    "Unable to handle non-xarray multi-dimensional quantiles"
                )
            elif quantiles_arr.ndim == 0:
                quantiles_arr = quantiles_arr[np.newaxis]
            quantiles = xr.DataArray(quantiles_arr, [("quantiles", quantiles_arr)])
        # if mixture dimension is not present, treat as one-component mixture
        if mix_dim not in a.dims:
            a = a.expand_dims(**{mix_dim: 1})
        if mix_dim not in b.dims:
            b = b.expand_dims(**{mix_dim: 1})
        return xr.apply_ufunc(
            bm.quantile,
            quantiles,
            a,
            b,
            input_core_dims=[[], [mix_dim], [mix_dim]],
            dask="allowed",
        )

    def bootstrap_quantile(
        self,
        quantiles: Union[xr.DataArray, Sequence[float]] = [0.1, 0.9],
    ) -> xr.DataArray:
        """Get quantiles of bootstrap posterior"""
        return self._compute_posterior_quantile(
            self.bootstrap_alpha, self.bootstrap_beta, quantiles=quantiles
        )

    def approximate_quantile(
        self,
        quantiles: Union[xr.DataArray, Sequence[float]] = [0.1, 0.9],
    ) -> xr.DataArray:
        """Get quantiles of approximate bootstrap posterior"""
        return self._compute_posterior_quantile(
            self.approximate_alpha, self.approximate_beta, quantiles=quantiles
        )

    @staticmethod
    def _compute_posterior_discretized_pmf(
        a: xr.DataArray,
        b: xr.DataArray,
        nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
        mix_dim: str = "bootstrap_replicate",
    ) -> xr.DataArray:
        endpoints = np.linspace(0, 1, 1 + nbins, dtype=a.dtype)
        dummy_bins = xr.DataArray(
            np.empty(nbins, dtype=a.dtype),
            {
                "pmf_bin_start": ("pmf_bin", endpoints[:-1]),
                "pmf_bin_end": ("pmf_bin", endpoints[1:]),
            },
            dims=["pmf_bin"],
        )
        # if mixture dimension is not present, treat as one-component mixture
        if mix_dim not in a.dims:
            a = a.expand_dims(**{mix_dim: 1})
        if mix_dim not in b.dims:
            b = b.expand_dims(**{mix_dim: 1})
        return xr.apply_ufunc(
            bm.pmf,
            a,
            b,
            dummy_bins,
            input_core_dims=[
                [mix_dim],
                [mix_dim],
                ["pmf_bin"],
            ],
            output_core_dims=[["pmf_bin"]],
            dask="allowed",
        )

    def bootstrap_discretized_pmf(
        self, nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS
    ) -> xr.DataArray:
        """Get discretized PMF with specified number of bins"""
        return self._compute_posterior_discretized_pmf(
            self.bootstrap_alpha, self.bootstrap_beta, nbins=nbins
        )

    def approximate_discretized_pmf(
        self, nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS
    ) -> xr.DataArray:
        """Get discretized PMF with specified number of bins"""
        return self._compute_posterior_discretized_pmf(
            self.approximate_alpha, self.approximate_beta, nbins=nbins
        )

    @cached_property
    def _raw_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        return self.raw_psi_mean.chunk({"prefix": None})

    @cached_property
    def _bootstrap_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        return self.bootstrap_psi_mean.chunk({"prefix": None})

    @cached_property
    def raw_psi_mean_population_median(self) -> xr.DataArray:
        return self._raw_psi_mean_core_prefix.median("prefix")

    @cached_property
    def bootstrap_psi_mean_population_median(self) -> xr.DataArray:
        return self._bootstrap_psi_mean_core_prefix.median("prefix")

    def raw_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
    ) -> xr.DataArray:
        """Get quantiles of psi mean over population (adds dim quantile)"""
        return self._raw_psi_mean_core_prefix.quantile(quantiles, "prefix").rename(
            quantile="population_quantile"
        )

    def bootstrap_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
    ) -> xr.DataArray:
        """Get quantiles of psi mean over population (adds dim quantile)"""
        return self._bootstrap_psi_mean_core_prefix.quantile(
            quantiles, "prefix"
        ).rename(quantile="population_quantile")

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
        # get whether individual connection passes thresholds
        passed = (events_coverage.numreads >= minreads) & (
            events_coverage.numbins >= minbins
        )
        # get whether any connection in event passed, per connection
        event_passed = _offsets.offset_logical_or(passed, offsets)
        # get total coverage per event, per connection
        raw_total = _offsets.offsetsum(
            events_coverage.numreads, offsets, axes=[0, -1, 0]
        )
        bootstrap_total = _offsets.offsetsum(
            events_coverage.bootstraps, offsets, axes=[0, -1, 0]
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
            ).expand_dims(prefix=[bam_experiment_name(events_coverage.bam_path)]),
            events_coverage.events.save_df,
        )

    @classmethod
    def from_zarr(cls, path: Union[str, Path, List[Union[str, Path]]]) -> "PsiCoverage":
        """Load one or more PsiCoverage files together at once

        Load one or more PsiCoverage files together at once. If they have
        overlapping prefixes, data will be loaded from the first file with the
        given prefix.
        """
        if not isinstance(path, list):
            path = [path]
        df = xr.open_mfdataset(
            path,
            engine="zarr",
            group=constants.NC_PSICOVERAGE,
            combine="nested",
            concat_dim="prefix",
            join="override",
            compat="override",
            coords="minimal",
            data_vars="minimal",
        )
        if len(path) > 1:
            # attributes are defined by path[0]. We'd rather just have none
            df.attrs.clear()
        events_df = xr.open_zarr(path[0], group=constants.NC_EVENTS)
        return cls(df, events_df)

    def updated(
        self,
        bootstrap_psi: Optional[xr.DataArray],
        raw_psi: Optional[xr.DataArray],
        **update_attrs,
    ) -> "PsiCoverage":
        """Updated PsiCoverage, replacing bootstrap_psi, raw_psi"""
        df = self.df
        # update psi arrays
        if bootstrap_psi is not None:
            if set(self.bootstrap_psi.dims) != set(bootstrap_psi.dims):
                raise ValueError("bootstrap_psi doesn't have same named axes")
            df = df.assign(bootstrap_psi=bootstrap_psi)
        if raw_psi is not None:
            if set(self.raw_psi.dims) != set(raw_psi.dims):
                raise ValueError("raw_psi doesn't have same named axes")
            df = df.assign(raw_psi=raw_psi)
        # update/add attributes
        df = df.assign_attrs(**update_attrs)
        # return resulting PsiCoverage object
        return PsiCoverage(df, self.events)

    def to_zarr(
        self,
        path: Union[str, Path],
        ec_chunksize: int = 8192,
        append: bool = False,
        consolidated: bool = True,
    ) -> None:
        """Save PSI coverage dataset as zarr

        Parameters
        ----------
        path: Union[str, Path]
            Path for output file in zarr format
        ec_chunksize: int
            How to chunk event connections to prevent memory from getting to
            large when loading many samples simultaneously
        append: bool
            Add to *existing* file (not checked). But if append is used on
            non-existing file, event information will not be saved.
        consolidated: bool
            When saving the file make sure that it is consolidated. In general,
            if you are appending a bunch of files together, it can make sense
            to set consolidated=False, and consolidate on the last write (only
            consolidate once). But, don't forget to consolidate at the end.

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
        save_df = self.df.drop_vars(["event_size", "lsv_idx"])
        if save_df.sizes["ec_idx"] > 0:
            save_df = save_df.chunk(USE_CHUNKS)  # type: ignore[arg-type]
        # clear any previous encodings from before, before saving
        for v in save_df.variables.values():
            v.encoding.clear()
        if append:
            # remove attributes referring to bam
            for k in [x for x in save_df.attrs.keys() if str(x).startswith("bam")]:
                save_df.attrs.pop(k, None)
            save_df.to_zarr(
                path,
                append_dim="prefix",
                group=constants.NC_PSICOVERAGE,
                consolidated=consolidated,
            )
        else:
            save_df.to_zarr(
                path,
                mode="w",
                group=constants.NC_PSICOVERAGE,
                consolidated=False,
            )
            self.events.to_zarr(
                path, mode="a", group=constants.NC_EVENTS, consolidated=consolidated
            )
        return

    @cached_property
    def num_passed(self) -> xr.DataArray:
        return self.event_passed.sum("prefix")

    def passed_min_experiments(
        self,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
    ) -> xr.DataArray:
        """Get boolean mask of events that pass enough experiments"""
        return self.num_passed >= min_experiments(min_experiments_f, self.num_prefixes)

    def sum(
        self,
        new_prefix: str,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
    ) -> "PsiCoverage":
        """Aggregate coverage/psi values over all prefixes"""
        event_passed = self.passed_min_experiments(min_experiments_f)
        raw_total = self.raw_total.sum("prefix")
        raw_coverage = (self.raw_total * self.raw_psi).sum("prefix")
        raw_psi = (raw_coverage / raw_total.where(raw_total > 0)).fillna(0)
        bootstrap_total = self.bootstrap_total.sum("prefix")
        bootstrap_coverage = (self.bootstrap_total * self.bootstrap_psi).sum("prefix")
        bootstrap_psi = (
            bootstrap_coverage / bootstrap_total.where(bootstrap_total > 0)
        ).fillna(0)
        df = xr.Dataset(
            data_vars=dict(
                event_passed=event_passed,
                raw_total=raw_total,
                raw_psi=raw_psi,
                bootstrap_total=bootstrap_total,
                bootstrap_psi=bootstrap_psi,
            ),
            coords=dict(
                lsv_offsets=self.lsv_offsets,
                event_size=self.event_size,
                lsv_idx=self.lsv_idx,
            ),
            attrs=dict(original_prefix=self.prefixes),
        ).expand_dims(prefix=[new_prefix])
        return PsiCoverage(df, self.events)

    def mask_events(self, passed: xr.DataArray) -> "PsiCoverage":
        """Return PsiCoverage passing only events that are passed in input

        Return PsiCoverage passing only events that are passed in input (and in
        the original object)
        """
        return PsiCoverage(
            self.df.assign(event_passed=self.event_passed & passed), self.events
        )

    def drop_unquantifiable(self) -> "PsiCoverage":
        """Drop all events that are not (passed in all prefixes)"""
        # what passed?
        ec_idx_passed = self.event_passed.all("prefix").load().reset_coords(drop=True)
        e_idx_passed = (
            ec_idx_passed.isel(ec_idx=self.lsv_offsets.values[:-1])
            .rename(ec_idx="e_idx")
            .reset_coords(drop=True)
        )
        # get subset of df, events corresponding to these, dropping variables
        # that require recalculation
        df_subset = self.df.drop_vars(["lsv_offsets", "lsv_idx"]).sel(
            ec_idx=ec_idx_passed
        )
        events_subset = self.events.drop_vars(["_offsets"]).sel(
            ec_idx=ec_idx_passed, e_idx=e_idx_passed
        )
        # recalculate offsets
        passed_sizes = (
            self.event_size.isel(ec_idx=self.lsv_offsets.values[:-1])
            .rename(ec_idx="e_idx")
            .sel(e_idx=e_idx_passed)
            .values
        )
        offsets = np.empty(len(passed_sizes) + 1, dtype=passed_sizes.dtype)
        offsets[0] = 0
        offsets[1:] = np.cumsum(passed_sizes)
        # add offsets back in
        df_subset = df_subset.assign_coords(lsv_offsets=("offset_idx", offsets))
        events_subset = events_subset.assign_coords(
            _offsets=("e_offsets_idx", offsets.astype(np.uint64))
        )
        # return subsetted PsiCoverage
        return PsiCoverage(df_subset, events_subset)

    def get_events(
        self,
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> Events:
        if self.events.intron_hash != introns.checksum():
            raise ValueError("GeneIntrons checksums do not match")
        if self.events.junction_hash != junctions.checksum():
            raise ValueError("GeneJunctions checksums do not match")
        return Events(
            _Events(
                introns._gene_introns,
                junctions._gene_junctions,
                self.events.ref_exon_idx,
                self.events.event_type,
                self.events._offsets,
                self.events.is_intron,
                self.events.connection_idx,
            )
        )

    def dataset(
        self,
        properties: Sequence[str] = [
            "raw_psi_mean",
            "raw_psi_std",
            "bootstrap_psi_mean",
            "bootstrap_psi_std",
            "raw_coverage",
        ],
        quantiles: Sequence[float] = list(),
        psibins: Optional[int] = None,
        use_posterior: str = "approximation",
    ) -> xr.Dataset:
        """Extract selected properties into single dataset

        Parameters
        ----------
        properties: Sequence[str]
            PsiCoverage properties to request.
        quantiles: Sequence[float]
            If non-empty, calculate quantiles of posterior distribution
        psibins: Optional[int]
            If specified, calculate discretized approximation to posterior
            distribution with this many bins
        use_posterior: str
            Compute quantiles/pmf with "bootstrap" or "approximation" posterior
            distribution (or "both"). Otherwise, raise error.
        """
        if use_posterior not in constants.PSICOV_POSTERIORS:
            raise ValueError(
                f"{use_posterior = } must be one of {constants.PSICOV_POSTERIORS}"
            )
        # initialize variables to return with noting if any experiment passed
        quantify_vars: Dict[str, xr.DataArray] = {
            "any_passed": self.event_passed.any("prefix")
        }
        # add properties
        for x in properties:
            quantify_vars[x] = getattr(self, x)
        if len(quantiles) or psibins:
            if len(quantiles):
                if use_posterior in constants.PSICOV_APPROX:
                    quantify_vars["approx_psi_quantile"] = self.approximate_quantile(
                        quantiles
                    )
                if use_posterior in constants.PSICOV_BOOTSTRAP:
                    quantify_vars["bootstrap_psi_quantile"] = self.bootstrap_quantile(
                        quantiles
                    )
            if psibins:
                if use_posterior in constants.PSICOV_APPROX:
                    quantify_vars["approx_psi_pmf"] = self.approximate_discretized_pmf(
                        psibins
                    )
                if use_posterior in constants.PSICOV_BOOTSTRAP:
                    quantify_vars["bootstrap_psi_pmf"] = self.bootstrap_discretized_pmf(
                        psibins
                    )
        return xr.Dataset(quantify_vars).reset_coords(drop=True)  # type: ignore[arg-type]
