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
from typing import (
    Any,
    Dict,
    Final,
    Hashable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

import dask.array as da
import numpy as np
import xarray as xr
from dask.delayed import Delayed
from dask.distributed import progress

import new_majiq._offsets as _offsets
import new_majiq.beta_mixture as bm
import new_majiq.constants as constants
from new_majiq._stats import nanmedian, nanquantile
from new_majiq.experiments import bam_experiment_name

from .Events import Events, _Events
from .EventsCoverage import EventsCoverage
from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions
from .SJExperiment import SJExperiment


def min_experiments(min_experiments_f: float, num_experiments: int) -> float:
    if min_experiments_f < 1:
        min_experiments_f *= num_experiments
    return max(1, min(min_experiments_f, num_experiments))


class PsiCoverage(object):
    """Coverage for Psi posterior distributions over raw or bootstrap coverage

    Coverage for Psi posterior distributions over raw or bootstrap coverage.
    Holds raw and bootstrap coverage over events per prefix (e.g. experiment,
    summary over experiments, etc.) and whether they passed quantification
    thresholds.
    Functions/attributes allow computation with underlying posterior
    distributions.

    Parameters
    ----------
    df: xr.Dataset
        Data variables:
            - event_passed[prefix, ec_idx]
            - raw_total[prefix, ec_idx]
            - raw_psi[prefix, ec_idx]
            - bootstrap_total[prefix, ec_idx, bootstrap_replicate]
            - bootstrap_psi[prefix, ec_idx, bootstrap_replicate]
        Coordinates:
            - lsv_offsets[offset_idx]
            - prefix[prefix]
        Derived (from _offsets):
            - event_size[ec_idx]
            - lsv_idx[ec_idx]
    events: xr.Dataset
        dataset that can be loaded along with matching introns/junctions as
        Events

    See Also
    --------
    PsiCoverage.from_sj_lsvs
    PsiCoverage.from_events_coverage
    PsiCoverage.from_zarr
    PsiCoverage.updated
    PsiCoverage.sum
    PsiCoverage.mask_events
    PsiCoverage.__getitem__
    """

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
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
        if self.num_prefixes > 2 * MAX_PREFIXES_END:
            print_prefixes_list = [
                *self.prefixes[:MAX_PREFIXES_END],
                *([] if self.num_prefixes <= 2 * MAX_PREFIXES_END else ["..."]),
                *self.prefixes[-MAX_PREFIXES_END:],
            ]
        else:
            print_prefixes_list = self.prefixes
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

    @cached_property
    def bootstrap_psi_mean_legacy(self) -> xr.DataArray:
        """old calculation of bootstrap psi mean where summarizing by median"""
        return xr.apply_ufunc(
            bm.means_median,
            self.bootstrap_alpha,
            self.bootstrap_beta,
            input_core_dims=[["bootstrap_replicate"], ["bootstrap_replicate"]],
            dask="parallelized",
        )

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
        result = self.raw_psi_mean
        if result.chunks:
            result = result.chunk({"prefix": None})
        return result

    @cached_property
    def _bootstrap_psi_mean_core_prefix(self) -> xr.DataArray:
        """For computing quantiles over a population of samples"""
        result = self.bootstrap_psi_mean
        if result.chunks:
            result = result.chunk({"prefix": None})
        return result

    @cached_property
    def raw_psi_mean_population_median(self) -> xr.DataArray:
        return xr.apply_ufunc(
            nanmedian,
            self._raw_psi_mean_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    @cached_property
    def bootstrap_psi_mean_population_median(self) -> xr.DataArray:
        return xr.apply_ufunc(
            nanmedian,
            self._bootstrap_psi_mean_core_prefix,
            input_core_dims=[["prefix"]],
            dask="allowed",
        )

    @staticmethod
    def _compute_population_quantile(
        x: xr.DataArray,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
    ) -> xr.DataArray:
        quantiles_xr = xr.DataArray(quantiles, [(quantile_dim_name, quantiles)])
        return xr.apply_ufunc(
            nanquantile,
            x,
            quantiles_xr,
            input_core_dims=[["prefix"], [quantile_dim_name]],
            output_core_dims=[[quantile_dim_name]],
            dask="allowed",
        )

    def raw_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
    ) -> xr.DataArray:
        """Get quantiles of psi mean over population (adds dim quantile)"""
        return self._compute_population_quantile(
            self._raw_psi_mean_core_prefix,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )

    def bootstrap_psi_mean_population_quantile(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        quantile_dim_name: str = "population_quantile",
    ) -> xr.DataArray:
        """Get quantiles of psi mean over population (adds dim quantile)"""
        return self._compute_population_quantile(
            self._bootstrap_psi_mean_core_prefix,
            quantiles,
            quantile_dim_name=quantile_dim_name,
        )

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
    def from_sj_lsvs(
        cls,
        sj: SJExperiment,
        lsvs: Events,
        minreads: float = constants.DEFAULT_QUANTIFY_MINREADS,
        minbins: float = constants.DEFAULT_QUANTIFY_MINBINS,
        num_bootstraps: int = constants.DEFAULT_COVERAGE_NUM_BOOTSTRAPS,
        pvalue_threshold: float = constants.DEFAULT_COVERAGE_STACK_PVALUE,
    ) -> "PsiCoverage":
        """Given SJ and Events information (generally LSVs), load PsiCoverage"""
        lsv_coverage = EventsCoverage.from_events_and_sj(
            lsvs, sj, num_bootstraps=num_bootstraps, pvalue_threshold=pvalue_threshold
        )
        return PsiCoverage.from_events_coverage(lsv_coverage, minreads, minbins)

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

    def _save_df(
        self,
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        remove_bam_attrs: bool = False,
    ) -> xr.Dataset:
        """Prepare dataset of psicoverage that will be saved

        Prepare dataset of psicoverage that will be saved. This sets the
        chunking, clears any encodings from before, removes attrs that have to
        do with BAM, etc.
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
        # remove attributes referring to bam?
        if remove_bam_attrs:
            for k in [x for x in save_df.attrs.keys() if str(x).startswith("bam")]:
                save_df.attrs.pop(k, None)
        return save_df

    def to_zarr(
        self,
        path: Union[str, Path],
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        consolidated: bool = True,
        show_progress: bool = False,
    ) -> None:
        """Save PSI coverage dataset as zarr

        Parameters
        ----------
        path: Union[str, Path]
            Path for output file in zarr format
        ec_chunksize: int
            How to chunk event connections to prevent memory from getting to
            large when loading many samples simultaneously
        consolidated: bool
            When saving the file make sure that it is consolidated. In general,
            if you are appending a bunch of files together, it can make sense
            to set consolidated=False, and consolidate on the last write (only
            consolidate once). But, don't forget to consolidate at the end.
        show_progress: bool
            Attempt to show progress on distributed cluster for Dask
        """
        save_df = self._save_df(ec_chunksize=ec_chunksize)
        save_df_future = cast(
            Delayed,
            save_df.to_zarr(
                path,
                mode="w",
                group=constants.NC_PSICOVERAGE,
                consolidated=False,
                compute=False,
            ),
        )
        if show_progress:
            save_df_future = save_df_future.persist()
            progress(save_df_future)
        else:
            save_df_future.compute()
        self.events.chunk(self.events.sizes).to_zarr(
            path, mode="a", group=constants.NC_EVENTS, consolidated=consolidated
        )
        return

    def to_zarr_slice(
        self,
        path: Union[str, Path],
        prefix_slice: slice,
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
    ) -> None:
        """Save PsiCoverage to specified path for specified slice on prefix

        Save PsiCoverage to specified path for specified slice. Typically run
        after PsiCoverage.to_zarr_slice_init()
        """
        self._save_df(ec_chunksize=ec_chunksize).drop_vars("prefix").pipe(
            lambda x: x.drop_vars(
                [k for k, v in x.variables.items() if "prefix" not in v.dims]
            )
        ).to_zarr(
            path, group=constants.NC_PSICOVERAGE, region=dict(prefix=prefix_slice)
        )
        return

    @classmethod
    def to_zarr_slice_init(
        cls,
        path: Union[str, Path],
        events_df: xr.Dataset,
        prefixes: List[str],
        num_bootstraps: int,
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        cov_dtype: type = np.float32,
        psicov_attrs: Dict[Hashable, Any] = dict(),
    ) -> None:
        """Init zarr for PsiCoverage over many prefixes for multithreaded write

        Initialize zarr for PsiCoverage over many prefixes. Saves all
        information except dimensions that are prefix-specific. This enables
        multithreaded (or multiprocess) write with to_zarr_slice

        Parameters
        ----------
        path: Union[str, Path]
            Path for output Zarr for psicoverage output
        events_df: xr.Dataset
            Dataset encoding psicoverage events (Events.save_df,
            PsiCoverage.events, etc.)
        prefixes: List[str]
            Values for the prefix dimension coordinate
        num_prefixes: int
            Number of prefixes that will be filled in
        num_bootstraps: int
            Number of bootstrap replicates that will be used
        ec_chunksize: int
            How to chunk event connections to prevent memory from getting to
            large when loading many samples simultaneously
        cov_dtype: type
            What type to use for psi/total_coverage arrays
        psicov_attrs: Dict[Hashable, Any]
            Attributes to include
        """
        # force events to be saved as single chunk (no benefit for chunking here)
        events_df = events_df.chunk(events_df.sizes)
        # save events
        events_df.to_zarr(path, mode="w", group=constants.NC_EVENTS, consolidated=False)
        # dims for skeleton
        raw_dims = ("prefix", "ec_idx")
        bootstrap_dims = (*raw_dims, "bootstrap_replicate")
        # shapes for skeleton
        raw_shape = (len(prefixes), events_df.sizes["ec_idx"])
        bootstrap_shape = (*raw_shape, num_bootstraps)
        # chunksizes for skeleton
        raw_chunks = (1, ec_chunksize)
        bootstrap_chunks = (*raw_chunks, None)
        # arrays for skeleton
        raw_arr = da.empty(raw_shape, dtype=cov_dtype, chunks=raw_chunks)
        passed_arr = da.empty(raw_shape, dtype=bool, chunks=raw_chunks)
        bootstrap_arr = da.empty(
            bootstrap_shape, dtype=cov_dtype, chunks=bootstrap_chunks
        )
        # save metadata for skeleton
        xr.Dataset(
            dict(
                bootstrap_psi=(bootstrap_dims, bootstrap_arr),
                bootstrap_total=(bootstrap_dims, bootstrap_arr),
                event_passed=(raw_dims, passed_arr),
                raw_psi=(raw_dims, raw_arr),
                raw_total=(raw_dims, raw_arr),
            ),
        ).to_zarr(
            path,
            mode="a",
            compute=False,
            group=constants.NC_PSICOVERAGE,
            consolidated=False,
        )
        # save offsets, prefixes, and attributes
        add_offsets = (
            events_df[["_offsets"]]
            .reset_coords()
            .astype(int)
            .set_coords("_offsets")
            .rename_dims(e_offsets_idx="offset_idx")
            .rename_vars(_offsets="lsv_offsets")
            .assign_coords(prefix=("prefix", prefixes))
        )
        add_offsets.attrs = psicov_attrs  # overwrite attrs with what we want
        add_offsets.to_zarr(
            path, mode="a", group=constants.NC_PSICOVERAGE, consolidated=True
        )

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
            "bootstrap_psi_mean_legacy",
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
