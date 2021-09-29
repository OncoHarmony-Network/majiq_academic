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
from typing import Final, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import xarray as xr

import new_majiq._offsets as _offsets
import new_majiq.beta_mixture as bm
import new_majiq.constants as constants
from new_majiq.Events import Events, _Events
from new_majiq.EventsCoverage import EventsCoverage
from new_majiq.experiments import bam_experiment_name
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.SpliceGraph import SpliceGraph


def min_experiments(min_experiments_f: float, num_experiments: int) -> float:
    if min_experiments_f < 1:
        min_experiments_f *= num_experiments
    return max(1, min(min_experiments_f, num_experiments))


def _silent_quantile(*args, **kwargs):
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return bm.quantile(*args, **kwargs)


def _silent_cdf(*args, **kwargs):
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return bm.cdf(*args, **kwargs)


class PsiCoverage(object):
    """PSI and total coverage (raw and bootstrapped) for arbitrary number of samples"""

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
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
        return np.reciprocal(self.event_size.astype(self.raw_psi.dtype))

    @cached_property
    def beta_prior(self) -> xr.DataArray:
        return np.subtract(1, self.alpha_prior)

    @cached_property
    def raw_alpha(self) -> xr.DataArray:
        return self.alpha_prior + self.raw_coverage

    @cached_property
    def bootstrap_alpha(self) -> xr.DataArray:
        return self.alpha_prior + self.bootstrap_coverage

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
        variance = mean * np.subtract(1, mean) / np.add(2, self.bootstrap_total)
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
            pmf_bin_start=("pmf_bin", endpoints.values[:-1]),
            pmf_bin_end=("pmf_bin", endpoints.values[1:]),
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

    @property
    def prefixes(self) -> List[str]:
        return self.df["prefix"].values.tolist()

    def sum(
        self,
        new_prefix: str,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
    ) -> "PsiCoverage":
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

    @staticmethod
    def _exons_formatted(
        exon_start: Sequence[int], exon_end: Sequence[int]
    ) -> List[str]:
        def format_coord(x):
            return x if x >= 0 else "na"

        return [
            f"{format_coord(a)}-{format_coord(b)}" for a, b in zip(exon_start, exon_end)
        ]

    def as_dataframe(self, sg: SpliceGraph) -> pd.DataFrame:
        """Basic quantifications annotated with splicegraph information"""
        # set up events object
        self.events.load()
        q_events = self.get_events(sg.introns, sg.junctions)
        # get event-level information
        event_size = np.diff(self.lsv_offsets.values)
        event_id = sg.exon_connections.event_id(
            q_events.ref_exon_idx, q_events.event_type
        )
        event_description = sg.exon_connections.event_description(
            q_events.ref_exon_idx, q_events.event_type
        )
        ref_exon_start = sg.exons.start[q_events.ref_exon_idx]
        ref_exon_end = sg.exons.end[q_events.ref_exon_idx]
        # connection information
        gene_idx = q_events.connection_gene_idx()
        gene_id = np.array(sg.genes.gene_id)[gene_idx]
        gene_name = np.array(sg.genes.gene_name)[gene_idx]
        strand = np.array([x.decode() for x in sg.genes.strand])[gene_idx]
        contig_idx = sg.genes.contig_idx[gene_idx]
        seqid = np.array(sg.contigs.seqid)[contig_idx]
        other_exon_idx = q_events.connection_other_exon_idx()
        other_exon_start = sg.exons.start[other_exon_idx]
        other_exon_end = sg.exons.end[other_exon_idx]
        # as Dataset
        df = xr.Dataset(
            {
                "raw_psi_mean": self.raw_posterior_mean,
                "raw_psi_std": np.sqrt(self.raw_posterior_variance),
                "bootstrap_psi_mean": self.bootstrap_posterior_mean,
                "bootstrap_psi_std": np.sqrt(self.bootstrap_posterior_variance),
                "seqid": ("ec_idx", seqid),
                "gene_id": ("ec_idx", gene_id),
                "ref_exon": (
                    "ec_idx",
                    np.repeat(
                        self._exons_formatted(ref_exon_start, ref_exon_end), event_size
                    ),
                ),
                "event_type": (
                    "ec_idx",
                    np.repeat([x.decode() for x in q_events.event_type], event_size),
                ),
                "is_intron": ("ec_idx", q_events.is_intron),
                "start": ("ec_idx", q_events.connection_start()),
                "end": ("ec_idx", q_events.connection_end()),
                "denovo": ("ec_idx", q_events.connection_denovo()),
                "gene_name": ("ec_idx", gene_name),
                "strand": ("ec_idx", strand),
                "other_exon": (
                    "ec_idx",
                    self._exons_formatted(other_exon_start, other_exon_end),
                ),
                "event_id": ("ec_idx", np.repeat(event_id, event_size)),
                "event_description": (
                    "ec_idx",
                    np.repeat(event_description, event_size),
                ),
            }
        ).reset_coords(drop=True)
        # if only one prefix, drop it
        idx_order = ["prefix", "ec_idx"]
        try:
            df = df.squeeze("prefix", drop=True)
            idx_order.remove("prefix")
        except ValueError:
            pass
        return df.to_dataframe(idx_order)  # type: ignore

    def quantifier_dataset(
        self,
        pmf_bins: Optional[int] = 40,
        quantiles: Optional[Sequence[float]] = None,
    ) -> xr.Dataset:
        """Default dataset for quantifications"""
        df = xr.Dataset(
            {
                "raw_psi_mean": self.raw_posterior_mean.load(),
                "raw_psi_std": np.sqrt(self.raw_posterior_variance.load()),
                "bootstrap_psi_mean": self.bootstrap_posterior_mean.load(),
                "bootstrap_psi_std": np.sqrt(self.bootstrap_posterior_variance.load()),
            },
            {},
            self.df.attrs,
        )
        if pmf_bins:
            df = df.assign(psi_pmf=self.bootstrap_discretized_pmf(pmf_bins))
        if quantiles:
            df = df.assign(psi_quantiles=self.bootstrap_quantile(quantiles))
        return df
