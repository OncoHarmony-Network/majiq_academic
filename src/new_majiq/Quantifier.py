"""
Quantifier.py

TODO: make this able to accept EventsCoverage from memory, not just saved to
disk?

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr
import pandas as pd

import new_majiq.constants as constants
import new_majiq.beta_mixture as bm

from new_majiq.SpliceGraph import SpliceGraph
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.Events import Events, _Events

from functools import cached_property
from typing import (
    Final,
    List,
    NamedTuple,
    Optional,
    Sequence,
    Tuple,
    Union,
)
from pathlib import Path


class QuantifierThresholds(NamedTuple):
    minreads: float = constants.DEFAULT_QUANTIFY_MINREADS
    minbins: float = constants.DEFAULT_QUANTIFY_MINBINS
    min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS


class ConnectionsChecksum(NamedTuple):
    introns: int
    junctions: int


def min_experiments(min_experiments_f: float, num_experiments: int) -> float:
    if min_experiments_f < 1:
        min_experiments_f *= num_experiments
    return max(1, min(min_experiments_f, num_experiments))


class QuantifiableEvents(object):
    """Process group of EventsCoverage for events that can be quantified"""

    def __init__(
        self,
        checksums: ConnectionsChecksum,
        offsets: np.ndarray,
        event_passed: np.ndarray,
    ):
        """Indicate events quantifiable

        Parameters
        ----------
        checksums: ConnectionsChecksum
            checksums from junctions and introns to quickly check that events
            are same
        offsets: np.ndarray
            Offsets for events into event connections (1 longer than events)
        event_passed: np.ndarray
            boolean mask over events indicating if the event is quantifiable
        """
        if offsets.ndim != 1:
            raise ValueError("offsets must be 1D")
        if event_passed.ndim != 1:
            raise ValueError("event_passed must be 1D")
        if len(offsets) != len(event_passed) + 1:
            raise ValueError("offsets and event_passed do not have concordant sizes")

        self._checksums: Final[ConnectionsChecksum] = checksums
        self._offsets: Final[np.ndarray] = offsets
        self._event_passed: Final[np.ndarray] = event_passed
        return

    @classmethod
    def from_quantifier_group(
        cls,
        experiments: Sequence[Union[Path, str]],
        thresholds: QuantifierThresholds = QuantifierThresholds(),
    ) -> "QuantifiableEvents":
        """Determine quantifiable events from input EventsCoverage paths"""
        if len(experiments) == 0:
            raise ValueError("No experiments passed into quantifiable events")
        checksums: ConnectionsChecksum
        offsets: np.ndarray
        passed_ct: np.ndarray
        for path in experiments:
            # load events to get offsets and checksums
            with xr.open_dataset(path, group=constants.NC_EVENTS) as df:
                cur_checksums: ConnectionsChecksum = ConnectionsChecksum(
                    introns=df.intron_hash, junctions=df.junction_hash
                )
                try:
                    if checksums != cur_checksums:
                        raise ValueError(
                            f"{path} has different intron/junction checksums"
                            f" than preceding inputs in {experiments}"
                        )
                except NameError:
                    checksums = cur_checksums
                    offsets = df._offsets.values.astype(int)
                # so checksums are good, and we have now defined offsets
            with xr.open_dataset(path, group=constants.NC_EVENTSCOVERAGE) as df:
                cur_summaries = df[["numreads", "numbins"]].load()
            # determine if each event connection passed experiment thresholds
            passed = (cur_summaries.numreads.values >= thresholds.minreads) & (
                cur_summaries.numbins.values >= thresholds.minbins
            )
            try:
                passed_ct += passed
            except NameError:  # this is the first experiment
                passed_ct = passed.astype(int)
        # which connections passed?
        connection_passed = passed_ct >= min_experiments(
            thresholds.min_experiments_f, len(experiments)
        )
        # event passes if any of its connections passed
        event_passed: np.ndarray = np.logical_or.reduceat(connection_passed, offsets[:-1])
        return QuantifiableEvents(checksums, offsets, event_passed)

    @property
    def checksums(self) -> ConnectionsChecksum:
        return self._checksums

    @property
    def offsets(self) -> np.ndarray:
        return self._offsets

    @property
    def quantifiable_offsets(self) -> np.ndarray:
        """If we subset the quantifiable events, these are the new offsets"""
        quantifiable_event_sizes = (
            self.offsets[1:][self.event_passed] - self.offsets[:-1][self.event_passed]
        )
        # offsets are 1 longer than this
        quantifiable_offsets = np.empty(1 + len(quantifiable_event_sizes), dtype=int)
        quantifiable_offsets[0] = 0
        quantifiable_offsets[1:] = np.cumsum(quantifiable_event_sizes)
        return quantifiable_offsets

    @property
    def event_passed(self) -> np.ndarray:
        return self._event_passed

    @cached_property
    def event_connection_passed_mask(self) -> np.ndarray:
        """Propagate passed events back to mask over event connections"""
        return np.repeat(self.event_passed, np.diff(self.offsets))

    def __and__(self, other: "QuantifiableEvents") -> "QuantifiableEvents":
        """Get events that are quantifiable in both"""
        if self.checksums != other.checksums:
            raise ValueError("Cannot AND quantifiable events that are not shared")
        return QuantifiableEvents(
            self.checksums, self.offsets, self.event_passed & other.event_passed
        )


class QuantifiableCoverage(object):
    """Coverage for a group of experiments at quantifiable events"""

    def __init__(
        self,
        offsets: np.ndarray,
        numreads: np.ndarray,
        bootstraps: np.ndarray,
        events: xr.Dataset,
    ):
        self._offsets: Final[np.ndarray] = offsets
        self._numreads: Final[np.ndarray] = numreads
        self._bootstraps: Final[np.ndarray] = bootstraps
        self._events: Final[xr.Dataset] = events
        # TODO check inputs
        return

    @property
    def offsets(self) -> np.ndarray:
        return self._offsets

    @cached_property
    def event_size(self) -> np.ndarray:
        return np.diff(self.offsets)

    @cached_property
    def connection_event_size(self) -> np.ndarray:
        return np.repeat(self.event_size, self.event_size)

    @cached_property
    def alpha_prior(self) -> np.ndarray:
        return np.reciprocal(self.connection_event_size, dtype=self.numreads.dtype)

    @cached_property
    def beta_prior(self) -> np.ndarray:
        return 1 - self.alpha_prior

    @property
    def events(self) -> xr.Dataset:
        return self._events

    @property
    def numreads(self) -> np.ndarray:
        return self._numreads

    @property
    def bootstraps(self) -> np.ndarray:
        return self._bootstraps

    def _total_coverage(self, coverage: np.ndarray) -> np.ndarray:
        event_total = np.add.reduceat(coverage, self.offsets[:-1], axis=0)
        return np.repeat(event_total, self.event_size, axis=0)

    @cached_property
    def total(self) -> np.ndarray:
        return self._total_coverage(self.numreads)

    @cached_property
    def bootstrap_total(self) -> np.ndarray:
        return self._total_coverage(self.bootstraps)

    @cached_property
    def other(self) -> np.ndarray:
        return self.total - self.numreads

    @cached_property
    def bootstrap_other(self) -> np.ndarray:
        return self.bootstrap_total - self.bootstraps

    @cached_property
    def alpha(self) -> np.ndarray:
        return self.numreads + self.alpha_prior

    @cached_property
    def beta(self) -> np.ndarray:
        return self.other + self.alpha_prior

    @cached_property
    def bootstrap_alpha(self) -> np.ndarray:
        return self.bootstraps + self.alpha_prior[:, np.newaxis]

    @cached_property
    def bootstrap_beta(self) -> np.ndarray:
        return self.bootstrap_other + self.beta_prior[:, np.newaxis]

    @staticmethod
    def _beta_mean(alpha, beta):
        return alpha / (alpha + beta)

    @staticmethod
    def _beta_var(mean, alpha, beta):
        return mean * mean / (1 + alpha + beta)

    @cached_property
    def posterior_mean(self) -> np.ndarray:
        return self._beta_mean(self.alpha, self.beta)

    @cached_property
    def _bootstrap_moments(self) -> Tuple[np.ndarray, np.ndarray]:
        return bm.moments(self.bootstrap_alpha, self.bootstrap_beta)

    @cached_property
    def bootstrap_posterior_mean(self) -> np.ndarray:
        return self._bootstrap_moments[0]

    @cached_property
    def posterior_variance(self) -> np.ndarray:
        return self._beta_var(self.posterior_mean, self.alpha, self.beta)

    @cached_property
    def bootstrap_posterior_variance(self) -> np.ndarray:
        return self._bootstrap_moments[1]

    def bootstrap_quantile(
        self,
        quantiles: Sequence[float] = [0.1, 0.9],
        nthreads: int = constants.DEFAULT_QUANTIFY_NTHREADS,
    ) -> np.ndarray:
        alpha = self.bootstrap_alpha
        beta = self.bootstrap_beta
        quantiles_arr = np.array(quantiles, dtype=alpha.dtype)
        result = np.empty((alpha.shape[0], len(quantiles)), dtype=alpha.dtype)

        def compute_slice(idx: slice) -> None:
            with np.errstate(divide="ignore"):
                bm.quantile(
                    quantiles_arr[np.newaxis],
                    alpha[idx, np.newaxis],
                    beta[idx, np.newaxis],
                    out=result[idx],
                )
            return

        from multiprocessing.dummy import Pool  # dummy is multithreading
        p = Pool(nthreads)
        WORKSIZE = 20000 // len(quantiles)
        p.map(
            compute_slice,
            [slice(x, x + WORKSIZE) for x in range(0, len(alpha), WORKSIZE)],
        )
        return result

    def bootstrap_discretized_pmf(
        self,
        nbins: int = 40,
        nthreads: int = constants.DEFAULT_QUANTIFY_NTHREADS,
    ) -> np.ndarray:
        if nbins < 2:
            raise ValueError(f"{nbins = } is invalid/trivial for discrete pmf")
        alpha = self.bootstrap_alpha
        beta = self.bootstrap_beta
        innerpoints = np.linspace(0, 1, 1 + nbins, dtype=alpha.dtype)[1:-1]
        result_cdf = np.empty((alpha.shape[0], 1 + nbins), dtype=alpha.dtype)
        result_cdf[:, 0] = 0
        result_cdf[:, -1] = 1
        result_cdf_inner = result_cdf[:, 1:-1]

        def compute_slice(idx: slice) -> None:
            with np.errstate(divide="ignore"):
                bm.cdf(
                    innerpoints[np.newaxis],
                    alpha[idx, np.newaxis],
                    beta[idx, np.newaxis],
                    out=result_cdf_inner[idx],
                )
            return

        from multiprocessing.dummy import Pool  # dummy is multithreading
        p = Pool(nthreads)
        WORKSIZE = 200000 // nbins
        p.map(
            compute_slice,
            [slice(x, x + WORKSIZE) for x in range(0, len(alpha), WORKSIZE)],
        )
        return np.diff(result_cdf, axis=1)

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
                self.offsets,
                self.events.is_intron,
                self.events.connection_idx,
            )
        )

    @staticmethod
    def _exons_formatted(
        exon_start: Sequence[int],
        exon_end: Sequence[int]
    ) -> List[str]:
        def format_coord(x):
            return x if x >= 0 else "na"

        return [
            f"{format_coord(a)}-{format_coord(b)}"
            for a, b in zip(exon_start, exon_end)
        ]

    def as_dataframe(self, sg: SpliceGraph) -> pd.DataFrame:
        # event information
        q_events = self.get_events(sg.introns, sg.junctions)
        event_id = sg.event_id(q_events.ref_exon_idx, q_events.event_type)
        event_description = sg.event_description(q_events.ref_exon_idx, q_events.event_type)
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
        return pd.DataFrame(
            {
                "seqid": seqid,
                "gene_id": gene_id,
                "ref_exon": np.repeat(
                    self._exons_formatted(ref_exon_start, ref_exon_end),
                    self.event_size
                ),
                "event_type": np.repeat([x.decode() for x in q_events.event_type], self.event_size),
                "is_intron": q_events.is_intron,
                "start": q_events.connection_start(),
                "end": q_events.connection_end(),
                "denovo": q_events.connection_denovo(),
                "psi_mean": self.posterior_mean,
                "bootstrap_psi_mean": self.bootstrap_posterior_mean,
                "bootstrap_psi_std": np.sqrt(self.bootstrap_posterior_variance),
                "total_reads": self.total,
                "gene_name": gene_name,
                "strand": strand,
                "other_exon": self._exons_formatted(
                    other_exon_start,
                    other_exon_end
                ),
                "event_id": np.repeat(event_id, self.event_size),
                "event_description": np.repeat(event_description, self.event_size),
            },
        ).set_index(["seqid", "gene_id", "ref_exon", "event_type", "is_intron", "start", "end"])

    def to_netcdf(
        self,
        path: Union[str, Path],
        pmf_bins: Optional[int] = 40,
        quantiles: Optional[Sequence[float]] = None,
        nthreads: int = 1,
    ) -> None:
        """Save to file with information so could be reopened as Events object"""
        if Path(path).exists():
            raise ValueError(f"Output {path} already exists")
        # save events
        (
            self.events
            .assign_coords(_offsets=("e_offsets_idx", self.offsets))
            .to_netcdf(path, "w", group=constants.NC_EVENTS)
        )
        # save quantifications
        df = xr.Dataset(
            {
                "psi_mean": ("ec_idx", self.posterior_mean),
                "bootstrap_psi_mean": ("ec_idx", self.bootstrap_posterior_mean),
                "bootstrap_psi_variance": ("ec_idx", self.bootstrap_posterior_variance),
            },
        )
        if pmf_bins:
            df = df.assign(
                psi_pmf=(
                    ("ec_idx", "psi_pmf_bin"),
                    self.bootstrap_discretized_pmf(pmf_bins, nthreads=nthreads),
                ),
            )
        if quantiles:
            df = df.assign(
                psi_quantiles=(
                    ("ec_idx", "quantile"),
                    self.bootstrap_quantile(quantiles, nthreads=nthreads),
                ),
            ).assign_coords(quantile=quantiles)
        df.to_netcdf(path, mode="a", group=constants.NC_EVENTSQUANTIFIED)
        return

    @classmethod
    def from_quantifier_group(
        cls,
        experiments: Sequence[Union[Path, str]],
        quantifiable: QuantifiableEvents,
    ) -> "QuantifiableCoverage":
        """Aggregate coverage over input experiments for quantifiable events"""
        if len(experiments) == 0:
            raise ValueError("No experiments passed into QuantifiableCoverage")
        checksums: Final[ConnectionsChecksum] = quantifiable.checksums
        quantifiable_offsets: Final[np.ndarray] = quantifiable.quantifiable_offsets
        coverage: xr.Dataset
        # get coverage first, checking checksums each time
        for x in experiments:
            with xr.open_dataset(x, group=constants.NC_EVENTS) as df:
                if checksums != ConnectionsChecksum(
                    introns=df.intron_hash, junctions=df.junction_hash
                ):
                    raise ValueError(
                        f"{x} has different intron/junction checksums"
                        f" than provided quantiable events"
                    )
            with xr.open_dataset(x, group=constants.NC_EVENTSCOVERAGE) as df:
                x_coverage = (
                    df[["numreads", "bootstraps"]]
                    .load()
                    .isel(ec_idx=quantifiable.event_connection_passed_mask)
                )
                try:
                    # this will fail if num-bootstraps is variable
                    # could handle but not for now
                    coverage += x_coverage
                except NameError:
                    coverage = x_coverage
        # get events
        events: xr.Dataset
        with xr.open_dataset(experiments[0], group=constants.NC_EVENTS) as df:
            events = (
                df.drop_dims("e_offsets_idx")
                .load()
                .isel(
                    e_idx=quantifiable.event_passed,
                    ec_idx=quantifiable.event_connection_passed_mask,
                )
            )
        return QuantifiableCoverage(
            quantifiable_offsets,
            coverage.numreads.values,
            coverage.bootstraps.values,
            events,
        )
