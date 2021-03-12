"""
Quantifier.py

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from functools import cached_property
from typing import (
    Final,
    NamedTuple,
    Sequence,
    Tuple,
    Union,
)
from pathlib import Path


class QuantifierThresholds(NamedTuple):
    minreads: float = 10
    minbins: float = 3
    min_experiments_f: float = 0.5


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
        return np.reciprocal(self.connection_event_size, dtype=np.float32)

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
    def bootstrap_posterior_mean(self) -> np.ndarray:
        return self._beta_mean(self.bootstrap_alpha, self.bootstrap_beta)

    @cached_property
    def posterior_variance(self) -> np.ndarray:
        return self._beta_var(self.posterior_mean, self.alpha, self.beta)

    @cached_property
    def bootstrap_posterior_variance(self) -> np.ndarray:
        return self._beta_var(
            self.bootstrap_posterior_mean,
            self.bootstrap_alpha,
            self.bootstrap_beta,
        )

    @cached_property
    def bootstrap_mean(self) -> np.ndarray:
        return np.array(self.bootstrap_posterior_mean.mean(axis=1))

    @cached_property
    def bootstrap_variance(self) -> np.ndarray:
        return np.array(
            # law of total variance
            self.bootstrap_posterior_mean.var(axis=1)
            + self.bootstrap_posterior_variance.mean(axis=1)
        )

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
