"""
DeltaPsiDataset.py

Pre-computed quantifications for DeltaPsi for use in VOILA visualizations

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import Final, List, Optional, cast

import numpy as np
import pandas as pd
import xarray as xr
from dask.distributed import progress

import new_majiq.constants as constants

from ..core.DeltaPsi import DeltaPsi
from ..core.DPsiPrior import DPsiPrior
from ..core.Events import Events, _Events
from ..core.GeneIntrons import GeneIntrons
from ..core.GeneJunctions import GeneJunctions
from ..core.PMFSummaries import PMFSummaries
from ..core.PsiCoverage import PsiCoverage
from ..core.SpliceGraph import SpliceGraph


class DeltaPsiDataset(object):
    """Precomputed DeltaPsi quantifications for use in VOILA visualizations"""

    EXPECTED_VARIABLES: Final = {
        # identity of groups
        "grp": ("grp",),
        "comparison_grp1": ("comparison",),
        "comparison_grp2": ("comparison",),
        # identity of bins
        "pmf_bin_start": ("pmf_bin",),
        "pmf_bin_end": ("pmf_bin",),
        # discretized posterior on deltapsi
        "discrete_bootstrap_logposterior": ("comparison", "ec_idx", "pmf_bin"),
        # prior on deltapsi
        "prior_a": ("comparison", "mixture_component"),
        "prior_pmix": ("comparison", "mixture_component"),
        # whether a comparison was passed
        "passed": ("comparison", "ec_idx"),
        # values for PSI
        "raw_psi_mean": ("grp", "ec_idx"),
        "raw_psi_std": ("grp", "ec_idx"),
        "bootstrap_psi_std": ("grp", "ec_idx"),
        "raw_coverage": ("grp", "ec_idx"),
        "approximate_alpha": ("grp", "ec_idx"),
        "approximate_beta": ("grp", "ec_idx"),
    }

    @classmethod
    def from_deltapsi(
        cls,
        dpsi: DeltaPsi,
        compute: bool = False,
        show_progress: bool = False,
    ) -> "DeltaPsiDataset":
        """Create :class:`DeltaPsiDataset` using :class:`DeltaPsi`

        Create :class:`DeltaPsiDataset` with minimal data required for
        visualization of DeltaPsi analysis in VOILA

        Parameters
        ----------
        dpsi: DeltaPsi
            Handle to coverage with full availability of possible DeltaPsi
            related computations
        compute: bool
            If true, eagerly perform inference of deltapsi posteriors
        show_progress: bool
            If `compute` is true, show progress bar in Dask
        """
        dpsi_ds = (
            xr.Dataset(
                {
                    "discrete_bootstrap_logposterior": dpsi.discrete_bootstrap_logposterior,
                    "prior_a": dpsi.prior.a,
                    "prior_pmix": dpsi.prior.pmix,
                    "passed": dpsi.passed,
                },
            )
            .expand_dims(comparison=1)
            .assign_coords(
                comparison_grp1=("comparison", [dpsi.name1]),
                comparison_grp2=("comparison", [dpsi.name2]),
            )
        )
        PSI_VARIABLES = [
            "raw_psi_mean",
            "raw_psi_std",
            "bootstrap_psi_std",
            "raw_coverage",
            "approximate_alpha",
            "approximate_beta",
        ]
        psi_ds = (
            xr.concat(
                [
                    dpsi.psi1.dataset(PSI_VARIABLES).drop_vars("any_passed"),
                    dpsi.psi2.dataset(PSI_VARIABLES).drop_vars("any_passed"),
                ],
                dim="prefix",
            )
            .rename_vars(prefix="grp")
            .swap_dims(prefix="grp")
        )
        df = xr.merge([dpsi_ds, psi_ds], compat="override", join="exact")
        if compute:
            if show_progress:
                df = df.persist()
                progress(*(x.data for x in df.variables.values() if x.chunks))
            df = df.load()
        return DeltaPsiDataset(df, dpsi.psi1.events)

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
        """Initialize :class:`DeltaPsiDataset` with specified underlying data

        Parameters
        ----------
        df: xr.Dataset
            Variables/coordinates matching DeltaPsiDataset.EXPECTED_VARIABLES
        events: xr.Dataset
            dataset that can be loaded along with matching introns/junctions as
            Events
        """
        # verify that every variable expected is present
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in DeltaPsiDataset df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(
                    f"DeltaPsiDataset df['{var}' must have dimensions {var_dims}"
                )
        # verify that all compared groups have PSI information available
        if not (
            df["comparison_grp1"].load().isin(df["grp"]).all()
            and df["comparison_grp2"].load().isin(df["grp"]).all()
        ):
            raise ValueError(
                "Not all compared groups"
                f" ({df['comparison_grp1'] = }, {df['comparison_grp2'] = })"
                f" have PSI information on dimension {df['grp'] = }"
            )
        # make these values available
        self.df: Final[xr.Dataset] = df.transpose(
            ..., "ec_idx", "mixture_component", "pmf_bin"
        )
        self.events: Final[xr.Dataset] = events
        return

    def __repr__(self) -> str:
        return (
            f"DeltaPsiDataset[{self.num_connections}]"
            f" for {self.num_comparisons} comparisons of {self.num_groups} groups"
        )

    @property
    def num_connections(self) -> int:
        """Number of event connections"""
        return self.df.sizes["ec_idx"]

    @property
    def num_comparisons(self) -> int:
        """Number of comparisons in this dataset"""
        return self.df.sizes["comparison"]

    @property
    def num_groups(self) -> int:
        """Number of groups in this dataset"""
        return self.df.sizes["grp"]

    @property
    def comparisons(self):
        """Enumerate which groups were compared in this dataset"""
        return tuple(
            zip(self.df["comparison_grp1"].values, self.df["comparison_grp2"].values)
        )

    @property
    def passed(self) -> xr.DataArray:
        """Boolean mask for passed events for each comparison"""
        return self.df["passed"]

    @property
    def raw_psi_mean(self) -> xr.DataArray:
        """PSI posterior means from raw coverage per group"""
        return self.df["raw_psi_mean"]

    @property
    def raw_psi_std(self) -> xr.DataArray:
        """PSI posterior standard deviations from raw coverage per group"""
        return self.df["raw_psi_std"]

    @property
    def bootstrap_psi_std(self) -> xr.DataArray:
        """PSI posterior standard deviations from bootstrap coverage per group"""
        return self.df["bootstrap_psi_std"]

    @property
    def raw_coverage(self) -> xr.DataArray:
        """Raw coverage for each connection"""
        return self.df["raw_coverage"]

    def psi_approximate_discretized_pmf(
        self, nbins: int = constants.DEFAULT_QUANTIFY_PSIBINS
    ) -> xr.DataArray:
        """Compute discretized PMF of approximate/smoothed bootstrap posterior

        Parameters
        ----------
        nbins: int
            Number of uniform bins on [0, 1] on which probability mass will be
            computed
        """
        return PsiCoverage._compute_posterior_discretized_pmf(
            self.df["approximate_alpha"], self.df["approximate_beta"], nbins=nbins
        )

    @cached_property
    def dpsi_prior(self) -> DPsiPrior:
        """Prior over deltapsi for each comparison"""
        return DPsiPrior(self.df["prior_a"], self.df["prior_pmix"])

    @cached_property
    def dpsi_posterior(self) -> PMFSummaries:
        """:class:`PMFSummaries` for average bootstrapped dpsi posteriors"""
        return PMFSummaries(
            cast(xr.DataArray, np.exp(self.df["discrete_bootstrap_logposterior"]))
        )

    def probability_changing(
        self, changing_threshold: float = constants.DEFAULT_DPSI_CHANGING_THRESHOLD
    ):
        return 1 - self.dpsi_posterior.interval_probability(
            -changing_threshold, changing_threshold
        )

    def probability_nonchanging(
        self,
        nonchanging_threshold: float = constants.DEFAULT_DPSI_NONCHANGING_THRESHOLD,
    ):
        return self.dpsi_posterior.interval_probability(
            -nonchanging_threshold, nonchanging_threshold
        )

    def get_events(
        self,
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> Events:
        """Construct :py:class:`Events` using saved dataset and introns, junctions

        Parameters
        ----------
        introns: GeneIntrons
        junctions: GeneJunctions

        Returns
        -------
        Events
        """
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

    def to_dataframe(
        self,
        sg: Optional[SpliceGraph] = None,
        changing_threshold: float = constants.DEFAULT_DPSI_CHANGING_THRESHOLD,
        nonchanging_threshold: float = constants.DEFAULT_DPSI_NONCHANGING_THRESHOLD,
        show_progress: bool = False,
    ) -> pd.DataFrame:
        """Return table of quantifications for TSV output

        Parameters
        ----------
        sg: Optional[SpliceGraph]
            If provided, splicegraph with introns/junctions consistent with
            events used to annotate resulting dataframe
        changing_threshold: float
            threshold t for P(abs(dPSI) >= t), the posterior probability that
            dPSI is changing by more than this amount
        nonchanging_threshold: float
            threshold t for P(abs(dPSI) <= t), the posterior probability that
            dPSI is changing by less than this amount
        show_progress: bool
            show progress bar in dask if enabled
        """
        # build tables with columns that will be concatenated together
        concat_df: List[pd.DataFrame] = list()
        # add dataframe with events annotations
        if sg is not None:
            concat_df.append(self.get_events(sg.introns, sg.junctions).ec_dataframe)
        # build dataset with quantifcations to load simultaneously
        ds = xr.Dataset(
            {
                "any_passed": self.passed.any("comparison"),
                "dpsi_mean": self.dpsi_posterior.mean,
                "dpsi_std": self.dpsi_posterior.standard_deviation,
                "probability_changing": 1
                - self.dpsi_posterior.interval_probability(
                    -changing_threshold, changing_threshold
                ),
                "probability_nonchanging": self.dpsi_posterior.interval_probability(
                    -nonchanging_threshold, nonchanging_threshold
                ),
                "raw_psi_mean": self.raw_psi_mean,
                "raw_psi_std": self.raw_psi_std,
                "bootstrap_psi_std": self.bootstrap_psi_std,
                "raw_coverage": self.raw_coverage,
            },
        )
        # load/compute into memory
        if show_progress:
            ds = ds.persist()
            progress(*(x.data for x in ds.variables.values() if x.chunks))
        ds = ds.load()
        # add deltapsi into table
        for i, (grp1, grp2) in enumerate(self.comparisons):
            comparison_prefix = ""
            if self.num_comparisons > 1:
                comparison_prefix = f"{grp1}-vs-{grp2}_"
            concat_df.append(
                ds[[name for name, v in ds.items() if "comparison" in v.dims]]
                .isel(comparison=i, drop=True)
                .to_dataframe()
                .add_prefix(comparison_prefix)
            )
        # add psi into table
        for i, grp in enumerate(self.df["grp"].values):
            concat_df.append(
                ds[[name for name, v in ds.items() if "grp" in v.dims]]
                .isel(grp=i, drop=True)
                .to_dataframe()
                .add_prefix(f"{grp}_")
            )
        return pd.concat(concat_df, axis=1, join="inner").loc[ds["any_passed"].values]
