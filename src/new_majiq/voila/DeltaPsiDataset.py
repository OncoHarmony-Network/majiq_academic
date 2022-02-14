"""
DeltaPsiDataset.py

Pre-computed quantifications for DeltaPsi for use in VOILA visualizations

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import Final, List, Optional, Union, cast

import numpy as np
import pandas as pd
import xarray as xr
from dask.delayed import Delayed
from dask.distributed import progress

import new_majiq.constants as constants

from ..core.DeltaPsi import DeltaPsi, DeltaPsiPMF
from ..core.DPsiPrior import DPsiPrior
from ..core.MixinHasEvents import MixinHasEvents
from ..core.MixinPsiInference import MixinApproximatePsi, MixinRawPsi
from ..core.SpliceGraph import SpliceGraph


class DeltaPsiGroupPsiStatistics(MixinRawPsi, MixinApproximatePsi):
    """Encapsulate group statistics around PSI for DeltaPsi groups

    Parameters
    ----------
    df: xr.Dataset
    """

    EXPECTED_VARIABLES: Final = {
        # values for PSI
        "raw_alpha": ("grp", "ec_idx"),
        "raw_beta": ("grp", "ec_idx"),
        "approximate_alpha": ("grp", "ec_idx"),
        "approximate_beta": ("grp", "ec_idx"),
        "raw_coverage": ("grp", "ec_idx"),
    }

    def __init__(self, df: xr.Dataset):
        # verify that every variable expected is present
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(f"df['{var}'] must have dimensions {var_dims}")
        self.df: Final[xr.Dataset] = df
        return

    @property
    def raw_coverage(self) -> xr.DataArray:
        """Raw coverage per event connection"""
        return self.df["raw_coverage"]

    @property
    def raw_alpha(self) -> xr.DataArray:
        """array(...) of alpha for raw posterior"""
        return self.df["raw_alpha"]

    @property
    def raw_beta(self) -> xr.DataArray:
        """array(...) of beta for raw posterior"""
        return self.df["raw_beta"]

    @property
    def approximate_alpha(self) -> xr.DataArray:
        """array(...) of alpha approximating bootstrap posterior"""
        return self.df["approximate_alpha"]

    @property
    def approximate_beta(self) -> xr.DataArray:
        """array(...) of beta approximating bootstrap posterior"""
        return self.df["approximate_beta"]


class DeltaPsiDataset(MixinHasEvents):
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
        "bootstrap_logposterior": ("comparison", "ec_idx", "pmf_bin"),
        # prior on deltapsi
        "prior_a": ("comparison", "mixture_component"),
        "prior_pmix": ("comparison", "mixture_component"),
        # whether a comparison was passed
        "passed": ("comparison", "ec_idx"),
        # values for PSI
        **DeltaPsiGroupPsiStatistics.EXPECTED_VARIABLES,
    }

    @classmethod
    def from_deltapsi(cls, dpsi: DeltaPsi) -> "DeltaPsiDataset":
        """Create :class:`DeltaPsiDataset` using :class:`DeltaPsi`

        Create :class:`DeltaPsiDataset` with minimal data required for
        visualization of DeltaPsi analysis in VOILA

        Parameters
        ----------
        dpsi: DeltaPsi
            Handle to coverage with full availability of possible DeltaPsi
            related computations
        """
        dpsi_ds = (
            xr.Dataset(
                {
                    "bootstrap_logposterior": dpsi.bootstrap_logposterior,
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
        psi_ds = (
            xr.concat(
                [
                    dpsi.psi1.dataset(
                        list(DeltaPsiGroupPsiStatistics.EXPECTED_VARIABLES.keys())
                    ).drop_vars("any_passed"),
                    dpsi.psi2.dataset(
                        list(DeltaPsiGroupPsiStatistics.EXPECTED_VARIABLES.keys())
                    ).drop_vars("any_passed"),
                ],
                dim="prefix",
            )
            .rename_vars(prefix="grp")
            .swap_dims(prefix="grp")
        )
        df = xr.merge([dpsi_ds, psi_ds], compat="override", join="exact")
        return DeltaPsiDataset(df, dpsi.psi1.events_df)

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
        # save events
        MixinHasEvents.__init__(self, df, events)
        # save df/groups after validating input
        # verify that every variable expected is present
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(f"df['{var}'] must have dimensions {var_dims}")
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
        self.groups: Final[DeltaPsiGroupPsiStatistics] = DeltaPsiGroupPsiStatistics(
            df.transpose(..., "ec_idx", "mixture_component", "pmf_bin")
        )
        return

    @property
    def df(self) -> xr.Dataset:
        """Underlying xarray dataset"""
        return self.groups.df

    @classmethod
    def from_zarr(
        cls, path: Union[str, Path, List[Union[str, Path]]]
    ) -> "DeltaPsiDataset":
        """Load :class:`DeltaPsiDataset` from one or more specified paths"""
        if not isinstance(path, list):
            path = [path]
        df = xr.open_mfdataset(
            path,
            engine="zarr",
            group=constants.NC_DELTAPSI,
            combine="nested",
            preprocess=lambda ds: ds.assign_coords(comparison=[ds.encoding["source"]]),
            join="outer",
            compat="no_conflicts",
            coords="minimal",
            data_vars="minimal",
        )
        if len(path) > 1:
            # attributes are defined by path[0]. We'd rather just have none
            df.attrs.clear()
        events_df = xr.open_zarr(path[0], group=constants.NC_EVENTS)
        return DeltaPsiDataset(df, events_df)

    def to_zarr(
        self,
        path: Union[str, Path],
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        consolidated: bool = True,
        show_progress: bool = False,
    ) -> None:
        """Save :class:`DeltaPsiDataset` to specified path"""
        save_df_future = cast(
            Delayed,
            self.df.to_zarr(
                path,
                mode="w",
                group=constants.NC_DELTAPSI,
                consolidated=False,
                compute=False,
            ),
        )
        if show_progress:
            save_df_future = save_df_future.persist()
            progress(save_df_future)
        else:
            save_df_future.compute()
        self.events_to_zarr(path, mode="a", consolidated=consolidated)
        return

    def __repr__(self) -> str:
        return (
            f"DeltaPsiDataset[{self.num_connections}]"
            f" for {self.num_comparisons} comparisons of {self.num_groups} groups"
        )

    @property
    def num_comparisons(self) -> int:
        """Number of comparisons in this dataset"""
        return self.df.sizes["comparison"]

    @property
    def num_groups(self) -> int:
        """Number of groups in this dataset"""
        return self.df.sizes["grp"]

    @property
    def psibins(self) -> int:
        """Number of bins used to bin PSI for deltapsi computation"""
        return self.df.sizes["pmf_bin"] // 2

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

    @cached_property
    def prior(self) -> DPsiPrior:
        """Prior over deltapsi for each comparison"""
        return DPsiPrior(self.df["prior_a"], self.df["prior_pmix"])

    @cached_property
    def discrete_logprior(self) -> xr.DataArray:
        return self.prior.discretized_logpmf(psibins=self.psibins)

    @property
    def bootstrap_logposterior(self) -> xr.DataArray:
        """Log-average bootstrap replicates after inference of deltapsi"""
        return self.df["bootstrap_logposterior"]

    @cached_property
    def discrete_prior(self) -> DeltaPsiPMF:
        """:class:`DeltaPsiPMF` for discretized deltapsi prior"""
        return DeltaPsiPMF(cast(xr.DataArray, np.exp(self.discrete_logprior)))

    @cached_property
    def bootstrap_posterior(self) -> DeltaPsiPMF:
        """:class:`DeltaPsiPMF` for average bootstrapped dpsi posteriors"""
        return DeltaPsiPMF(cast(xr.DataArray, np.exp(self.bootstrap_logposterior)))

    def to_dataframe(
        self,
        sg: Optional[SpliceGraph] = None,
        annotated: Optional[SpliceGraph] = None,
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
        annotated: Optional[SpliceGraph]
            If specified, override denovo definitions by comparing
            exons/introns/junctions to annotated splicegraph.
            Only used if `sg` also provided.
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
            concat_df.append(
                self.get_events(sg.introns, sg.junctions).ec_dataframe(
                    annotated=annotated
                )
            )
        # build dataset with quantifcations to load simultaneously
        ds = xr.Dataset(
            {
                "any_passed": self.passed.any("comparison"),
                "dpsi_mean": self.bootstrap_posterior.mean,
                "dpsi_std": self.bootstrap_posterior.standard_deviation,
                "probability_changing": self.bootstrap_posterior.probability_changing(
                    changing_threshold
                ),
                "probability_nonchanging": self.bootstrap_posterior.probability_nonchanging(
                    nonchanging_threshold
                ),
                "raw_psi_mean": self.groups.raw_psi_mean,
                "raw_psi_std": self.groups.raw_psi_std,
                "bootstrap_psi_std": self.groups.bootstrap_psi_std,
                "raw_coverage": self.groups.raw_coverage,
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
