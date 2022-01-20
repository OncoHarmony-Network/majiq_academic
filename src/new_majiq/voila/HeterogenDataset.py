"""
HeterogenDataset.py

Intermediate quantifications for HET for use in VOILA visualizations

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import Collection, Dict, Final, List, Optional, Sequence, Union, cast

import numpy as np
import pandas as pd
import xarray as xr
from dask.delayed import Delayed
from dask.distributed import progress

import new_majiq.constants as constants

from ..core.Events import Events, _Events
from ..core.GeneIntrons import GeneIntrons
from ..core.GeneJunctions import GeneJunctions
from ..core.Heterogen import Heterogen
from ..core.MixinPsiInference import MixinRawPsiMeanPopulation
from ..core.SpliceGraph import SpliceGraph


class HeterogenGroupPsi(MixinRawPsiMeanPopulation):
    """Encapsulate statistics on raw_psi_mean for Heterogen groups

    Parameters
    ----------
    raw_psi_mean: xr.DataArray
        array[..., prefix] with raw_psi_mean for each experiment/prefix in the
        group
    """

    def __init__(self, raw_psi_mean: xr.DataArray):
        if "prefix" not in raw_psi_mean.dims:
            raise ValueError("raw_psi_mean does not have required dimension 'prefix'")
        if "ec_idx" not in raw_psi_mean.dims:
            raise ValueError("raw_psi_mean does not have required dimension 'ec_idx'")
        self._raw_psi_mean: Final[xr.DataArray] = raw_psi_mean
        return

    def __repr__(self) -> str:
        return (
            f"HeterogenGroupPsi[{self.num_connections}]"
            f" for {self.num_prefixes} prefixes"
        )

    @property
    def raw_psi_mean(self) -> xr.DataArray:
        return self._raw_psi_mean

    @cached_property
    def num_passed(self) -> xr.DataArray:
        """Number of prefixes for which an event was passed"""
        return self.raw_psi_mean.count("prefix")

    @property
    def num_connections(self) -> int:
        """Number of event connections"""
        return self.raw_psi_mean.sizes["ec_idx"]

    @property
    def num_prefixes(self) -> int:
        return self.raw_psi_mean.sizes["prefix"]

    @property
    def prefixes(self) -> List[str]:
        return self.raw_psi_mean["prefix"].values.tolist()


class HeterogenDataset(object):
    """Precomputed Heterogen quantifications for use in VOILA visualizations"""

    EXPECTED_VARIABLES: Final = {
        # identity of comparisons
        "comparison_grp1": ("comparison",),
        "comparison_grp2": ("comparison",),
        # identity of groups
        "grp": ("grp",),
        "grp_size": ("grp",),  # check vs prefix_grp
        # identity of prefixes
        "prefix": ("prefix",),
        "prefix_grp": ("prefix",),
        # identity of stats
        "stats": ("stats",),
        # bins for groups
        "pmf_bin_start": ("pmf_bin",),
        "pmf_bin_end": ("pmf_bin",),
        "psi_pmf": ("grp", "ec_idx", "pmf_bin"),
        # psi for prefixes
        "raw_psi_mean": ("ec_idx", "prefix"),
        # psisamples for comparisons
        "psisamples": ("comparison",),
        "pval_quantile": ("pval_quantile",),
        # stats
        "raw_pvalue": ("comparison", "ec_idx", "stats"),
        "approximate_pvalue": ("comparison", "ec_idx", "stats"),
        "approximate_pvalue_quantiles": (
            "comparison",
            "ec_idx",
            "stats",
            "pval_quantile",
        ),
    }

    @classmethod
    def from_heterogen(
        cls,
        het: Heterogen,
        pvalue_quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        psibins: int = constants.DEFAULT_QUANTIFY_PSIBINS,
    ) -> "HeterogenDataset":
        # define stats vector to input
        if isinstance(use_stats, str):
            use_stats = [use_stats]  # make sure it's always a collection
        else:
            use_stats = sorted(set(use_stats))
        if any(unknown := [x for x in use_stats if x not in constants.STATS_AVAILABLE]):
            raise ValueError(
                f"Input statistics {unknown} unavailable"
                f" (available {set(constants.STATS_AVAILABLE.keys())})"
            )
        # construct psi_pmf (TODO: make stochastic if too slow)
        psi_pmf = xr.concat(
            [
                het.psi1.approximate_discretized_pmf(psibins).mean("prefix"),
                het.psi2.approximate_discretized_pmf(psibins).mean("prefix"),
            ],
            dim="grp",
        )
        # compute stats
        raw_stats = het.raw_stats(use_stats=use_stats).expand_dims(comparison=1)
        approximate_stats = het.approximate_stats(
            quantiles=pvalue_quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
            drop=False,
        ).expand_dims(comparison=1)
        # construct raw_psi_mean
        raw_psi_mean = xr.concat(
            [het.psi1.raw_psi_mean, het.psi2.raw_psi_mean], dim="prefix"
        )
        # construct dataset
        df = xr.Dataset(
            {
                "psi_pmf": psi_pmf,
                "raw_psi_mean": raw_psi_mean,
                "raw_pvalue": raw_stats["pvalue"],
                "approximate_pvalue": approximate_stats["pvalue"],
                "approximate_pvalue_quantiles": approximate_stats["pvalue_quantiles"],
            },
            {
                "comparison_grp1": ("comparison", [het.name1]),
                "comparison_grp2": ("comparison", [het.name2]),
                "grp": [het.name1, het.name2],
                "grp_size": ("grp", [het.psi1.num_prefixes, het.psi2.num_prefixes]),
                # "prefix": het.psi1.prefixes + het.psi2.prefixes,
                "prefix_grp": (
                    "prefix",
                    [het.name1] * het.psi1.num_prefixes
                    + [het.name2] * het.psi2.num_prefixes,
                ),
                # "stats": use_stats,
                "psisamples": ("comparison", [psisamples]),
            },
        )
        return HeterogenDataset(df, het.psi1.events)

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
        """Initialize :class:`HeterogenDataset` with specified underlying data

        Parameters
        ----------
        df: xr.Dataset
            Variables/coordinates matching HeterogenDataset.EXPECTED_VARIABLES
        events: xr.Dataset
            dataset that can be loaded along with matching introns/junctions as
            Events
        """
        # verify that every variable expected is present
        for var, var_dims in self.EXPECTED_VARIABLES.items():
            if var not in df.variables:
                raise ValueError(f"{var} must be in df variables")
            if set(df[var].dims) != set(var_dims):
                raise ValueError(f"df['{var}'] must have dimensions {var_dims}")
        # verify that all compared groups are recognized
        if not (
            df["comparison_grp1"].load().isin(df["grp"]).all()
            and df["comparison_grp2"].load().isin(df["grp"]).all()
        ):
            raise ValueError(
                "Not all compared groups"
                f" ({df['comparison_grp1'] = }, {df['comparison_grp2'] = })"
                f" have group information on dimension {df['grp'] = }"
            )
        # verify that prefixes/groups match up appropriately
        if (
            not df["prefix_grp"]
            .to_series()
            .value_counts()
            .equals(df["grp_size"].to_series())
        ):
            raise ValueError("Mismatch between prefix_grp numbers and grp_size")
        # make these values available
        self.events: Final[xr.Dataset] = events
        self.df: Final[xr.Dataset] = df.transpose(
            "comparison", "grp", "ec_idx", "prefix", "stats", "pval_quantile", "pmf_bin"
        )
        self.groups: Final[Dict[str, HeterogenGroupPsi]] = {
            grp: HeterogenGroupPsi(
                self.df["raw_psi_mean"].sel(prefix=self.df["prefix_grp"].load() == grp)
            )
            for grp in self.df["grp"].values
        }
        return

    def to_zarr(
        self,
        path: Union[str, Path],
        ec_chunksize: int = constants.DEFAULT_COVERAGE_CHUNKS,
        consolidated: bool = True,
        show_progress: bool = False,
    ) -> None:
        """Save :class:`HeterogenDataset` to specified path"""
        save_df_future = cast(
            Delayed,
            self.df.to_zarr(
                path,
                mode="w",
                group=constants.NC_HETEROGEN,
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

    @classmethod
    def from_zarr(
        cls, path: Union[str, Path, List[Union[str, Path]]]
    ) -> "HeterogenDataset":
        """Load :class:`HeterogenDataset` from one or more specified paths"""
        if not isinstance(path, list):
            path = [path]
        df = xr.open_mfdataset(
            path,
            engine="zarr",
            group=constants.NC_HETEROGEN,
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
        return HeterogenDataset(df, events_df)

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
    def num_prefixes(self) -> int:
        return self.df.sizes["prefix"]

    @property
    def psibins(self) -> int:
        """Number of bins used to bin PSI for deltapsi computation"""
        return self.df.sizes["pmf_bin"]

    @property
    def comparisons(self):
        """Enumerate which groups were compared in this dataset"""
        return tuple(
            zip(self.df["comparison_grp1"].values, self.df["comparison_grp2"].values)
        )

    def __repr__(self) -> str:
        return (
            f"HeterogenDataset[{self.num_connections}]"
            f" for {self.num_comparisons} comparisons of"
            f" {self.num_groups} groups, {self.num_prefixes} total prefixes"
        )

    @property
    def psi_pmf(self) -> xr.DataArray:
        """Average of PSI approximate bootstrap posterior distributions per group"""
        return self.df["psi_pmf"]

    @property
    def psisamples(self) -> xr.DataArray:
        return self.df["psisamples"]

    @property
    def num_stats(self) -> int:
        return self.df.sizes["stats"]

    @property
    def stats(self) -> xr.DataArray:
        return self.df["stats"]

    @property
    def pval_quantile(self) -> xr.DataArray:
        return self.df["pval_quantile"]

    @property
    def raw_pvalue(self) -> xr.DataArray:
        return self.df["raw_pvalue"]

    @property
    def approximate_pvalue(self) -> xr.DataArray:
        return self.df["approximate_pvalue"]

    @property
    def approximate_pvalue_quantiles(self) -> xr.DataArray:
        return self.df["approximate_pvalue_quantiles"]

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

    @cached_property
    def any_passed(self) -> xr.DataArray:
        return self.df["raw_psi_mean"].notnull().any("prefix")

    def to_dataframe(
        self,
        sg: Optional[SpliceGraph] = None,
        population_quantiles: Sequence[
            float
        ] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        show_progress: bool = False,
    ) -> pd.DataFrame:
        """Return table of quantifications for TSV output

        Parameters
        ----------
        sg: Optional[SpliceGraph]
            If provided, splicegraph with introns/junctions consistent with
            events used to annotate resulting dataframe
        population_quantiles: Sequence[float]
            quantiles of raw_psi_mean to evaluate per group
        show_progress: bool
            show progress bar in dask if enabled
        """
        # build tables with columns that will be concatenated together
        concat_df: List[pd.DataFrame] = list()
        # population quantiles to compute (add median, unique up to 3 decimal places)
        population_quantiles = sorted({0.5, *np.round(population_quantiles, 3)})
        # add dataframe with events annotations
        if sg is not None:
            concat_df.append(self.get_events(sg.introns, sg.junctions).ec_dataframe)
        # build dataset with quantifications to load simultaneously
        ds: xr.Dataset = (
            self.df.drop_dims("pmf_bin")
            .drop_vars("raw_psi_mean")
            .assign(
                any_passed=self.any_passed,
                raw_psi_quantile=xr.concat(
                    [
                        x.raw_psi_mean_population_quantile(
                            population_quantiles
                        ).expand_dims(grp=[grp])
                        for grp, x in self.groups.items()
                    ],
                    dim="grp",
                ),
                num_passed=xr.concat(
                    [
                        x.num_passed.expand_dims(grp=[grp])
                        for grp, x in self.groups.items()
                    ],
                    dim="grp",
                ),
            )
        )
        # load/compute into memory
        if show_progress:
            ds = ds.persist()
            progress(*(x.data for x in ds.variables.values() if x.chunks))
        ds = ds.load()
        # add comparisons into table
        for i, (grp1, grp2) in enumerate(self.comparisons):
            comparison_prefix = f"{grp1}-vs-{grp2}-" if self.num_comparisons > 1 else ""
            for j, statistic in enumerate(self.stats.values):
                stat_prefix = f"{statistic}-" if self.num_stats > 1 else ""
                concat_df.append(
                    ds[["raw_pvalue", "approximate_pvalue"]]
                    .isel(comparison=i, stats=j, drop=True)
                    .to_dataframe()
                    .assign(
                        **{
                            f"approximate_pvalue_quantiles_{q:0.3f}": ds[
                                "approximate_pvalue_quantiles"
                            ]
                            .isel(comparison=i, stats=j, pval_quantile=k, drop=True)
                            .to_series()
                            for k, q in enumerate(self.pval_quantile.values)
                        }
                    )
                    .add_prefix(f"{comparison_prefix}{stat_prefix}")
                )
        # add raw_psi_quantile into table
        df_quantile = (
            ds["raw_psi_quantile"].to_series().unstack(["grp", "population_quantile"])
        )
        df_quantile.columns = [
            f"{grp}-raw_psi_quantile_{q:0.3f}" for grp, q in df_quantile.columns
        ]
        concat_df.append(df_quantile)
        del df_quantile

        return pd.concat(concat_df, axis=1, join="inner").loc[ds["any_passed"].values]
