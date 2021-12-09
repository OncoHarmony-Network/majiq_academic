"""
PsiControlsSummary.py

Summary of independent experiments/prefixes quantifications from PsiCoverage
for use as controls

Author: Joseph K Aicher
"""

from functools import cached_property
from pathlib import Path
from typing import Final, Hashable, List, Sequence, Set, Union, cast

import xarray as xr
from dask.delayed import Delayed
from dask.distributed import progress

import new_majiq.constants as constants
from new_majiq._workarounds import _load_zerodim_variables
from new_majiq.Events import Events, _Events
from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.PsiCoverage import PsiCoverage, min_experiments


def _q_from_alpha(alpha: xr.DataArray) -> xr.DataArray:
    """Two-sided quantiles from outer product of alpha, is_lb

    Two-sided quantiles from outer product of alpha, is_lb.
    When is_lb == True, q = 0.5 * alpha. Otherwise, q = 1 - 0.5 * alpha.
    """
    is_lb = xr.DataArray(_ := [False, True], [("is_lb", _)], name="is_lb")
    q_lb = 0.5 * alpha  # lower bound quantiles
    q_ub = 1 - q_lb  # upper bound quantiles
    return q_lb.where(is_lb, q_ub).rename("q")


def _psirange_from_psiquantiles(quantiles: xr.DataArray) -> xr.DataArray:
    """Range between lower/upper quantile (quantiles[..., is_lb])"""
    return (quantiles.sel(is_lb=False) - quantiles.sel(is_lb=True)).rename("psi_range")


class PsiControlsSummary(object):
    """Summary of PSI posterior means over large group of controls

    Summary of PSI posterior means over large group of controls noting the
    number of experiments which passed and quantiles over the experiments which
    passed

    Parameters
    ----------
    df: xr.Dataset
        Dimensions:
            controls_alpha[controls_alpha]
                threshold used for quantiles in both directions (motivated by
                two-sided tests interpretation on CDFs for null-hypothesis
                testing). Quantiles are 0.5 * alpha, 1 - 0.5 * alpha
            is_lb[bool] = [False, True]
                Should calculated quantile be (0.5 * alpha) (True) or
                (1 - 0.5 * alpha) (False)
        Data variables:
            num_passed[ec_idx]
                Number of prefixes that passed quantification thresholds
            psi_median[ec_idx]
                Median of raw_psi_mean over prefixes
            psi_quantile[ec_idx, controls_alpha, is_lb]
                Quantiles of raw_psi_mean over prefixes
        Attributes:
            prefixes: List[str]
                prefixes that were summarized over
    events: xr.Dataset
        dataset that can be loaded along with matching introns/junctions as
        Events

    See Also
    --------
    PsiControlsSummary.from_psicov
    PsiControlsSummary.from_zarr
    """

    def __init__(self, df: xr.Dataset, events: xr.Dataset):
        REQUIRED_INDEXES: Set[Hashable] = {"controls_alpha", "is_lb"}
        REQUIRED_VARS: Set[Hashable] = {"num_passed", "psi_median", "psi_quantile"}
        if missing := REQUIRED_INDEXES - df.indexes.keys():
            raise ValueError(f"df missing requires indexed dimensions ({missing := })")
        if missing := REQUIRED_VARS - df.variables.keys():
            raise ValueError(f"df missing required data variables ({missing := })")
        if missing := {v for v in REQUIRED_VARS if "ec_idx" not in df[v].dims}:
            raise ValueError(
                f"df variables missing required ec_idx dimension ({missing = })"
            )
        if missing := REQUIRED_INDEXES - set(df["psi_quantile"].dims):
            raise ValueError(
                f"df.psi_quantile missing required dimensions ({missing = })"
            )
        if "q" not in df.variables:
            df = df.assign_coords(controls_q=_q_from_alpha(df["controls_alpha"]))
        if "prefixes" not in df.attrs or not isinstance(df.attrs["prefixes"], list):
            raise ValueError("df.attrs missing required list attribute 'prefixes'")
        if df.sizes["ec_idx"] != events.sizes["ec_idx"]:
            raise ValueError("df and events do not have same number of junctions")
        self.df: Final[xr.Dataset] = df
        self.events: Final[xr.Dataset] = events
        return

    @property
    def controls_alpha(self) -> xr.DataArray:
        return self.df["controls_alpha"]

    @property
    def alpha(self) -> xr.DataArray:
        """alias for controls_alpha"""
        return self.controls_alpha

    @property
    def is_lb(self) -> xr.DataArray:
        return self.df["is_lb"]

    @property
    def controls_q(self) -> xr.DataArray:
        return self.df["controls_q"]

    @property
    def q(self) -> xr.DataArray:
        """alias for controls_q"""
        return self.controls_q

    @property
    def num_passed(self) -> xr.DataArray:
        return self.df["num_passed"]

    @property
    def prefixes(self) -> List[str]:
        return self.df.attrs["prefixes"]

    @property
    def num_prefixes(self) -> int:
        return len(self.prefixes)

    @property
    def num_connections(self) -> int:
        return self.df.sizes["ec_idx"]

    def passed_min_experiments(
        self,
        min_experiments_f: Union[
            float, Sequence[float]
        ] = constants.DEFAULT_OUTLIERS_MINEXPERIMENTS,
    ) -> xr.DataArray:
        """Get boolean mask of events that pass enough experiments"""
        if isinstance(min_experiments_f, float):
            min_experiments_f = [min_experiments_f]
        self_min_experiments = xr.DataArray(
            [min_experiments(x, self.num_prefixes) for x in min_experiments_f],
            [("min_experiments_f", min_experiments_f)],
            name="min_experiments",
        )
        return self.num_passed >= self_min_experiments

    @property
    def psi_median(self) -> xr.DataArray:
        return self.df["psi_median"]

    @property
    def psi_quantile(self) -> xr.DataArray:
        return self.df["psi_quantile"]

    @cached_property
    def psi_range(self) -> xr.DataArray:
        """For each controls_alpha, range between lower/upper quantiles (scale)"""
        return _psirange_from_psiquantiles(self.psi_quantile)

    def load(self, show_progress: bool = False) -> "PsiControlsSummary":
        """Ensure that all data is computed/loaded in memory"""
        if show_progress:
            self.df.persist()
            progress(*(x.data for x in self.df.variables.values() if x.chunks))
        self.df.load()
        return self

    @classmethod
    def from_psicov(
        cls,
        psicov: PsiCoverage,
        alpha: Union[float, Sequence[float]] = constants.DEFAULT_OUTLIERS_ALPHA,
    ) -> "PsiControlsSummary":
        if isinstance(alpha, float):
            alpha = [alpha]
        else:
            # make sure alpha is sorted and has unique elements
            alpha = sorted(set(alpha))
        # get everything for df but the quantiles
        df = xr.Dataset(
            dict(
                psi_median=psicov.raw_psi_mean_population_median,
            ),
            dict(
                num_passed=psicov.num_passed,
                controls_alpha=("controls_alpha", alpha),
            ),
            dict(prefixes=psicov.prefixes),
        ).assign_coords(
            controls_q=lambda df: _q_from_alpha(df["controls_alpha"]),
        )
        stacked_q = df["controls_q"].stack(_idx=["controls_alpha", "is_lb"])
        df = df.assign(
            psi_quantile=(
                # get quantiles on stacked_q (flattened)
                psicov.raw_psi_mean_population_quantile(
                    stacked_q.values.tolist(), "controls_q"
                )
                # make controls_q indexed by dimension _idx, set its values
                .swap_dims(controls_q="_idx").assign_coords(_idx=stacked_q["_idx"])
                # unstack multiindex we created, get original dimensions back
                .unstack("_idx")
            )
        ).drop_vars(["event_size", "lsv_idx"], errors="ignore")
        return PsiControlsSummary(df, psicov.events)

    @classmethod
    def from_zarr(cls, path: Union[str, Path]) -> "PsiControlsSummary":
        return PsiControlsSummary(
            xr.open_zarr(path, group=constants.NC_PSICONTROLS),
            xr.open_zarr(path, group=constants.NC_EVENTS),
        )

    def to_zarr(
        self,
        path: Union[str, Path],
        consolidated: bool = True,
        show_progress: bool = False,
    ) -> None:
        """Save PSI coverage dataset as zarr

        Parameters
        ----------
        path: Union[str, Path]
            Path for output file in zarr format
        consolidated: bool
            When saving the file make sure that it is consolidated. In general,
            if you are appending a bunch of files together, it can make sense
            to set consolidated=False, and consolidate on the last write (only
            consolidate once). But, don't forget to consolidate at the end.
        show_progress: bool
            Attempt to show progress on distributed cluster for Dask
        """
        # save df
        save_df_future = cast(
            Delayed,
            self.df.chunk(self.df.sizes)
            .pipe(_load_zerodim_variables)
            .to_zarr(
                path,
                mode="w",
                group=constants.NC_PSICONTROLS,
                consolidated=False,
                compute=False,
            ),
        )
        if show_progress:
            save_df_future = save_df_future.persist()
            progress(save_df_future)
        else:
            save_df_future.compute()
        # save events
        self.events.chunk(self.events.sizes).to_zarr(
            path, mode="a", group=constants.NC_EVENTS, consolidated=consolidated
        )
        return

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
            f"PsiControlsSummary[{self.num_connections}]"
            f" for {self.num_prefixes} experiments [{print_prefixes}]"
        )
