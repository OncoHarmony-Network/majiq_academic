"""
Heterogen.py

Quantify/test differences between two groups of independent experiments

Author: Joseph K Aicher
"""

from functools import cached_property
from typing import Collection, Final, List, Sequence, Union

import numpy as np
import xarray as xr

import new_majiq.beta_mixture as bm
import new_majiq.constants as constants
from new_majiq.logger import get_logger
from new_majiq.PsiCoverage import PsiCoverage


class Heterogen(object):
    """Heterogeneous analysis between two groups of independent experiments

    Parameters
    ----------
    psi1, psi2: PsiCoverage
        PsiCoverage files to compare
    min_experiments_f
        Number or proportion of experiments required to pass in each group for
        quantification
    stats: Union[str, Collection[str]]
        Statistics to compute with (must be in constants.STATS_AVAILABLE)
    name1, name2: str
        Names to indicate group identity
    """

    def __init__(
        self,
        psi1: PsiCoverage,
        psi2: PsiCoverage,
        min_experiments_f: float = constants.DEFAULT_QUANTIFY_MINEXPERIMENTS,
        name1: str = "grp1",
        name2: str = "grp2",
    ):
        # normalize inputs, check for errors
        if psi1.num_connections != psi2.num_connections:
            raise ValueError(
                "psi1 and psi2 must have the same number of connections"
                f" ({psi1.num_connections=}, {psi2.num_connections=})"
            )
        if prefix_overlap := set(psi1.prefixes) & set(psi2.prefixes):
            # warning if there is overlap between prefixes
            get_logger().warning(
                f"Heterogen input groups have overlapping prefixes ({prefix_overlap})"
            )
        # get events that passed both psi1 and psi2 in enough experiments
        passed = psi1.passed_min_experiments(
            min_experiments_f
        ) & psi2.passed_min_experiments(min_experiments_f)
        self.psi1: Final[PsiCoverage] = psi1.mask_events(passed)
        self.psi2: Final[PsiCoverage] = psi2.mask_events(passed)
        self.passed: Final[xr.DataArray] = passed
        self.name1: Final[str] = name1
        self.name2: Final[str] = name2
        return

    @cached_property
    def labels(self) -> xr.DataArray:
        return xr.DataArray(
            np.repeat([True, False], [self.psi1.num_prefixes, self.psi2.num_prefixes]),
            dims=["prefix"],
        )

    @staticmethod
    def _compute_stats(
        a1: xr.DataArray,
        b1: xr.DataArray,
        a2: xr.DataArray,
        b2: xr.DataArray,
        labels: xr.DataArray,
        quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
        mix_dim: str = "bootstrap_replicate",
        name1: str = "grp1",
        name2: str = "grp2",
    ) -> xr.Dataset:
        """get pvalues on distribution means and pvalue quantiles on samples"""
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
        use_stats_xr: xr.DataArray = xr.DataArray(
            [constants.STATS_AVAILABLE[x] for x in use_stats], [("stats", use_stats)]
        )
        # quantiles to input
        quantiles_xr = xr.DataArray(quantiles, [("pval_quantile", quantiles)])
        # concatenated distribution parameters
        a = xr.concat(
            [a1.reset_index("prefix", drop=True), a2.reset_index("prefix", drop=True)],
            "prefix",
        ).chunk({"prefix": None})
        b = xr.concat(
            [b1.reset_index("prefix", drop=True), b2.reset_index("prefix", drop=True)],
            "prefix",
        ).chunk({"prefix": None})
        # make sure they have mixture dimension
        if mix_dim not in a.dims:
            a = a.expand_dims(**{mix_dim: 1})
        if mix_dim not in b.dims:
            b = b.expand_dims(**{mix_dim: 1})
        # return result
        result = xr.apply_ufunc(
            bm.stats_sample,
            a,
            b,
            labels,
            quantiles_xr,
            psisamples,
            use_stats_xr,
            input_core_dims=[
                ["prefix", mix_dim],
                ["prefix", mix_dim],
                ["prefix"],
                ["pval_quantile"],
                [],
                ["stats"],
            ],
            output_core_dims=[["stats"], ["stats", "pval_quantile"]],
            output_dtypes=[np.float64, np.float64],
            dask="parallelized",
        )
        ds = xr.Dataset(
            {
                "pvalue": result[0].assign_attrs(
                    name1=name1,
                    name2=name2,
                    grp1_prefix=a1.prefix.values.tolist(),
                    grp2_prefix=a2.prefix.values.tolist(),
                ),
                "pvalue_quantiles": result[1].assign_attrs(
                    name1=name1,
                    name2=name2,
                    grp1_prefix=a1.prefix.values.tolist(),
                    grp2_prefix=a2.prefix.values.tolist(),
                    psisamples=psisamples,
                ),
            },
        )
        if not (psisamples > 0 and quantiles):
            ds = ds.drop_dims("pval_quantile")
        return ds

    def raw_stats(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
    ) -> xr.Dataset:
        """Statistics on means, samples from raw posteriors"""
        return self._compute_stats(
            self.psi1.raw_alpha,
            self.psi1.raw_beta,
            self.psi2.raw_alpha,
            self.psi2.raw_beta,
            self.labels,
            quantiles=quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
        )

    def bootstrap_stats(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
    ) -> xr.Dataset:
        """Statistics on means, samples from bootstrap posteriors"""
        return self._compute_stats(
            self.psi1.bootstrap_alpha,
            self.psi1.bootstrap_beta,
            self.psi2.bootstrap_alpha,
            self.psi2.bootstrap_beta,
            self.labels,
            quantiles=quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
        )

    def approximate_stats(
        self,
        quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
    ) -> xr.Dataset:
        """Statistics on means, samples from approximate posteriors"""
        return self._compute_stats(
            self.psi1.approximate_alpha,
            self.psi1.approximate_beta,
            self.psi2.approximate_alpha,
            self.psi2.approximate_beta,
            self.labels,
            quantiles=quantiles,
            psisamples=psisamples,
            use_stats=use_stats,
        )

    def dataset(
        self,
        raw_psi: bool = True,
        bootstrap_psi: bool = True,
        raw_stats: bool = False,
        bootstrap_stats: bool = False,
        approximate_stats: bool = True,
        population_quantiles: Sequence[
            float
        ] = constants.DEFAULT_HET_POPULATION_QUANTILES,
        pvalue_quantiles: Sequence[float] = constants.DEFAULT_HET_PVALUE_QUANTILES,
        psisamples: int = constants.DEFAULT_HET_PSISAMPLES,
        use_stats: Union[str, Collection[str]] = constants.DEFAULT_HET_USESTATS,
    ) -> xr.Dataset:
        combine_ds: List[Union[xr.DataArray, xr.Dataset]] = [
            self.passed.rename("passed"),
            xr.concat(
                [
                    self.psi1.num_passed.expand_dims(grp=[self.name1]),
                    self.psi2.num_passed.expand_dims(grp=[self.name2]),
                ],
                dim="grp",
            ).rename("num_passed"),
        ]
        if raw_stats:
            combine_ds.append(
                self.raw_stats(
                    quantiles=pvalue_quantiles,
                    psisamples=psisamples,
                    use_stats=use_stats,
                ).pipe(lambda x: x.rename_vars(**{y: f"raw_{y}" for y in x.data_vars}))
            )
        if bootstrap_stats:
            combine_ds.append(
                self.bootstrap_stats(
                    quantiles=pvalue_quantiles,
                    psisamples=psisamples,
                    use_stats=use_stats,
                ).pipe(
                    lambda x: x.rename_vars(
                        **{y: f"bootstrap_{y}" for y in x.data_vars}
                    )
                )
            )
        if approximate_stats:
            combine_ds.append(
                self.approximate_stats(
                    quantiles=pvalue_quantiles,
                    psisamples=psisamples,
                    use_stats=use_stats,
                ).pipe(
                    lambda x: x.rename_vars(**{y: f"approx_{y}" for y in x.data_vars})
                )
            )
        if raw_psi:
            combine_ds.append(
                xr.concat(
                    [
                        self.psi1.raw_psi_mean_population_median.expand_dims(
                            grp=[self.name1]
                        ),
                        self.psi2.raw_psi_mean_population_median.expand_dims(
                            grp=[self.name2]
                        ),
                    ],
                    dim="grp",
                ).rename("raw_psi_median")
            )
            if population_quantiles:
                combine_ds.append(
                    xr.concat(
                        [
                            self.psi1.raw_psi_mean_population_quantile(
                                quantiles=population_quantiles
                            ).expand_dims(grp=[self.name1]),
                            self.psi2.raw_psi_mean_population_quantile(
                                quantiles=population_quantiles
                            ).expand_dims(grp=[self.name2]),
                        ],
                        dim="grp",
                    ).rename("raw_psi_quantile")
                )
        if bootstrap_psi:
            combine_ds.append(
                xr.concat(
                    [
                        self.psi1.bootstrap_psi_mean_population_median.expand_dims(
                            grp=[self.name1]
                        ),
                        self.psi2.bootstrap_psi_mean_population_median.expand_dims(
                            grp=[self.name2]
                        ),
                    ],
                    dim="grp",
                ).rename("bootstrap_psi_median")
            )
            if population_quantiles:
                combine_ds.append(
                    xr.concat(
                        [
                            self.psi1.bootstrap_psi_mean_population_quantile(
                                quantiles=population_quantiles
                            ).expand_dims(grp=[self.name1]),
                            self.psi2.bootstrap_psi_mean_population_quantile(
                                quantiles=population_quantiles
                            ).expand_dims(grp=[self.name2]),
                        ],
                        dim="grp",
                    ).rename("bootstrap_psi_quantile")
                )
        return xr.merge(combine_ds, compat="override", join="exact").reset_coords(
            drop=True
        )
