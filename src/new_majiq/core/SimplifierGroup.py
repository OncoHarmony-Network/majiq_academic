"""
SimplifierGroup.py

Author: Joseph K Aicher
"""

from typing import TYPE_CHECKING, Final

import numpy as np
import xarray as xr

import new_majiq.constants as constants
from new_majiq.internals import SimplifierGroup as _SimplifierGroup

from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions
from .SpliceGraphReads import SpliceGraphReads

if TYPE_CHECKING:
    from .SJExperiment import SJExperiment


class SimplifierGroup(object):
    """Accumulate groups of experiments to determine connections for simplification"""

    def __init__(self, group: _SimplifierGroup):
        self._group: Final[_SimplifierGroup] = group
        return

    @property
    def num_experiments(self) -> int:
        """Number of experiments added to group"""
        return self._group.num_experiments

    @property
    def introns(self) -> GeneIntrons:
        return GeneIntrons(self._group._exon_connections._introns)

    @property
    def introns_passed_src(self) -> np.ndarray:
        return self._group.introns_passed_src

    @property
    def introns_passed_dst(self) -> np.ndarray:
        return self._group.introns_passed_dst

    @property
    def junctions(self) -> GeneJunctions:
        return GeneJunctions(self._group._exon_connections._junctions)

    @property
    def junctions_passed_src(self) -> np.ndarray:
        return self._group.junctions_passed_src

    @property
    def junctions_passed_dst(self) -> np.ndarray:
        return self._group.junctions_passed_dst

    @property
    def df_introns(self) -> xr.Dataset:
        return self.introns.df.assign(
            passed_src=("gi_idx", self.introns_passed_src),
            passed_dst=("gi_idx", self.introns_passed_dst),
        ).assign_attrs(num_experiments=self.num_experiments)

    @property
    def df_junctions(self) -> xr.Dataset:
        return self.junctions.df.assign(
            passed_src=("gj_idx", self.junctions_passed_src),
            passed_dst=("gj_idx", self.junctions_passed_dst),
        ).assign_attrs(num_experiments=self.num_experiments)

    def add_experiment(
        self,
        sj: "SJExperiment",
        min_psi: float = constants.DEFAULT_SIMPLIFIER_MINPSI,
        minreads_annotated: float = constants.DEFAULT_SIMPLIFIER_MINREADS_ANNOTATED,
        minreads_denovo: float = constants.DEFAULT_SIMPLIFIER_MINREADS_DENOVO,
        minreads_introns: float = constants.DEFAULT_SIMPLIFIER_MINREADS_INTRON,
    ) -> "SimplifierGroup":
        """Add reads from experiment to group for simplification"""
        self._group.add_experiment(
            SpliceGraphReads._internals_from_connections_and_sj(
                self.introns, self.junctions, sj
            ),
            min_psi,
            minreads_annotated,
            minreads_denovo,
            minreads_introns,
        )
        return self

    def update_connections(
        self, min_experiments: float = constants.DEFAULT_SIMPLIFIER_MINEXPERIMENTS
    ) -> None:
        """Unsimplify introns/junctions with enough evidence, reset counts"""
        self._group.update_connections(min_experiments)
        return
