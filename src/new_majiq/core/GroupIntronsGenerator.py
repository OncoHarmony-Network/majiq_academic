"""
GroupIntronsGenerator.py

Build group for passing gene introns

Author: Joseph K Aicher
"""

from typing import Final

import numpy as np
import xarray as xr

import new_majiq.constants as constants
from new_majiq.internals import ExperimentThresholds
from new_majiq.internals import GroupIntronsGenerator as _GroupIntronsGenerator

from .GeneIntrons import GeneIntrons
from .SJIntronsBins import SJIntronsBins


class GroupIntronsGenerator(object):
    """Accumulator of :py:class:`SJIntronsBins` that pass per-experiment thresholds to update :py:class:`GeneIntrons` flags"""

    def __init__(self, introns: GeneIntrons):
        self._group: Final[_GroupIntronsGenerator] = _GroupIntronsGenerator(
            introns._gene_introns
        )
        return

    @property
    def num_experiments(self) -> int:
        """Number of experiments in current group"""
        return self._group.num_experiments

    @property
    def introns(self) -> GeneIntrons:
        return GeneIntrons(self._group._introns)

    @property
    def num_passed(self) -> np.ndarray:
        return self._group.num_passed

    @property
    def df(self) -> xr.Dataset:
        return self.introns.df.assign(
            num_passed=("gi_idx", self.num_passed)
        ).assign_attrs(num_experiments=self.num_experiments)

    def add_experiment(
        self,
        sj_introns: SJIntronsBins,
        thresholds: ExperimentThresholds = constants.DEFAULT_BUILD_EXP_THRESHOLDS,
    ) -> "GroupIntronsGenerator":
        """Add :py:class:`SJIntronsBins` experiment to build group

        Parameters
        ----------
        sj_introns: SJIntronsBins
            Per-bin read coverage over introns
        thresholds: ExperimentThresholds
            Per-experiment thresholds to determine if there is enough evidence
            for each intron
        """
        self._group.add_experiment(sj_introns._sj_intronsbins, thresholds)
        return self

    def update_introns(
        self, min_experiments: float = constants.DEFAULT_BUILD_MINEXPERIMENTS
    ) -> None:
        """In-place update of original :py:class:`GeneIntrons` flags passing group filters

        Parameters
        ----------
        min_experiments: float
            Threshold for group filters. This specifies the fraction (value <
            1) or absolute number (value >= 1) of experiments that must pass
            individually in the build group that must pass in order for the
            intron to be updated as being reliable
        """
        self._group.update_introns(min_experiments)
        return
