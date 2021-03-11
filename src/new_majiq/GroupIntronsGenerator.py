"""
GroupIntronsGenerator.py

Build group for passing gene introns

Author: Joseph K Aicher
"""

import numpy as np
import xarray as xr

import new_majiq.constants as constants

from new_majiq.GeneIntrons import GeneIntrons
from new_majiq.SJIntronsBins import SJIntronsBins
from new_majiq.internals import GroupIntronsGenerator as _GroupIntronsGenerator
from new_majiq.internals import ExperimentThresholds

from typing import (
    Final,
)


class GroupIntronsGenerator(object):
    """Count number of experiments that have passed each input intron"""

    def __init__(self, introns: GeneIntrons):
        self._group: Final[_GroupIntronsGenerator] = (
            _GroupIntronsGenerator(introns._gene_introns)
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
        return (
            self.introns.df
            .assign(num_passed=("gi_idx", self.num_passed))
            .assign_attrs(num_experiments=self.num_experiments)
        )

    def add_experiment(
        self,
        sj_introns: SJIntronsBins,
        thresholds: ExperimentThresholds = constants.DEFAULT_BUILD_EXP_THRESHOLDS,
    ) -> "GroupIntronsGenerator":
        """Add input experiment to build group"""
        self._group.add_experiment(sj_introns._sj_intronsbins, thresholds)
        return self

    def update_introns(
        self, min_experiments: float = constants.DEFAULT_BUILD_MINEXPERIMENTS
    ) -> None:
        """Pass introns in place that have passed enough experiments, reset counts"""
        self._group.update_introns(min_experiments)
        return
