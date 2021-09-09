"""
PassedJunctions.py

Accumulate experiment junctions towards passing them (GroupJunctionsGenerator)
Accumulate groups of junctions, create new junctions that have passed
(PassedJunctionsGenerator)

Author: Joseph K Aicher
"""

from typing import Final

import new_majiq.constants as constants
from new_majiq.Exons import Exons
from new_majiq.GeneJunctions import GeneJunctions
from new_majiq.internals import ExperimentThresholds
from new_majiq.internals import GroupJunctionsGenerator as _GroupJunctionsGenerator
from new_majiq.internals import PassedJunctionsGenerator as _PassedJunctionsGenerator
from new_majiq.SJJunctionsBins import SJJunctionsBins


class GroupJunctionsGenerator(object):
    def __init__(self, junctions: GeneJunctions, exons: Exons):
        self._group: Final[_GroupJunctionsGenerator] = _GroupJunctionsGenerator(
            junctions._gene_junctions, exons._exons
        )
        return

    @property
    def num_experiments(self) -> int:
        return self._group.num_experiments

    @property
    def num_known(self) -> int:
        """Number of 'known' junctions (originally passed in)"""
        return self._group.num_known

    @property
    def num_denovo(self) -> int:
        """Number of potential denovo junctions (excluding known denovos)"""
        return self._group.num_denovo

    def add_experiment(
        self,
        sj_junctions: SJJunctionsBins,
        thresholds: ExperimentThresholds = constants.DEFAULT_BUILD_EXP_THRESHOLDS,
        add_denovo: bool = constants.DEFAULT_BUILD_DENOVO_JUNCTIONS,
    ) -> "GroupJunctionsGenerator":
        """Add experiment to build group"""
        self._group.add_experiment(
            sj_junctions._sj_junctionsbins, thresholds, add_denovo
        )
        return self


class PassedJunctionsGenerator(object):
    def __init__(self, junctions: GeneJunctions):
        self._passed: Final[_PassedJunctionsGenerator] = _PassedJunctionsGenerator(
            junctions._gene_junctions
        )
        return

    @property
    def num_known(self) -> int:
        """Number of 'known' junctions (originally passed in)"""
        return self._passed.num_known

    @property
    def num_denovo(self) -> int:
        """Number of denovo junctions passing filters (excluding known denovos)"""
        return self._passed.num_denovo

    def add_group(
        self,
        group: GroupJunctionsGenerator,
        min_experiments: float = constants.DEFAULT_BUILD_MINEXPERIMENTS,
    ) -> "PassedJunctionsGenerator":
        """Add group towards passing junctions with enough evidence"""
        self._passed.add_group(group._group, min_experiments)
        return self

    def get_passed(
        self, denovo_simplified: bool = constants.DEFAULT_BUILD_DENOVO_SIMPLIFIED
    ) -> GeneJunctions:
        """Get new GeneJunctions taking into account the groups that were added"""
        return GeneJunctions(self._passed.get_passed(denovo_simplified))
