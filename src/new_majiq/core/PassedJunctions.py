"""
PassedJunctions.py

Accumulate experiment junctions towards passing them (GroupJunctionsGenerator)
Accumulate groups of junctions, create new junctions that have passed
(PassedJunctionsGenerator)

Author: Joseph K Aicher
"""

from typing import Final

import new_majiq.constants as constants
from new_majiq.internals import ExperimentThresholds
from new_majiq.internals import GroupJunctionsGenerator as _GroupJunctionsGenerator
from new_majiq.internals import PassedJunctionsGenerator as _PassedJunctionsGenerator

from .Exons import Exons
from .GeneJunctions import GeneJunctions
from .SJJunctionsBins import SJJunctionsBins


class GroupJunctionsGenerator(object):
    """Accumulator of :py:class:`SJJunctionsBins` that pass per-experiment thresholds"""

    def __init__(self, junctions: GeneJunctions, exons: Exons):
        self._group: Final[_GroupJunctionsGenerator] = _GroupJunctionsGenerator(
            junctions._gene_junctions, exons._exons
        )
        return

    @property
    def _known(self) -> GeneJunctions:
        """Underlying known :class:`GeneJunctions`"""
        return GeneJunctions(self._group._known)

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
        """Add :py:class:`SJJunctionsBins` experiment to build group

        Parameters
        ----------
        sj_junctions: SJJunctionsBins
            Per-bin read coverage over junctions
        thresholds: ExperimentThresholds
            Per-experiment thresholds to determine if there is enough evidence
            for each junction
        add_denovo: bool
            Indicate whether to evaluate novel junctions from sj_junctions
        """
        self._group.add_experiment(
            sj_junctions._sj_junctionsbins, thresholds, add_denovo
        )
        return self


class PassedJunctionsGenerator(object):
    """Accumulator of :py:class:`GroupJunctionsGenerator` to create updated :py:class:`GeneJunctions`"""

    def __init__(self, junctions: GeneJunctions):
        self._passed: Final[_PassedJunctionsGenerator] = _PassedJunctionsGenerator(
            junctions._gene_junctions
        )
        return

    @property
    def _known(self) -> GeneJunctions:
        """Underlying known :class:`GeneJunctions`"""
        return GeneJunctions(self._passed._known)

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
        """Update passed junctions with :py:class:`GroupJunctionsGenerator` build group

        Parameters
        ----------
        group: GroupJunctionsGenerator
            Accumulator of how many times each junction passed per-experiment
            thresholds in a build group
        min_experiments: float
            Threshold for group filters. This specifies the fraction (value <
            1) or absolute number (value >= 1) of experiments that must pass
            individually in the build group that must pass in order for the
            junction to be updated as being reliable
        """
        self._passed.add_group(group._group, min_experiments)
        return self

    def add_junction(
        self, gene_idx: int, start: int, end: int
    ) -> "PassedJunctionsGenerator":
        """Pass the specified junction

        Parameters
        ----------
        gene_idx: int
            gene index for junction to be passed
        start, end: int
            Coordinates for junction to be passed
        """
        num_genes = len(self._known.genes)
        if not 0 <= gene_idx < num_genes:
            raise ValueError(
                f"Invalid gene_idx ({gene_idx}, must be in [0, {num_genes}))"
            )
        gene_start = self._known.genes.start[gene_idx]
        gene_end = self._known.genes.end[gene_idx]
        if start < gene_start or gene_end < end:
            raise ValueError(
                f"Invalid coordinates ({start}, {end});"
                f" must be in [{gene_start}, {gene_end}]"
            )
        self._passed.add_junction(gene_idx, start, end)
        return self

    def get_passed(
        self, denovo_simplified: bool = constants.DEFAULT_BUILD_DENOVO_SIMPLIFIED
    ) -> GeneJunctions:
        """Return :py:class:`GeneJunctions` with updated flags and novel junctions

        Parameters
        ----------
        denovo_simplified: bool
            Indicate whether novel junctions should start in the simplified
            state (requires subsequently performing simplification with
            :py:class:`SimplifierGroup` to identify which junctions pass
            simplification thresholds)
        """
        return GeneJunctions(self._passed.get_passed(denovo_simplified))
