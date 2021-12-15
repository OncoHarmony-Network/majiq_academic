"""
GeneJunctionsAccumulator.py

Accumulator of :class:`GeneJunctions` objects to create single one with unique
junctions, combining flags from inputs that had the junctions.

Author: Joseph K Aicher
"""

from typing import Final

from new_majiq.internals import GeneJunctionsAccumulator as _GeneJunctionsAccumulator

from .GeneJunctions import GeneJunctions
from .Genes import Genes


class GeneJunctionsAccumulator(object):
    """Accumulate :class:`GeneJunctions` objects into combined :class:`GeneJunctions`

    Parameters
    ----------
    genes: Genes
        Underlying genes that must be shared by all GeneJunctions input into
        the accumulator
    """

    def __init__(self, genes: Genes):
        self._accumulator: Final[_GeneJunctionsAccumulator] = _GeneJunctionsAccumulator(
            genes._genes
        )
        return

    @property
    def genes(self) -> Genes:
        """:class:`Genes` shared by input and accumulated junctions"""
        return Genes(self._accumulator._genes)

    def add(self, junctions: GeneJunctions, make_annotated: bool = False) -> None:
        """Add :class:`GeneJunctions` to accumulator

        Add :class:`GeneJunctions` to accumulator, optionally treating them all
        as annotated. Junctions that have been seen before are updated such
        that a junction passes the build if it passed in any of the inputs and
        a junction is denovo or simplified only if it is so in all of the
        inputs.

        Parameters
        ----------
        junctions: GeneJunctions
            junctions to add
        make_annotated: bool
            If True, treat all input junctions as annotated
        """
        self._accumulator.add(junctions._gene_junctions, make_annotated=make_annotated)
        return

    def accumulated(self) -> GeneJunctions:
        """Return :class:`GeneJunctions` with previously added junctions

        Returns
        -------
        GeneJunctions
            Junctions that were added (:meth:`GeneJunctions.add`) with merged
            flags
        """
        return GeneJunctions(self._accumulator.accumulated())
