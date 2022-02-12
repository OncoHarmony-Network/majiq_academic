"""
GeneModules.py

Identification of and functions on splicegraph modules.

Author: Joseph K Aicher
"""

import numpy as np
import numpy.typing as npt
import xarray as xr

from new_majiq.internals import GeneModules as _GeneModules

from .ExonConnections import ExonConnections
from .GeneRegions import GeneRegions
from .SpliceGraphMask import SpliceGraphMask


class GeneModules(GeneRegions):
    """Collection of splicing modules per gene.

    Collection of splicing modules per gene.
    Coordinates reflect extreme values of start and end for connections in
    module.
    Additional data on top of `GeneRegions` for first and last exons per module
    (inclusive) and index of first intron/junction and number of unmasked
    introns/junctions.

    Parameters
    ----------
    gene_modules: _GeneModules
        Internals representation of gene modules
    """

    def __init__(self, gene_modules: _GeneModules):
        super().__init__(gene_modules)
        return

    @classmethod
    def from_connections_and_mask(
        self, exon_connections: ExonConnections, sg_mask: SpliceGraphMask
    ) -> "GeneModules":
        """Create :class:`GeneModules` from :class:`ExonConnections` and :class:`SpliceGraphMask`"""
        return GeneModules(
            _GeneModules(exon_connections._exon_connections, sg_mask._sg_mask)
        )

    @property
    def _gene_modules(self) -> _GeneModules:
        """expose underlying internals representation"""
        return self._gene_regions

    @property
    def exon_connections(self) -> ExonConnections:
        return ExonConnections(self._gene_modules.exon_connections)

    @property
    def mask(self) -> SpliceGraphMask:
        return SpliceGraphMask(self._gene_modules.mask)

    @property
    def start_exon_idx(self) -> npt.NDArray[np.uint64]:
        """Index of exon the splicing module starts from."""
        return self._gene_modules.start_exon_idx

    @property
    def end_exon_idx(self) -> npt.NDArray[np.uint64]:
        """Index of exon the splicing module ends at."""
        return self._gene_modules.end_exon_idx

    @property
    def start_intron_idx(self) -> npt.NDArray[np.uint64]:
        """Index of first intron in the splicing module."""
        return self._gene_modules.start_intron_idx

    @property
    def num_introns(self) -> npt.NDArray[np.uint64]:
        """Number of unmasked introns that are part of the splicing module"""
        return self._gene_modules.num_introns

    @property
    def start_junction_idx(self) -> npt.NDArray[np.uint64]:
        """Index of first junction in the splicing module."""
        return self._gene_modules.start_junction_idx

    @property
    def num_junctions(self) -> npt.NDArray[np.uint64]:
        """Number of unmasked junctions that are part of the splicing module"""
        return self._gene_modules.num_junctions

    @property
    def introns_module_idx(self) -> npt.NDArray[np.uint64]:
        """Index of module to which each intron belongs (if not masked)"""
        return self._gene_modules.introns_module_idx

    @property
    def junctions_module_idx(self) -> npt.NDArray[np.uint64]:
        """Index of module to which each junction belongs (if not masked)"""
        return self._gene_modules.junctions_module_idx

    @property
    def df(self) -> xr.Dataset:
        """view on underlying modules as xarray Dataset"""
        return xr.Dataset(
            {},
            {
                "gene_idx": ("gm_idx", self.gene_idx),
                "start": ("gm_idx", self.start),
                "end": ("gm_idx", self.end),
                "start_exon_idx": ("gm_idx", self.start_exon_idx),
                "end_exon_idx": ("gm_idx", self.end_exon_idx),
                "start_intron_idx": ("gm_idx", self.start_intron_idx),
                "num_introns": ("gm_idx", self.num_introns),
                "start_junction_idx": ("gm_idx", self.start_junction_idx),
                "num_junctions": ("gm_idx", self.num_junctions),
                "introns_gm_idx": ("gi_idx", self.introns_module_idx),
                "junctions_gm_idx": ("gj_idx", self.junctions_module_idx),
            },
        )

    @property
    def df_introns(self) -> xr.Dataset:
        """Dataset of introns annotated by mask and module status"""
        return self.mask.df_introns.assign(
            module_idx=("gi_idx", self.introns_module_idx)
        )

    @property
    def df_junctions(self) -> xr.Dataset:
        """Dataset of junctions annotated by mask and module status"""
        return self.mask.df_junctions.assign(
            module_idx=("gj_idx", self.junctions_module_idx)
        )

    def event_is_masked(
        self, ref_exon_idx: npt._ArrayLikeInt_co, event_type: npt._ArrayLikeStr_co
    ) -> npt.NDArray[np.bool_]:
        """Identify if selected events would be entirely masked."""
        return self._gene_modules.event_is_masked(
            ref_exon_idx, ExonConnections._event_type_is_source(event_type)
        )

    def event_module_idx(
        self, ref_exon_idx: npt._ArrayLikeInt_co, event_type: npt._ArrayLikeStr_co
    ) -> npt.NDArray[np.uint64]:
        """Identify module to which selected events belong."""
        return self._gene_modules.event_module_idx(
            ref_exon_idx, ExonConnections._event_type_is_source(event_type)
        )
