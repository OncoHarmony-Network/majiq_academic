"""
Exons.py

Exons for genes

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union, cast

import numpy as np
import numpy.typing as npt
import xarray as xr

import new_majiq.constants as constants
from new_majiq.internals import Exons as _Exons

from ._workarounds import _load_zerodim_variables
from .GeneRegions import GeneRegions
from .Genes import Genes

if TYPE_CHECKING:
    from .GeneIntrons import GeneIntrons
    from .GeneJunctions import GeneJunctions


class Exons(GeneRegions):
    """Collection of exons per gene and their annotated/updated coordinates"""

    def __init__(self, exons: _Exons):
        super().__init__(exons)
        return

    def infer_with_junctions(self, junctions: "GeneJunctions") -> "Exons":
        """Return updated :py:class:`Exons` accommodating novel junctions per gene

        Parameters
        ----------
        junctions: GeneJunctions
            Junctions over same genes. Novel exons will be added or have
            boundaries extended compared to original annotated exons to match
            the junctions
        """
        from new_majiq.internals import SpliceGraph as _SpliceGraph

        return Exons(_SpliceGraph.infer_exons(self._exons, junctions._gene_junctions))

    def checksum(self):
        return self._exons.checksum()

    @property
    def _exons(self) -> _Exons:
        """expose underlying internals representation of Exons"""
        return self._gene_regions

    @property
    def exon_idx(self) -> npt.NDArray[np.int64]:
        return self._region_idx

    def get_annotated(self) -> "Exons":
        """Annotated :class:`Exons` with original coordinates, excluding novel exons"""
        return Exons(self._exons.get_annotated())

    def potential_introns(
        self, make_simplified: bool = constants.DEFAULT_BUILD_DENOVO_SIMPLIFIED
    ) -> "GeneIntrons":
        """:py:class:`GeneIntrons` enumerating all possible introns between exons

        Return :py:class:`GeneIntrons` enumerating all possible introns between
        exons. All introns will be initially marked as de novo and not passing
        filters (use :py:meth:`GeneIntrons.update_flags_from` to update these
        flags)

        Parameters
        ----------
        make_simplified: bool
            Indicate whether introns should start in the simplified state
            (requires subsequently performing simplification with
            :py:class:`SimplifierGroup` to identify which introns pass
            simplification thresholds)
        """
        from .GeneIntrons import GeneIntrons

        return GeneIntrons(self._exons.potential_introns(make_simplified))

    def empty_introns(self) -> "GeneIntrons":
        """Return empty :py:class:`GeneIntrons` that match these exons"""
        from .GeneIntrons import GeneIntrons

        return GeneIntrons.from_genes(self.genes, self)

    @property
    def annotated_start(self) -> npt.NDArray[np.int64]:
        """Annotated coordinates start. If denovo, -1

        Note that denovo behavior is different than previous versions of MAJIQ
        """
        return self._exons.annotated_start

    @property
    def annotated_end(self) -> npt.NDArray[np.int64]:
        """Annotated coordinates end. If denovo, -1

        Note that denovo behavior is different than previous versions of MAJIQ
        """
        return self._exons.annotated_end

    def is_denovo(
        self,
        exon_idx: Optional[npt._ArrayLikeInt_co] = None,
        annotated_exons: Optional["Exons"] = None,
    ) -> npt.NDArray[np.bool_]:
        """Return denovo status of selected exons

        Parameters
        ----------
        exon_idx: Optional[array_like[int]]
            Index into exons for to get denovo status for.
            If None, get denovo status for all exons.
        annotated_exons: Optional[Exons]
            If specified, use exons that are found in `annotated_exons` to
            reduce the number of exons called annotated by treating exons
            overlapping with `annotated_exons` as annotated.
        """
        if exon_idx is None:
            exon_idx = self.exon_idx
        denovos = self._exons.is_denovo(exon_idx)
        if annotated_exons:
            denovos = denovos & ~self.overlaps(annotated_exons, exon_idx)
        return denovos

    def is_exon_extension(
        self, exon_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.bool_]:
        """Return if exon(s) have exon extension"""
        if exon_idx is None:
            exon_idx = self.exon_idx
        return self._exons.is_exon_extension(exon_idx)

    def is_full_exon(
        self, exon_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.bool_]:
        """Return if exon(s) are full exons"""
        if exon_idx is None:
            exon_idx = self.exon_idx
        return self._exons.is_full_exon(exon_idx)

    def is_half_exon(
        self, exon_idx: Optional[npt._ArrayLikeInt_co] = None
    ) -> npt.NDArray[np.bool_]:
        """Return if exon(s) are half exons"""
        if exon_idx is None:
            exon_idx = self.exon_idx
        return self._exons.is_half_exon(exon_idx)

    @property
    def df(self) -> xr.Dataset:
        """view on underlying exons as xarray Dataset"""
        return xr.Dataset(
            {},
            {
                "exon_idx": self.exon_idx,
                "gene_idx": ("exon_idx", self.gene_idx),
                "start": ("exon_idx", self.start),
                "end": ("exon_idx", self.end),
                "annotated_start": ("exon_idx", self.annotated_start),
                "annotated_end": ("exon_idx", self.annotated_end),
            },
        )

    def to_zarr(
        self, path: Union[str, Path], mode: str, consolidated: bool = True
    ) -> None:
        """Serialize to zarr format. Note genes need to be saved separately"""
        self.df.drop_vars("exon_idx").pipe(lambda x: x.chunk(x.sizes)).pipe(
            _load_zerodim_variables
        ).to_zarr(path, mode=mode, group=constants.NC_EXONS, consolidated=consolidated)
        return

    @classmethod
    def from_arrays(
        cls,
        genes: Genes,
        gene_idx: npt._ArrayLikeInt_co,
        start: npt._ArrayLikeInt_co,
        end: npt._ArrayLikeInt_co,
        annotated_start: Optional[npt._ArrayLikeInt_co] = None,
        annotated_end: Optional[npt._ArrayLikeInt_co] = None,
    ) -> "Exons":
        """Create :class:`Exons` from :class:`Genes` and input arrays"""
        start = np.array(start, copy=False)
        end = np.array(end, copy=False)
        if annotated_start is None and annotated_end is None:
            # assume annotated unless not full exon
            full_exon = (start >= 0) & (end >= 0)
            annotated_start = np.where(full_exon, start, -1)
            annotated_end = np.where(full_exon, end, -1)
        elif (annotated_start is None) != (annotated_end is None):
            raise ValueError("annotated_{start,end} must both or neither be passed")
        annotated_start = cast(npt._ArrayLikeInt_co, annotated_start)
        annotated_end = cast(npt._ArrayLikeInt_co, annotated_end)
        return Exons(
            _Exons(genes._genes, gene_idx, start, end, annotated_start, annotated_end)
        )

    @classmethod
    def from_zarr(
        self,
        path: Union[str, Path],
        genes: Optional[Genes] = None,
    ) -> "Exons":
        """Read exons from zarr file

        Parameters
        ----------
        path: Union[str, Path]
            path to zarr file
        genes: Optional[Genes]
            genes on which the exons are defined. If None, try loading from
            zarr file. Note that new_majiq checks if objects refer to the
            same genes (not that they are identical), so it is usually
            desired to provide the variable than using the default behavior
        """
        if genes is None:
            genes = Genes.from_zarr(path)
        with xr.open_zarr(path, group=constants.NC_EXONS) as df:
            df.load()
            return Exons.from_arrays(
                genes,
                df.gene_idx.values,
                df.start.values,
                df.end.values,
                annotated_start=df.annotated_start.values,
                annotated_end=df.annotated_end.values,
            )
