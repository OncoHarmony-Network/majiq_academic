"""
mplSpliceGraph.py

Convenience functions for plotting splicegraphs with matplotlib

Author: Joseph K Aicher
"""

from typing import Callable, Dict, Final, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.collections import PatchCollection
from matplotlib.patches import Arc, Rectangle
from scipy.interpolate import interp1d

from ..core.SpliceGraph import SpliceGraph
from ..core.SpliceGraphReads import SpliceGraphReads

DEFAULT_EXON_HEIGHT = 1.0
DEFAULT_INTRON_HEIGHT = 0.3
DEFAULT_HALF_EXON_WIDTH = 0
DENOVO_COLORS = {True: "green", False: "red"}


class SpliceGraphGeneView(object):
    """Prepare tables for plotting information about a specific gene

    Parameters
    ----------
    sg: SpliceGraph
        SpliceGraph with gene to be plotted.
    gene_idx: int
        Which gene to plot.
    sg_reads: Optional[SpliceGraphReads]
        Information on per-intron/junction coverage. If there are multiple
        prefixes, use reads from the first prefix (select the reads you want to
        plot beforehand).
    show_simplified: bool
        If True, plot all introns/junctions, including those that were
        simplified.
    """

    @classmethod
    def with_gene_name(
        cls, sg: SpliceGraph, gene_name: str, **kwargs
    ) -> "SpliceGraphGeneView":
        """Create :class:`SpliceGraphGeneView` for gene with specified gene_name

        Notes
        -----
        It is possible that sg has multiple genes with the same name. This will
        plot the first gene it finds matching the specified name.
        """
        return SpliceGraphGeneView(sg, sg.genes.gene_name.index(gene_name), **kwargs)

    @classmethod
    def with_gene_id(
        cls, sg: SpliceGraph, gene_id: str, **kwargs
    ) -> "SpliceGraphGeneView":
        """Create :class:`SpliceGraphGeneView` for gene with specified gene_id"""
        return SpliceGraphGeneView(sg, sg.genes.gene_id.index(gene_id), **kwargs)

    def __init__(
        self,
        sg: SpliceGraph,
        gene_idx: int,
        sg_reads: Optional[SpliceGraphReads] = None,
        show_simplified: bool = False,
    ):
        """Create :class:`SpliceGraphGeneView` for specified gene"""
        # store relevant information about the gene
        self.seqid: Final[str] = sg.contigs.seqid[sg.genes.contig_idx[gene_idx]]
        self.strand: Final[str] = sg.genes.strand[gene_idx].decode()
        self.gene_id: Final[str] = sg.genes.gene_id[gene_idx]
        self.gene_name: Final[str] = sg.genes.gene_name[gene_idx]
        # get junction coordinates
        self.junctions: Final[pd.DataFrame] = (
            sg.junctions.df.isel(gj_idx=sg.junctions.slice_for_gene(gene_idx))
            .to_dataframe()
            .drop(columns=["gene_idx", "start_exon_idx", "end_exon_idx"])
            .pipe(lambda df: df if show_simplified else df.loc[~df.simplified])
            .assign(
                height=lambda df: _junction_heights(df.start.values, df.end.values),
            )
            .pipe(
                lambda df: df
                if sg_reads is None
                else df.assign(
                    junction_reads=sg_reads.junctions_reads.isel(
                        gj_idx=df.index, prefix=0
                    ).values
                )
            )
        )
        self.introns: Final[pd.DataFrame] = (
            sg.introns.df.isel(gi_idx=sg.introns.slice_for_gene(gene_idx))
            .to_dataframe()
            .drop(columns=["gene_idx", "start_exon_idx", "end_exon_idx"])
            .pipe(lambda df: df if show_simplified else df.loc[~df.simplified])
            .pipe(
                lambda df: df
                if sg_reads is None
                else df.assign(
                    intron_reads=sg_reads.introns_reads.isel(
                        gi_idx=df.index, prefix=0
                    ).values
                )
            )
        )
        self.exons: Final[pd.DataFrame] = (
            sg.exons.df.isel(exon_idx=sg.exons.slice_for_gene(gene_idx))
            .to_dataframe()
            .assign(
                denovo=lambda df: sg.exons.is_denovo(df.index),
                half_exon=lambda df: sg.exons.is_half_exon(df.index),
                plot_start=lambda df: df.start.where(df.start >= 0, df.end),
                plot_end=lambda df: df.end.where(df.end >= 0, df.start),
                annotated_start=lambda df: df.annotated_start.where(
                    df.annotated_start >= 0, df.plot_start
                ),
                annotated_end=lambda df: df.annotated_end.where(
                    df.annotated_end >= 0, df.plot_end
                ),
            )
        )
        return

    def plot_mpl(
        self,
        ax: Optional[Axes] = None,
        view_start: Optional[float] = None,
        view_end: Optional[float] = None,
        rescale_kwargs: Optional[Dict] = None,
        add_title: bool = True,
    ) -> Axes:
        """Plot splicegraph over given coordinates, optionally rescaling

        Parameters
        ----------
        use_ax: Optional[Axes]
            Axes on which exons will be plotted. If None, uses pyplot.gca() to
            get current axes (creating one if necessary)
        view_start, view_end: Optional[float]
            Range of coordinates we want to plot. If None, use start/end
            coordinates of the gene.
        rescale_kwargs: Optional[Dict]
            Override keyword arguments to `fit_reparameterize`. If None, do not
            reparameterize; if Dict, do reparameterize.
        add_title: bool
            Add title describing gene to plot
        """
        _ax: Axes = plt.gca() if ax is None else ax
        rescaler = (
            _no_rescaler
            if rescale_kwargs is None
            else fit_reparameterize(
                self.exons.start.values,
                self.exons.end.values,
                **rescale_kwargs,
            )
        )
        plot_junctions(
            _ax,
            rescaler(self.junctions.start.values),
            rescaler(self.junctions.end.values),
            self.junctions.height.values,
            self.junctions.denovo.values,
            (
                self.junctions.junction_reads > 0
                if "junction_reads" in self.junctions.columns
                else self.junctions.passed_build
            ).values,
        )
        junction_height = 0.5
        if len(self.junctions):
            junction_height *= self.junctions.height.max()
        exon_height = min(DEFAULT_EXON_HEIGHT, 0.3 * junction_height)
        plot_introns(
            _ax,
            rescaler(self.introns.start.values - 1),
            rescaler(self.introns.end.values + 1),
            self.introns.denovo.values,
            (
                self.introns.intron_reads > 0
                if "intron_reads" in self.introns.columns
                else self.introns.passed_build
            ).values,
            exon_height=exon_height,
            intron_height=exon_height * DEFAULT_INTRON_HEIGHT,
        )
        plot_exons(
            _ax,
            rescaler(self.exons.plot_start.values),
            rescaler(self.exons.plot_end.values),
            self.exons.denovo,
            rescaler(self.exons.annotated_start.values),
            rescaler(self.exons.annotated_end.values),
            exon_height=exon_height,
        )
        start: float = view_start if view_start else self.exons.plot_start.min()
        end: float = view_end if view_end else self.exons.plot_end.max()
        start, end = rescaler([start, end])
        buffer = 0.02 * (end - start)
        _ax.set_xlim(start - buffer, end + buffer)
        _ax.set_ylim(-1.5 * exon_height, max(1, 1.2 * junction_height))
        if add_title:
            _ax.set_title(
                f"{self.gene_name} (contig={self.seqid}, gene_id={self.gene_id})"
            )
        return _ax


def plot_exons(
    ax: Axes,
    start: npt.NDArray[Union[np.integer, np.floating]],
    end: npt.NDArray[Union[np.integer, np.floating]],
    denovo: npt.NDArray[np.bool_],
    annotated_start: npt.NDArray[Union[np.integer, np.floating]],
    annotated_end: npt.NDArray[Union[np.integer, np.floating]],
    exon_height: float = DEFAULT_EXON_HEIGHT,
    **collection_kwargs,
) -> None:
    """Plot specified exons

    Parameters
    ----------
    ax: Axes
        Axes on which exons will be plotted
    start, end: array
        Coordinates for each exon
    denovo: array
        Indicate if junction is denovo vs not
    annotated_start, annotated_end: array
        For annotated exons, the original coordinates (without extension)
    exon_height: float
        Height of exon on canvas
    **collection_kwargs
        Additional kwargs to pass into Axes.add_collection().
    """

    def _add_rectangles(x_start, x_end, ec, fc, alpha):
        width = x_end - x_start
        ax.add_collection(
            PatchCollection(
                [
                    Rectangle(
                        (s, -exon_height),
                        w,
                        exon_height,
                        ec=ec,
                        fc=fc,
                        alpha=alpha,
                    )
                    for s, w in zip(x_start, width)
                ],
                match_original=True,
                **collection_kwargs,
            ),
        )
        return

    # plot exons
    _add_rectangles(start, end, ec="black", fc="lightgray", alpha=0.7)
    # add overlay for denovo exons
    _add_rectangles(
        start[denovo],
        end[denovo],
        ec=DENOVO_COLORS[True],
        fc=DENOVO_COLORS[True],
        alpha=0.4,
    )
    # add overlay for exon extension
    _add_rectangles(
        start[annotated_start > start],
        annotated_start[annotated_start > start],
        ec=DENOVO_COLORS[True],
        fc=DENOVO_COLORS[True],
        alpha=0.4,
    )
    _add_rectangles(
        annotated_end[end > annotated_end],
        end[end > annotated_end],
        ec=DENOVO_COLORS[True],
        fc=DENOVO_COLORS[True],
        alpha=0.4,
    )
    return


def plot_introns(
    ax: Axes,
    start: npt.NDArray[Union[np.integer, np.floating]],
    end: npt.NDArray[Union[np.integer, np.floating]],
    denovo: npt.NDArray[np.bool_],
    evidence: npt.NDArray[np.bool_],
    exon_height: float = DEFAULT_EXON_HEIGHT,
    intron_height: float = DEFAULT_INTRON_HEIGHT,
    **collection_kwargs,
) -> None:
    """Plot specified introns

    Parameters
    ----------
    ax: Axes
        Axes on which introns will be plotted
    start, end: array
        Coordinates for each intron
    denovo: array
        Indicate if junction is denovo vs not
    evidence: array
        Indicate if there was RNA-seq evidence for the junction
    exon_height, intron_height: float
        Height of exons, introns on canvas
    **collection_kwargs
        Additional kwargs to pass into Axes.add_collection().
    """
    width = end - start
    intron_y = -0.5 * (exon_height + intron_height)
    ax.add_collection(
        PatchCollection(
            [
                Rectangle(
                    (s, intron_y),
                    w,
                    intron_height,
                    alpha=0.6,
                    fc=DENOVO_COLORS[d],
                    ec=DENOVO_COLORS[d],
                    ls="-" if e else "--",
                )
                for s, w, d, e in zip(start, width, denovo, evidence)
            ],
            match_original=True,
            **collection_kwargs,
        ),
    )
    return None


def plot_junctions(
    ax: Axes,
    start: npt.NDArray[Union[np.integer, np.floating]],
    end: npt.NDArray[Union[np.integer, np.floating]],
    height: npt.NDArray[Union[np.integer, np.floating]],
    denovo: npt.NDArray[np.bool_],
    evidence: npt.NDArray[np.bool_],
    **collection_kwargs,
) -> None:
    """Plot specified junctions

    Parameters
    ----------
    ax: Axes
        Axes on which junctions will be plotted
    start, end: array
        Coordinates for each junction
    height: array
        Height parameter to use for each junction
    denovo: array
        Indicate if junction is denovo vs not
    evidence: array
        Indicate if there was RNA-seq evidence for the junction
    **collection_kwargs
        Additional kwargs to pass into each Arc
    """
    midpoint = 0.5 * (start + end)
    width = end - start
    for m, w, h, d, e in zip(midpoint, width, height, denovo, evidence):
        ax.add_patch(
            Arc(
                (m, 0),
                w,
                h,
                theta1=0,
                theta2=180,
                ec=DENOVO_COLORS[d],
                ls="-" if e else "--",
                **collection_kwargs,
            )
        )
    return


def _no_rescaler(x: npt._ArrayLikeNumber_co) -> npt.NDArray[np.float64]:
    """Don't reparameterize x (except casting to float)"""
    return np.array(x, dtype=np.float64)


def fit_reparameterize(
    exon_start: npt.NDArray[Union[np.integer, np.floating]],
    exon_end: npt.NDArray[Union[np.integer, np.floating]],
    min_scale: float = 0.2,
    intron_length_scale: float = 0.3,
    intron_rescale: float = 0.8,
) -> Callable[[npt._ArrayLikeNumber_co], npt.NDArray[np.float64]]:
    """Reparameterize coordinates to have less variable feature widths.

    Reparameterize coordinates to have less variable feature widths.
    Given exon coordinates in a gene, return a piecewise linear function valid
    within the range of the gene for which exons and introns have limited
    variability in size.

    Parameters
    ----------
    exon_start, exon_end: array
        Sorted coordinates for exon starts/ends
    min_scale: float
        Set minimum final width of exons in reparameterized coordinates vs
        maximum width of 1. Must be a value between 0 and 1.
    intron_length_scale: float
        Relative scale of introns vs exons when identifying maximum length for
        reparameterization.
        Reparameterization is done with respect to exon and intron lengths; the
        longest feature will have length 1 afterwards.
        Introns are often much longer than exons, this adjusts their length
        beforehand so that more real estate in a plot can be given to exons.
        Must be a value greater than 0 (and typically no greater than 1, too).
    intron_rescale: float
        Additional rescaling on introns after setting the length scale.
        This makes it so that final feature lengths for introns are between
        [min_scale * intron_rescale, intron_rescale].
        Must be a value greater than 0 (and typically no greater than 1, too).

    Returns
    -------
    Callable[[array_like], array[float]]
        Function that takes values between minimum/maximum exon coordinates to
        piecewise linear interpolation of coordinates
    """
    if not (0 < min_scale < 1):
        raise ValueError(f"min_scale must be between 0 and 1 ({min_scale = })")
    if intron_length_scale <= 0:
        raise ValueError(
            f"intron_length_scale must be positive-valued ({intron_length_scale = })"
        )
    if intron_rescale <= 0:
        raise ValueError(
            f"intron_rescale must be positive-valued ({intron_rescale = })"
        )
    # cast to arrays
    exon_start = np.array(exon_start, copy=False)
    exon_end = np.array(exon_end, copy=False)
    # filter out half exons from this process
    full_mask = (exon_start >= 0) & (exon_end >= 0)
    exon_start = exon_start[full_mask]
    exon_end = exon_end[full_mask]
    if len(exon_start) < 2:
        return _no_rescaler
    # breakpoints for linear interpolation
    breakpoints = np.empty(2 * exon_start.shape[0])
    breakpoints[::2] = exon_start
    breakpoints[1::2] = exon_end
    # raw lengths of exons/introns
    raw_lengths = 1 + np.diff(breakpoints)
    # treat intron lengths as a fraction of exon lengths
    scaled_lengths = (
        raw_lengths * np.tile([1, intron_length_scale], len(exon_start))[:-1]
    )
    # normalize scaled lengths so minimum value is min_scale, max value is 1
    minL = scaled_lengths.min()
    maxL = scaled_lengths.max()
    norm_lengths = min_scale + (scaled_lengths - minL) * (1 - min_scale) / (maxL - minL)
    # apply additional rescaling of intron lengths
    norm_lengths *= np.tile([1, intron_rescale], len(exon_start))[:-1]
    # so, reparameterized breakpoints will be
    new_breakpoints = np.empty_like(breakpoints)
    new_breakpoints[0] = 0
    np.cumsum(norm_lengths, out=new_breakpoints[1:])
    # return object that performs desired interpolation as callable
    return interp1d(breakpoints, new_breakpoints, fill_value="extrapolate")


def _junction_heights(
    start: npt.NDArray[Union[np.integer, np.floating]],
    end: npt.NDArray[Union[np.integer, np.floating]],
) -> npt.NDArray[np.int64]:
    """Determine heights of junctions in plotted splicegraph

    Notes
    -----
    Heights are determined by identifying chains of nested junctions (e.g., A
    is enclosed by B is enclosed by C...), and using the size of the largest
    chain starting with each junction.

    Parameters
    ----------
    start, end: Sequence[Union[int, float]]
        Junction coordinates in lexsorted order (start, end)

    Returns
    -------
    array[int]
        Heights that can be used for the introns
    """
    # using separate updated array in case iterative updates of height are
    # desired (adjustments to heights, etc.)
    updated = np.zeros(len(start), dtype=np.bool_)
    height = np.ones(len(start), dtype=np.int64)

    def get_height(i: int) -> int:
        """Get the height recursively, using enclosed junctions"""
        if updated[i]:
            return height[i]
        start_i = start[i]
        end_i = end[i]
        cur_height = height[i]
        # junctions behind `i`: subset if and only if same start
        for j in range(i - 1, -1, -1):
            if start[j] == start_i:
                cur_height = max(cur_height, 1 + get_height(j))
            else:
                # all remaining junctions to left cannot be enclosed
                break
        # junctions after `i`: subset if end less than or equal to end_i
        for j in range(i + 1, len(start)):
            if end[j] <= end_i:
                cur_height = max(cur_height, 1 + get_height(j))
            elif start[j] >= end_i:
                # all remaining junctions to right cannot be enclosed
                break
        # update height/updated
        height[i] = cur_height
        updated[i] = True
        return cur_height

    # apply get_height to every value of i
    for i in range(len(height)):
        get_height(i)

    return height
