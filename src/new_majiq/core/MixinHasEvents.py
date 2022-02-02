"""
MixinHasEvents.py

Mixin class for objects that hold on to a dataset referencing
splicegraph events that can be loaded using GeneIntrons and GeneJunctions

Author: Joseph K Aicher
"""

from pathlib import Path
from typing import Final, Optional, Union

import xarray as xr

import new_majiq.constants as constants

from .Events import Events, _Events
from .GeneIntrons import GeneIntrons
from .GeneJunctions import GeneJunctions


class MixinHasEvents(object):
    """Mixin for dataset defining associated Events

    Mixin for dataset defining associated Events. Several dataset classes are
    defined over e_idx or ec_idx for events and need to keep information about
    the Events that they are related to. However, they frequently do not need
    to have them loaded except in rare occasions. This tracks the bare-minimum
    information needed and enables them to load the full Events object as
    needed.

    Parameters
    ----------
    df: Optional[xr.Dataset]
        If not None, passed-through dataset to validate sizes against
        (i.e.  ec_idx, e_idx)
    events_df: xr.Dataset
        dataset that can be loaded along with matching introns/junctions as
        Events
    """

    def __init__(self, df: Optional[xr.Dataset], events_df: xr.Dataset):
        if df is not None:
            try:
                df_ec_size = df.sizes["ec_idx"]
            except KeyError:
                pass
            else:
                if df_ec_size != events_df.sizes["ec_idx"]:
                    raise ValueError("df/events ec_idx are not same size")
            try:
                df_e_size = df.sizes["e_idx"]
            except KeyError:
                pass
            else:
                if df_e_size != events_df.sizes["e_idx"]:
                    raise ValueError("df/events e_idx are not same size")

        self.events_df: Final[xr.Dataset] = events_df
        return

    @property
    def num_connections(self) -> int:
        """Total number of connections over all events"""
        return self.events_df.sizes["ec_idx"]

    @property
    def num_events(self) -> int:
        """Total number of events"""
        return self.events_df.sizes["e_idx"]

    def events_to_zarr(
        self, path: Union[str, Path], mode: str, consolidated: bool = True
    ) -> None:
        """Save events information to specified path"""
        self.events_df.chunk(self.events_df.sizes).to_zarr(
            path, mode=mode, group=constants.NC_EVENTS, consolidated=consolidated
        )
        return

    def get_events(
        self,
        introns: GeneIntrons,
        junctions: GeneJunctions,
    ) -> Events:
        """Construct :py:class:`Events` using saved dataset and introns, junctions

        Parameters
        ----------
        introns: GeneIntrons
        junctions: GeneJunctions

        Returns
        -------
        Events
        """
        if self.events_df.intron_hash != introns.checksum():
            raise ValueError("GeneIntrons checksums do not match")
        if self.events_df.junction_hash != junctions.checksum():
            raise ValueError("GeneJunctions checksums do not match")
        return Events(
            _Events(
                introns._gene_introns,
                junctions._gene_junctions,
                self.events_df.ref_exon_idx,
                self.events_df.event_type,
                self.events_df._offsets,
                self.events_df.is_intron,
                self.events_df.connection_idx,
            )
        )
