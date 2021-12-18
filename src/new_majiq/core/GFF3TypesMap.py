"""
GFF3Types.py

Specify how GFF3 types are used to build transcript models that are translated
into :class:`SpliceGraph`

Author: Joseph K Aicher
"""

from typing import Dict, List, Optional, Union

from new_majiq.internals import GFF3FeatureType, _default_gff3_types


class GFF3TypesMap(object):
    """Define how GFF3 types (column 3) parsed in gene/transcript/exon hierarchy

    Parameters
    ----------
    gff3_types_map: Optional[Dict[str, Union[str, GFF3FeatureType]]]
        Mapping from GFF3 type (column 3) to how it is parsed in the
        gene/transcript/exon hierarchy.
        Valid values are specified by keys or values of
        :attr:`GFF3FeatureType.__members__`
        If None, use default value from
        :func:`new_majiq.internals._default_gff3_types`.

    Notes
    -----
    MAJIQ identifies exons and matches them to their ancestors.
    The immediate parent of an exon is identified as its "transcript".
    If it is not found in the map as either ACCEPT_GENE or ACCEPT_TRANSCRIPT,
    it is not added to the splicegraph. In this case, if it is not explicitly
    rejected (REJECT_SILENT), the type of the parent is accounted for as a
    potential error in parsing (see `log_function` in
    :meth:`SpliceGraph.from_gff3`).
    The first ancestor that is found in the map as ACCEPT_GENE is considered
    to be the exon/transcript's gene. If no such ancestor is found, the type
    of the last (top-level) ancestor is accounted for as a potential error in
    parsing (as done for transcripts).

    It is possible to use HARD_SKIP to ignore records altogether from the
    hierarchy that have no chance of having exons as a child. This could
    potentially improve performance (smaller) but will raise an exception if an
    exon would have the skipped feature as a parent.

    REJECT_OTHER is one of the enumeration values of :class:`GFF3FeatureType`,
    but it can be generally be excluded. If a type is not found in the map,
    REJECT_OTHER is used as the default value.
    """

    def __init__(
        self, gff3_types_map: Optional[Dict[str, Union[str, GFF3FeatureType]]] = None
    ) -> None:
        self.current_map: Dict[str, GFF3FeatureType]
        if gff3_types_map is None:
            self.current_map = _default_gff3_types()
        else:
            self.current_map = {}
            for k, v in gff3_types_map.items():
                self[k] = v
        return

    def __repr__(self) -> str:
        # current_map but to names
        str_map = {k: v.name for k, v in self.current_map.items()}
        return f"GFF3TypesMap({str_map})"

    def __str__(self) -> str:
        # reversed map from how it's processed to lists of GFF3 types
        reversed_map: Dict[str, List[str]] = {}
        # fill reversed_map
        for k, v in self.current_map.items():
            try:
                reversed_map[v.name].append(k)
            except KeyError:
                reversed_map[v.name] = [k]
        # put it in sorted order
        reversed_str = ", ".join(
            f"{rk}<-{sorted(reversed_map[rk])}" for rk in sorted(reversed_map)
        )
        return f"GFF3TypesMap({reversed_str})"

    def __getitem__(self, k):
        return self.current_map[k]

    def __setitem__(self, k, feature_type):
        """Add or modify mapping from GFF3 type to how its treated in hierarchy

        Parameters
        ----------
        k: str
            Key for GFF3 type to set
        feature_type: Union[str, GFF3FeatureType]
            How the GFF3 type will be processed. If string, mapped to
            appropriate value of GFF3FeatureType using
            :attr:`GFF3FeatureType.__members__` (not case-sensitive).
            Otherwise, treated as GFF3FeatureType after confirming that it is a
            valid enumeration value
        """
        if not isinstance(k, str):
            raise ValueError(f"Keys in GFF3TypesMap must be string (key = {k})")
        try:
            if isinstance(feature_type, str):
                feature_type = GFF3FeatureType.__members__[feature_type.upper()]
            else:
                # if feature_type is not int or GFF3FeatureType, TypeError
                # but it accepts non-enumerated ints, so further checks name
                # against __members__, yielding either enumerated
                # GFF3FeatureType or KeyError
                feature_type = GFF3FeatureType.__members__[
                    GFF3FeatureType(feature_type).name
                ]
        except (TypeError, KeyError):
            raise ValueError(
                f"Invalid parsing option {feature_type} for GFF3 type {k}"
                f" (valid choices: {sorted(GFF3FeatureType.__members__)})"
            )
        self.current_map[k] = feature_type
        return

    def exon_types(self) -> List[str]:
        """List of values for GFF3 type (column 3) recognized as EXON"""
        return sorted(
            k for k, v in self.current_map.items() if v == GFF3FeatureType.EXON
        )

    def transcript_types(self) -> List[str]:
        """List of values for GFF3 type (column 3) recognized as transcript

        List of values for GFF3 type (column 3) recognized as transcript when
        the immediate parent of accepted exons. Note this includes both
        ACCEPT_TRANSCRIPT and ACCEPT_GENE.
        """
        return sorted(
            k
            for k, v in self.current_map.items()
            if v in (GFF3FeatureType.ACCEPT_GENE, GFF3FeatureType.ACCEPT_TRANSCRIPT)
        )

    def gene_types(self) -> List[str]:
        """List of values for GFF3 type (column 3) recognized as gene (ACCEPT_GENE)

        List of values for GFF3 type (column 3) recognized as gene when
        ancestor of exon (if multiple ancestors are recognized, the most recent
        ancestor in hierarchy is the assigned gene)
        """
        return sorted(
            k for k, v in self.current_map.items() if v == GFF3FeatureType.ACCEPT_GENE
        )

    def silent_types(self) -> List[str]:
        """List of values for GFF3 type (column 3) that will be silently ignored

        List of values for GFF3 type (column 3) that will be ignored, but
        silently. Non-accepted features that are not included here will be
        accounted for.
        """
        return sorted(
            k for k, v in self.current_map.items() if v == GFF3FeatureType.REJECT_SILENT
        )

    def hard_skip_types(self) -> List[str]:
        """List of values for GFF3 type (column 3) hard excluded from hierarchy

        List of values for GFF3 type (column 3) hard excluded from hierarchy.
        If an exon has one of these as an ancestor, an exception will be raised.
        """
        return sorted(
            k for k, v in self.current_map.items() if v == GFF3FeatureType.HARD_SKIP
        )
