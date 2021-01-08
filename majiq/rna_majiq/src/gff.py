"""
gff.py

Functions for basic parsing of GFF3 format files
"""

from collections import namedtuple
import urllib.parse as urllib
import gzip

gffInfoFields = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)


def __parse_gff_attributes(attribute_string):
    """
    Parse the GFF3 attribute column and return a dict
    :param attribute_string:
    """
    if attribute_string == ".":
        return {}
    ret = {}
    for attribute in (x for x in attribute_string.split(";") if x):
        try:
            key, value = attribute.split("=")
        except ValueError as e:
            print(f"Unable to unpack attribute {attribute} not in GFF3 style")
            raise e
        key = urllib.unquote(key)
        if key in ret:
            key = "extra_%s" % key
            if key not in ret:
                ret[key] = []
            ret[key].append(urllib.unquote(value))
        else:
            ret[key] = urllib.unquote(value)
    return ret


def parse_gff3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    # Parse with transparent decompression
    open_func = gzip.open if filename.endswith(".gz") else open
    with open_func(filename, mode="rt") as infile:

        for line in infile:

            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            # If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            # Normalize data
            normalized_info = {
                "seqid": None if parts[0] == "." else urllib.unquote(parts[0]),
                "source": None if parts[1] == "." else urllib.unquote(parts[0]),
                "type": None if parts[2] == "." else urllib.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.unquote(parts[7]),
                "attributes": __parse_gff_attributes(parts[8]),
            }
            # Alternatively, you can emit the dictionary here, if you need mutabwility:
            # yield normalized_info
            yield GFFRecord(**normalized_info)
