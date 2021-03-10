"""
version.py

Access string representation of new_majiq version
"""


def version() -> str:
    from new_majiq.internals import __version__
    return __version__
