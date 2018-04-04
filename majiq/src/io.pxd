
from libcpp.vector cimport vector
from majiq.src.internals.grimoire cimport  Gene

cdef int _read_gff(str filename, vector[Gene *] glist, object logging) except -1