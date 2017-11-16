from majiq.grimoire.junction cimport Junction
cimport numpy as np

cdef class Exon:
    cdef public int start
    cdef public int end
    cdef public bint annot

    cdef public tuple db_coords
    cdef public set ib
    cdef public set ob
    cdef public bint intron
    cdef int db_idx

cdef class Intron:
    cdef public int nchunks
    cdef public int start
    cdef public int end
    cdef public int chunk_len
    cdef public bint annot

    cdef public Junction junc1
    cdef public Junction junc2
    cdef public list junc1_cov
    cdef public list junc2_cov
    cdef public np.ndarray parts
    cdef int db_idx


#(156720572, 156721087, True), (156721031, 156721087, True),