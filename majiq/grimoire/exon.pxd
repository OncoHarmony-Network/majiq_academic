cdef class Exon:
    cdef public int start
    cdef public int end
    cdef bint annot

    cdef public tuple db_coords
    cdef public set ib
    cdef public set ob
    cdef public bint intron
