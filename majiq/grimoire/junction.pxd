#from majiq.grimoire.exon cimport Exon

cdef class Junction:
    cdef public int start
    cdef public int end
    cdef str gene_id
    cdef public int index

    cdef int nreads
    cdef public object acceptor
    cdef public object donor
    cdef public bint annot
    cdef public bint intronic
