#from majiq.grimoire.exon cimport Exon

cdef class Junction:
    cdef public int start
    cdef public int end
    cdef public str gene_id
    cdef public int index
    cdef public int lsv_index

    cdef public int nreads
    cdef public object acceptor
    cdef public object donor
    cdef public bint annot
    cdef public bint intronic
