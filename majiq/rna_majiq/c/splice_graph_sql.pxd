from libcpp.string cimport string
from libcpp.vector cimport vector
cdef extern from "sqlite3.h":
    struct sqlite3

cdef int open_db(string file_name, sqlite3 ** db) nogil
cdef int close_db(sqlite3 *db) nogil
cdef int gene(sqlite3 *db, string id, string name, string strand, string chromosome) nogil
cdef int gene_overlap(sqlite3 *db, string id, vector[string] overlapping) nogil
# cdef int experiment(sqlite3 *db, string name) nogil
cdef int exon(sqlite3 *db, string gene_id, int start, int end, int annotated_start, int annotated_end,
               bint annotated) nogil
cdef int transcript_exon(sqlite3 *db, string gene_id, string transcript_id, int start, int end) nogil
cdef int junction(sqlite3 *db, string gene_id, int start, int end, bint annotated, bint is_simplified,
                  bint is_constitutive, bint has_flag) nogil
cdef int junction_reads(sqlite3 *db, int reads, string exp_name, string junc_gene_id, int junc_start,
                         int junc_end) nogil
cdef int intron_retention(sqlite3 *db, string gene_id, int start, int end, bint annotated, bint is_simplified,
                          bint is_constitutive, bint has_flag) nogil
cdef int intron_retention_reads(sqlite3 *db, int reads, string exp_name, string ir_gene_id, int ir_start,
                                 int ir_end) nogil
cdef int alt_start(sqlite3 *db, string gene_id, int coordinate) nogil
cdef int alt_end(sqlite3 *db, string gene_id, int coordinate) nogil
