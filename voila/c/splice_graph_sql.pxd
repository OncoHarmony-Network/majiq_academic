from libcpp.string cimport string
cdef extern from "sqlite3/sqlite3.h":
    struct sqlite3

cdef sqlite3 *open_db(string file_name) nogil
cdef void close_db(sqlite3 *db) nogil
cdef void gene(sqlite3 *db, string id, string name, string strand, string chromosome) nogil
cdef void experiment(sqlite3 *db, string name) nogil
cdef void exon(sqlite3 *db, string gene_id, int start, int end, int annotated_start, int annotated_end, bint annotated) nogil
cdef void junction(sqlite3 *db, string gene_id, int start, int end, bint annotated) nogil
cdef void junction_reads(sqlite3 *db, int reads, string exp_name, string junc_gene_id, int junc_start, int junc_end) nogil
cdef void intron_retention(sqlite3 *db, string gene_id, int start, int end, bint annotated) nogil
cdef void intron_retention_reads(sqlite3 *db, int reads, string exp_name, string ir_gene_id, int ir_start, int ir_end) nogil
cdef void alt_start(sqlite3 *db, string gene_id, int coordinate) nogil
cdef void alt_end(sqlite3 *db, string gene_id, int coordinate) nogil
