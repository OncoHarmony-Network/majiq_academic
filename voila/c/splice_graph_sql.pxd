cdef extern from "sqlite3/sqlite3.h":
    struct sqlite3

cdef sqlite3 *open_db(char *file_name) nogil
cdef void close_db(sqlite3 *db) nogil
cdef void gene(sqlite3 *db, char *id, char *name, char *strand, char *chromosome) nogil
cdef void experiment(sqlite3 *db, char *name) nogil
cdef void exon(sqlite3 *db, char *gene_id, int start, int end, int annotated) nogil
cdef void junction(sqlite3 *db, char *gene_id, int start, int end, int annotated) nogil
cdef void junction_reads(sqlite3 *db, int reads, char *exp_name, char *junc_gene_id, int junc_start, int junc_end) nogil
cdef void intron_retention(sqlite3 *db, char *gene_id, int start, int end, int annotated) nogil
cdef void intron_retention_reads(sqlite3 *db, int reads, char *exp_name, char *ir_gene_id, int ir_start, int ir_end) nogil