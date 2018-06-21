from libc.stdio cimport printf, sprintf, fprintf, stderr
from libc.stdlib cimport malloc, free, abort, abs
from libc.string cimport strlen

from libcpp.string cimport string


cdef extern from "sqlite3/sqlite3.h":
    struct sqlite3
    int sqlite3_open(const char *filename, sqlite3 ** ppDb) nogil
    int sqlite3_exec(sqlite3*, const char *sql, int (*callback)(void*, int, char**, char**), void *,
                     char ** errmsg) nogil
    int sqlite3_close(sqlite3*) nogil
    void sqlite3_free(void*) nogil

# cdef extern from "<string.h>" nogil:
#     size_t strlen(const char *s)

cdef extern from "<string>" nogil:
    string to_string (int val);

cdef:
    char rc
    char *zErrMsg = <char *> 0

    char *begin_trans = "PRAGMA foreign_keys = ON;BEGIN TRANSACTION;"
    char *commit = "COMMIT;"
    char *exp_insert = "INSERT INTO experiment (name) VALUES ('%s');"
    char *gene_insert = "INSERT INTO gene " \
                        "('id','name','strand','chromosome') " \
                        "VALUES ('%s','%s','%s','%s');"
    char *exon_insert = "INSERT INTO exon " \
                        "(gene_id, start, end, annotated_start, annotated_end, annotated) " \
                        "VALUES ('%s',%d,%d,%d,%d,%d);"
    char *junc_insert = "INSERT INTO junction " \
                        "(gene_id, start, end, has_reads, annotated) " \
                        "VALUES ('%s',%d,%d,0,%d);"
    char *junc_reads_insert = "INSERT INTO junction_reads " \
                              "(reads, experiment_name, junction_gene_id, junction_start, junction_end) " \
                              "VALUES (%d,'%s','%s',%d,%d);" \
                              "UPDATE junction SET has_reads=1 " \
                              "WHERE gene_id='%s' AND start=%d AND 'end'=%d;"
    char *ir_insert = "INSERT INTO intron_retention " \
                      "(gene_id, start, end, has_reads, annotated) " \
                      "VALUES ('%s',%d,%d,0,%d);"
    char *ir_reads_insert = "INSERT INTO intron_retention_reads " \
                            "(reads, experiment_name, intron_retention_gene_id, intron_retention_start, intron_retention_end) " \
                            "VALUES (%d,'%s','%s',%d,%d); " \
                            "UPDATE intron_retention SET has_reads=1 " \
                            "WHERE gene_id='%s' AND start=%d AND 'end'=%d;"

cdef int int_len(int value) nogil:
    cdef int l = (not value) + (value < 0)

    value = abs(value)

    while value:
        l += 1
        value /= 10

    return l


cdef int callback(void *NotUsed, int argc, char ** argv, char ** azColName) nogil:
    cdef:
        int i
    for i in range(argc):
        printf("%s = %s\n", azColName[i], argv[i])
    printf("\n")
    return 0

cdef int exec_db(sqlite3 *db, string sql) nogil:
    rc = sqlite3_exec(db, sql.c_str(), callback, <void *> 0, &zErrMsg)
    if rc != 0:
        fprintf(stderr, "%s\n", zErrMsg)
        sqlite3_free(zErrMsg)
    return rc

cdef sqlite3 *open_db(string file_name) nogil:
    cdef:
        sqlite3 *db

    sqlite3_open(file_name.c_str(), &db)
    exec_db(db, begin_trans)

    return db

cdef void close_db(sqlite3 *db) nogil:
    exec_db(db, commit)
    sqlite3_close(db)

cdef void gene(sqlite3 *db, string id, string name, string strand, string chromosome) nogil:
    cdef:
        int arg_len
        char *sql
        int rm_chars_len

    arg_len = id.length() + name.length() + strand.length() + chromosome.length()
    rm_chars_len = 4 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(gene_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, gene_insert, id.c_str(), name.c_str(), strand.c_str(), chromosome.c_str())
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting gene:\n")
        fprintf(stderr, "\tid: %s\n", id.c_str())
        fprintf(stderr, "\tname: %s\n", name.c_str())
        fprintf(stderr, "\tstrand: %s\n", strand.c_str())
        fprintf(stderr, "\tchromosome: %s\n", chromosome.c_str())
        abort()

    free(sql)

cdef void exon(sqlite3 *db, string gene_id, int start, int end, int annotated_start, int annotated_end, int annotated) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + int_len(start) + int_len(end) + int_len(annotated_start) + int_len(annotated_end) + int_len(annotated)
    rm_chars_len = 6 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(exon_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, exon_insert, gene_id.c_str(), start, end, annotated_start, annotated_end, annotated)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting exon:\n")
        fprintf(stderr, "\tgene_id: %s\n", gene_id.c_str())
        fprintf(stderr, "\tstart: %d\n", start)
        fprintf(stderr, "\tend: %d\n", end)
        fprintf(stderr, "\tannotated start: %d\n", annotated_start)
        fprintf(stderr, "\tannotated end: %d\n", annotated_end)
        fprintf(stderr, "\tannotated: %d\n", annotated)
        abort()

    free(sql)



cdef void junction(sqlite3 *db, string gene_id, int start, int end, int annotated) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + int_len(start) + int_len(end) + int_len(annotated)
    rm_chars_len = 4 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(junc_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, junc_insert, gene_id.c_str(), start, end, annotated)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting junction:\n", )
        fprintf(stderr, "\tgene_id: %s\n", gene_id.c_str())
        fprintf(stderr, "\tstart: %d\n", start)
        fprintf(stderr, "\tend: %d\n", end)
        fprintf(stderr, "\tannotated: %d\n", annotated)
        abort()

    free(sql)

cdef void intron_retention(sqlite3 *db, string gene_id, int start, int end, int annotated) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + int_len(start) + int_len(end) + int_len(annotated)
    rm_chars_len = 4 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(ir_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, ir_insert, gene_id.c_str(), start, end, annotated)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting intron retention:\n", )
        fprintf(stderr, "\tgene_id: %s\n", gene_id.c_str())
        fprintf(stderr, "\tstart: %d\n", start)
        fprintf(stderr, "\tend: %d\n", end)
        fprintf(stderr, "\tannotated: %d\n", annotated)
        abort()

    free(sql)

cdef void experiment(sqlite3 *db, string name) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = name.length()
    rm_chars_len = 2
    sql = <char *> malloc(sizeof(char) * (strlen(exp_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, exp_insert, name.c_str())

    if exec_db(db, sql):
        fprintf(stderr, "Error inserting experiment:\n")
        fprintf(stderr, "\tname: %s\n", name.c_str())
        abort()

    free(sql)

cdef void junction_reads(sqlite3 *db, int reads, string exp_name, string junc_gene_id, int junc_start,
                         int junc_end) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = int_len(reads) + exp_name.length() + junc_gene_id.length() + int_len(junc_start) + int_len(
        junc_end) + junc_gene_id.length() + int_len(junc_start) + int_len(junc_end)
    rm_chars_len = 8 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(junc_reads_insert) + arg_len - rm_chars_len + 1))

    sprintf(sql, junc_reads_insert, reads, exp_name.c_str(), junc_gene_id.c_str(), junc_start, junc_end, junc_gene_id,
            junc_start,
            junc_end)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting junc reads:\n")
        fprintf(stderr, "\treads: %d\n", reads)
        fprintf(stderr, "\texperiment name: %s\n", exp_name.c_str())
        fprintf(stderr, "\tjunction gene id: %s\n", junc_gene_id.c_str())
        fprintf(stderr, "\tjunction start: %d\n", junc_start)
        fprintf(stderr, "\tjunction end: %d\n", junc_end)
        abort()

    free(sql)

cdef void intron_retention_reads(sqlite3 *db, int reads, string exp_name, string ir_gene_id, int ir_start,
                                 int ir_end) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = int_len(reads) + exp_name.length() + (ir_gene_id.length() * 2) + (int_len(ir_start) * 2) + (
            int_len(ir_end) * 2)
    rm_chars_len = 8 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(ir_reads_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, ir_reads_insert, reads, exp_name.c_str(), ir_gene_id.c_str(), ir_start, ir_end, ir_gene_id.c_str(),
            ir_start, ir_end)

    if exec_db(db, sql):
        fprintf(stderr, "Error inserting ir reads:\n")
        fprintf(stderr, "\treads: %d\n", reads)
        fprintf(stderr, "\texperiment name: %s\n", exp_name.c_str())
        fprintf(stderr, "\tjunction gene id: %s\n", ir_gene_id.c_str())
        fprintf(stderr, "\tjunction start: %d\n", ir_start)
        fprintf(stderr, "\tjunction end: %d\n", ir_end)
        abort()

    free(sql)
