from libc.stdio cimport sprintf, fprintf, stderr
from libc.stdlib cimport malloc, free, abort, abs
from libc.string cimport strlen
from libcpp.string cimport string

cdef extern from "sqlite3.h":
    int SQLITE_BUSY
    struct sqlite3
    int sqlite3_open(const char *, sqlite3 **) nogil
    int sqlite3_exec(sqlite3*, const char *, int (*)(void*, int, char**, char**), void *, char **) nogil
    int sqlite3_close(sqlite3*) nogil
    void sqlite3_free(void*) nogil
    int sqlite3_busy_timeout(sqlite3*, int) nogil

cdef:
    int rc
    char *zErrMsg = <char *> 0

    char *begin_trans = "PRAGMA foreign_keys = ON;BEGIN DEFERRED TRANSACTION;"
    char *commit = "COMMIT;"
    char *exp_insert = "INSERT INTO experiment (name) VALUES ('%s');"
    char *gene_insert = "INSERT INTO gene " \
                        "(id,name,strand,chromosome) " \
                        "VALUES ('%s','%s','%s','%s');"
    char *alt_start_insert = "INSERT INTO alt_start (gene_id,coordinate) VALUES ('%s',%d);"
    char *alt_end_insert = "INSERT INTO alt_end (gene_id,coordinate) VALUES ('%s',%d);"
    char *exon_insert = "INSERT INTO exon " \
                        "(gene_id,start,end,annotated_start,annotated_end,annotated) " \
                        "VALUES ('%s',%d,%d,%d,%d,%d);"
    char *junc_insert = "INSERT INTO junction " \
                        "(gene_id,start,end,has_reads,annotated) " \
                        "VALUES ('%s',%d,%d,0,%d);"
    char *junc_reads_insert = "INSERT INTO junction_reads " \
                              "(reads,experiment_name,junction_gene_id,junction_start,junction_end) " \
                              "VALUES (%d,'%s','%s',%d,%d);" \
                              "UPDATE junction SET has_reads=1 " \
                              "WHERE gene_id='%s' AND start=%d AND end=%d;"
    char *ir_insert = "INSERT INTO intron_retention " \
                      "(gene_id,start,end,has_reads,annotated) " \
                      "VALUES ('%s',%d,%d,0,%d);"
    char *ir_reads_insert = "INSERT INTO intron_retention_reads " \
                            "(reads,experiment_name,intron_retention_gene_id,intron_retention_start,intron_retention_end) " \
                            "VALUES (%d,'%s','%s',%d,%d); " \
                            "UPDATE intron_retention SET has_reads=1 " \
                            "WHERE gene_id='%s' AND start=%d AND end=%d;"

cdef int int_len(int value) nogil:
    cdef int l = (not value) + (value < 0)

    value = abs(value)

    while value:
        l += 1
        value /= 10

    return l

cdef int callback(void *NotUsed, int argc, char ** argv, char ** azColName) nogil:
    cdef int i

    for i in range(argc):
        fprintf(stderr, "%s = %s\n", azColName[i], argv[i])
    fprintf(stderr, "\n")
    return 0

cdef int exec_db(sqlite3 *db, char *sql) nogil:
    while True:
        rc = sqlite3_exec(db, sql, callback, <void *> 0, &zErrMsg)
        if rc != SQLITE_BUSY:
            break

    if rc:
        fprintf(stderr, "exec_db: %s: %d\n", zErrMsg, rc)
        fprintf(stderr, '%s\n', sql)
        sqlite3_free(zErrMsg)
        abort()

    return rc

cdef sqlite3 *open_db(string file_name) nogil:
    cdef sqlite3 *db

    rc = sqlite3_open(file_name.c_str(), &db)
    if rc:
        fprintf(stderr, 'open_db: %s: %d\n', file_name.c_str(), rc)
        abort()

    rc = sqlite3_busy_timeout(db, 120 * 1000)
    if rc:
        fprintf(stderr, 'busy_timeout: %s: %d\n', file_name.c_str(), rc)
        abort()

    exec_db(db, begin_trans)

    return db

cdef void close_db(sqlite3 *db) nogil:
    exec_db(db, commit)
    rc = sqlite3_close(db)
    if rc:
        fprintf(stderr, 'close: %d\n', rc)

cdef void experiment(sqlite3 *db, string name) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = name.length()
    rm_chars_len = 2
    sql = <char *> malloc(sizeof(char) * (strlen(exp_insert) + arg_len - rm_chars_len + 1))
    exec_db(db, sql)
    free(sql)

cdef void gene(sqlite3 *db, string id, string name, string strand, string chromosome) nogil:
    cdef:
        int arg_len
        char *sql
        int rm_chars_len

    arg_len = id.length() + name.length() + strand.length() + chromosome.length()
    rm_chars_len = 4 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(gene_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, gene_insert, id.c_str(), name.c_str(), strand.c_str(), chromosome.c_str())

    exec_db(db, sql)
    free(sql)

cdef void alt_gene(sqlite3 *db, string gene_id, int coordinate, char *sql_string) nogil:
    cdef:
        int arg_len
        char *sql
        int rm_chars_len

    arg_len = gene_id.length() + int_len(coordinate)
    rm_chars_len = 2 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(sql_string) + arg_len - rm_chars_len + 1))
    sprintf(sql, sql_string, gene_id.c_str(), coordinate)
    exec_db(db, sql)
    free(sql)

cdef void alt_start(sqlite3 *db, string gene_id, int coordinate) nogil:
    alt_gene(db, gene_id, coordinate, alt_start_insert)

cdef void alt_end(sqlite3 *db, string gene_id, int coordinate) nogil:
    alt_gene(db, gene_id, coordinate, alt_end_insert)

cdef void exon(sqlite3 *db, string gene_id, int start, int end, int annotated_start, int annotated_end,
               bint annotated) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + int_len(start) + int_len(end) + int_len(annotated_start) + int_len(
        annotated_end) + int_len(annotated)
    rm_chars_len = 6 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(exon_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, exon_insert, gene_id.c_str(), start, end, annotated_start, annotated_end, annotated)
    exec_db(db, sql)
    free(sql)

cdef void inter_exon(sqlite3 *db, string gene_id, int start, int end, bint annotated, char *sql_string) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + int_len(start) + int_len(end) + int_len(annotated)
    rm_chars_len = 4 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(sql_string) + arg_len - rm_chars_len + 1))
    sprintf(sql, sql_string, gene_id.c_str(), start, end, annotated)
    exec_db(db, sql)
    free(sql)

cdef void junction(sqlite3 *db, string gene_id, int start, int end, bint annotated) nogil:
    inter_exon(db, gene_id, start, end, annotated, junc_insert)

cdef void intron_retention(sqlite3 *db, string gene_id, int start, int end, bint annotated) nogil:
    inter_exon(db, gene_id, start, end, annotated, ir_insert)

cdef void reads_update(sqlite3 *db, int reads, string exp_name, string gene_id, int start, int end,
                       const char *sql_string) nogil:
    if not reads:
        return

    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = int_len(reads) + exp_name.length() + ((gene_id.length() + int_len(start) + int_len(end)) * 2)
    rm_chars_len = 8 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(sql_string) + arg_len - rm_chars_len + 1))
    sprintf(sql, sql_string, reads, exp_name.c_str(), gene_id.c_str(), start, end, gene_id.c_str(), start, end)
    exec_db(db, sql)
    free(sql)

cdef void junction_reads(sqlite3 *db, int reads, string exp_name, string junc_gene_id, int junc_start,
                         int junc_end) nogil:
    reads_update(db, reads, exp_name, junc_gene_id, junc_start, junc_end, junc_reads_insert)

cdef void intron_retention_reads(sqlite3 *db, int reads, string exp_name, string ir_gene_id, int ir_start,
                                 int ir_end) nogil:
    reads_update(db, reads, exp_name, ir_gene_id, ir_start, ir_end, ir_reads_insert)
