from libc.stdio cimport printf, sprintf, fprintf, stderr
from libc.stdlib cimport malloc, free, abort

cdef extern from "sqlite3/sqlite3.h":
    struct sqlite3
    int sqlite3_open(const char *filename, sqlite3 ** ppDb) nogil
    int sqlite3_exec(sqlite3*, const char *sql, int (*callback)(void*, int, char**, char**), void *,
                     char ** errmsg) nogil
    int sqlite3_close(sqlite3*) nogil
    void sqlite3_free(void*)nogil

cdef extern from "<string.h>" nogil:
    size_t strlen(const char *s)

cdef:
    char rc
    char *zErrMsg = <char *> 0

    char *begin_trans = "BEGIN TRANSACTION;"
    char *commit = "COMMIT;"
    char *exp_insert = "INSERT INTO experiment (name) VALUES ('%s');"
    char *gene_insert = "INSERT INTO gene " \
                        "('id','name','strand','chromosome') " \
                        "VALUES ('%s','%s','%s','%s');"
    char *exon_insert = "INSERT INTO exon " \
                        "(gene_id, start, end, annotated) " \
                        "VALUES ('%s',%d,%d,%d);"
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
    cdef:
        int l = not value

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

cdef int exec_db(sqlite3 *db, char *sql) nogil:
    rc = sqlite3_exec(db, sql, callback, <void *> 0, &zErrMsg)
    if rc != 0:
        fprintf(stderr, "%s\n", zErrMsg)
        sqlite3_free(zErrMsg)
    return rc

cdef sqlite3 *open_db(char *file_name) nogil:
    cdef:
        sqlite3 *db

    sqlite3_open(file_name, &db)
    exec_db(db, "PRAGMA foreign_keys = ON;")
    exec_db(db, begin_trans)

    return db

cdef void close_db(sqlite3 *db) nogil:
    exec_db(db, commit)
    sqlite3_close(db)

cdef void gene(sqlite3 *db, char *id, char *name, char *strand, char *chromosome) nogil:
    cdef:
        int arg_len
        char *sql
        int rm_chars_len

    arg_len = strlen(id) + strlen(name) + strlen(strand) + strlen(chromosome)
    rm_chars_len = 4 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(gene_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, gene_insert, id, name, strand, chromosome)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting gene:\n")
        fprintf(stderr, "\tid: %s\n", id)
        fprintf(stderr, "\tname: %s\n", name)
        fprintf(stderr, "\tstrand: %s\n", strand)
        fprintf(stderr, "\tchromosome: %s\n", chromosome)
        abort()

    free(sql)

cdef void exon(sqlite3 *db, char *gene_id, int start, int end, int annotated) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = strlen(gene_id) + int_len(start) + int_len(end) + int_len(annotated)
    rm_chars_len = 4 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(exon_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, exon_insert, gene_id, start, end, annotated)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting exon:\n")
        fprintf(stderr, "\tgene_id: %s\n", gene_id)
        fprintf(stderr, "\tstart: %d\n", start)
        fprintf(stderr, "\tend: %d\n", end)
        fprintf(stderr, "\tannotated: %d\n", annotated)
        abort()

    free(sql)

cdef void junction(sqlite3 *db, char *gene_id, int start, int end, int annotated) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = strlen(gene_id) + int_len(start) + int_len(end) + int_len(annotated)
    rm_chars_len = 4 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(junc_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, junc_insert, gene_id, start, end, annotated)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting junction:\n", )
        fprintf(stderr, "\tgene_id: %s\n", gene_id)
        fprintf(stderr, "\tstart: %d\n", start)
        fprintf(stderr, "\tend: %d\n", end)
        fprintf(stderr, "\tannotated: %d\n", annotated)
        abort()

    free(sql)

cdef void intron_retention(sqlite3 *db, char *gene_id, int start, int end, int annotated) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = strlen(gene_id) + int_len(start)  + int_len(end) + int_len(annotated)
    rm_chars_len = 4 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(ir_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, ir_insert, gene_id, start, end, annotated)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting intron retention:\n", )
        fprintf(stderr, "\tgene_id: %s\n", gene_id)
        fprintf(stderr, "\tstart: %d\n", start)
        fprintf(stderr, "\tend: %d\n", end)
        fprintf(stderr, "\tannotated: %d\n", annotated)
        abort()

    free(sql)

cdef void experiment(sqlite3 *db, char *name) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = strlen(name)
    rm_chars_len = 2
    sql = <char *> malloc(sizeof(char) * (strlen(exp_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, exp_insert, name)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting experiment:\n")
        fprintf(stderr, "\tname: %s\n", name)
        abort()

    free(sql)

cdef void junction_reads(sqlite3 *db, int reads, char *exp_name, char *junc_gene_id, int junc_start,
                         int junc_end) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = int_len(reads) + strlen(exp_name) + strlen(junc_gene_id) + int_len(junc_start) + int_len(junc_end) + strlen(
        junc_gene_id) + int_len(junc_start) + int_len(junc_end)
    rm_chars_len = 8 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(junc_reads_insert) + arg_len - rm_chars_len + 1))

    sprintf(sql, junc_reads_insert, reads, exp_name, junc_gene_id, junc_start, junc_end, junc_gene_id, junc_start, junc_end)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting junc reads:\n")
        fprintf(stderr, "\treads: %d\n", reads)
        fprintf(stderr, "\texperiment name: %s\n", exp_name)
        fprintf(stderr, "\tjunction gene id: %s\n", junc_gene_id)
        fprintf(stderr, "\tjunction start: %d\n", junc_start)
        fprintf(stderr, "\tjunction end: %d\n", junc_end)
        abort()

    free(sql)

cdef void intron_retention_reads(sqlite3 *db, int reads, char *exp_name, char *ir_gene_id, int ir_start,
                                 int ir_end) nogil:
    cdef:
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = int_len(reads) + strlen(exp_name) + (strlen(ir_gene_id) * 2) + (int_len(ir_start) * 2) + (
                int_len(ir_end) * 2)
    rm_chars_len = 8 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(ir_reads_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, ir_reads_insert, reads, exp_name, ir_gene_id, ir_start, ir_end, ir_gene_id, ir_start, ir_end)

    if exec_db(db, sql):
        fprintf(stderr, "Error inserting ir reads:\n")
        fprintf(stderr, "\treads: %d\n", reads)
        fprintf(stderr, "\texperiment name: %s\n", exp_name)
        fprintf(stderr, "\tjunction gene id: %s\n", ir_gene_id)
        fprintf(stderr, "\tjunction start: %d\n", ir_start)
        fprintf(stderr, "\tjunction end: %d\n", ir_end)
        abort()

    free(sql)
