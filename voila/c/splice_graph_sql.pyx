from libc.stdio cimport printf, sprintf, fprintf, stderr
from libc.stdlib cimport malloc, free, abort

cdef extern from "<sqlite3.h>":
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
                        "(gene_id, start, end, intron_retention, annotated) " \
                        "VALUES ('%s','%s','%s','%s','%s');"
    char *junc_insert = "INSERT INTO junction " \
                        "(gene_id, start, end, has_reads, intron_retention, annotated) " \
                        "VALUES ('%s','%s','%s','0','%s', '%s');"
    char *reads_insert = "INSERT INTO reads " \
                         "(reads, experiment_name, junction_gene_id, junction_start, junction_end) " \
                         "VALUES ('%s','%s','%s','%s','%s');" \
                         "UPDATE junction SET has_reads='1' " \
                         "WHERE gene_id='%s' AND start='%s' AND end='%s';"

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

cdef void exon(sqlite3 *db, char *gene_id, int start, int end, int intron_retention, int annotated) nogil:
    cdef:
        char *start_str
        char *end_str
        char *ir_str
        char *annot_str

        int start_len
        int end_len
        int ir_len
        int annot_len

        int arg_len
        int rm_chars_len
        char *sql

    start_len = int_len(start)
    end_len = int_len(end)
    ir_len = int_len(intron_retention)
    annot_len = int_len(annotated)

    start_str = <char *> malloc(sizeof(char) * (start_len + 1))
    end_str = <char *> malloc(sizeof(char) * (end_len + 1))
    ir_str = <char *> malloc(sizeof(char) * (ir_len + 1))
    annot_str = <char *> malloc(sizeof(char) * (annot_len + 1))

    sprintf(start_str, "%d", start)
    sprintf(end_str, "%d", end)
    sprintf(ir_str, "%d", intron_retention)
    sprintf(annot_str, "%d", annotated)

    arg_len = strlen(gene_id) + start_len + end_len + ir_len + annot_len
    rm_chars_len = 5 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(exon_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, exon_insert, gene_id, start_str, end_str, ir_str, annot_str)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting exon:\n")
        fprintf(stderr, "\tgene_id: %s\n", gene_id)
        fprintf(stderr, "\tstart: %d\n", start)
        fprintf(stderr, "\tend: %d\n", end)
        fprintf(stderr, "\tintron retention: %d\n", intron_retention)
        fprintf(stderr, "\tannotated: %d\n", annotated)
        abort()

    free(sql)
    free(start_str)
    free(end_str)
    free(ir_str)
    free(annot_str)

cdef void junction(sqlite3 *db, char *gene_id, int start, int end, int intron_retention, int annotated) nogil:
    cdef:
        char *start_str
        char *end_str
        char *ir_str
        char *annot_str

        int start_len
        int end_len
        int ir_len
        int annot_len

        int arg_len
        int rm_chars_len
        char *sql

    start_len = int_len(start)
    end_len = int_len(end)
    ir_len = int_len(intron_retention)
    annot_len = int_len(annotated)

    start_str = <char *> malloc(sizeof(char) * (start_len + 1))
    end_str = <char *> malloc(sizeof(char) * (end_len + 1))
    ir_str = <char *> malloc(sizeof(char) * (ir_len + 1))
    annot_str = <char *> malloc(sizeof(char) * (annot_len + 1))

    sprintf(start_str, "%d", start)
    sprintf(end_str, "%d", end)
    sprintf(ir_str, "%d", intron_retention)
    sprintf(annot_str, "%d", annotated)

    arg_len = strlen(gene_id) + start_len + end_len + ir_len + annot_len
    rm_chars_len = 5 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(junc_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, junc_insert, gene_id, start_str, end_str, ir_str, annot_str)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting junction:\n",)
        fprintf(stderr, "\tgene_id: %s\n", gene_id)
        fprintf(stderr, "\tstart: %d\n", start)
        fprintf(stderr, "\tend: %d\n", end)
        fprintf(stderr, "\tintron retention: %d\n", intron_retention)
        fprintf(stderr, "\tannotated: %d\n", annotated)
        abort()

    free(sql)
    free(start_str)
    free(end_str)
    free(ir_str)
    free(annot_str)

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

cdef void reads(sqlite3 *db, int reads, char *exp_name, char *junc_gene_id, int junc_start, int junc_end) nogil:
    cdef:
        char *reads_str
        char *start_str
        char *end_str

        int reads_len
        int start_len
        int end_len

        int arg_len
        int rm_chars_len
        char *sql

    reads_len = int_len(reads)
    start_len = int_len(junc_start)
    end_len = int_len(junc_end)

    reads_str = <char *> malloc(sizeof(char) * (reads_len + 1))
    start_str = <char *> malloc(sizeof(char) * (start_len + 1))
    end_str = <char *> malloc(sizeof(char) * (end_len + 1))

    sprintf(reads_str, "%d", reads)
    sprintf(start_str, "%d", junc_start)
    sprintf(end_str, "%d", junc_end)

    arg_len = reads_len + strlen(exp_name) + strlen(junc_gene_id) + start_len + end_len
    rm_chars_len = 5 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(reads_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, reads_insert, reads_str, exp_name, junc_gene_id, start_str, end_str, junc_gene_id, start_str, end_str)
    if exec_db(db, sql):
        fprintf(stderr, "Error inserting reads:\n")
        fprintf(stderr, "\treads: %d\n", reads)
        fprintf(stderr, "\texperiment name: %s\n", exp_name)
        fprintf(stderr, "\tjunction gene id: %s\n", junc_gene_id)
        fprintf(stderr, "\tjunction start: %d\n", junc_start)
        fprintf(stderr, "\tjunction end: %d\n", junc_end)
        abort()

    free(start_str)
    free(end_str)
    free(sql)
