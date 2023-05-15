from libc.stdio cimport sprintf, fprintf, stderr, stdout
from libc.stdlib cimport malloc, free, abs
from libc.string cimport strlen
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "sqlite3.h":
    int SQLITE_BUSY
    int SQLITE_LOCKED
    struct sqlite3
    int sqlite3_open(const char *, sqlite3 **) nogil
    int sqlite3_exec(sqlite3*, const char *, int (*)(void*, int, char**, char**), void *, char **) nogil
    int sqlite3_close(sqlite3*) nogil
    void sqlite3_free(void*) nogil
    int sqlite3_busy_timeout(sqlite3*, int) nogil

cdef:
    char *begin_trans = "PRAGMA foreign_keys = ON;BEGIN DEFERRED TRANSACTION;"
    char *commit = "COMMIT;"
    char *exp_insert = "INSERT INTO experiment (name) VALUES ('%s');"
    char *gene_insert = "INSERT INTO gene " \
                        "(id,name,strand,chromosome) " \
                        "VALUES ('%s','%s','%s','%s');"
    char *gene_overlap_insert = "INSERT INTO gene_overlap " \
                                "(gene_id_1,gene_id_2) " \
                                "VALUES ('%s','%s');"
    char *alt_start_insert = "INSERT INTO alt_start (gene_id,coordinate) VALUES ('%s',%d);"
    char *alt_end_insert = "INSERT INTO alt_end (gene_id,coordinate) VALUES ('%s',%d);"
    char *exon_insert = "INSERT INTO exon " \
                        "(gene_id,start,end,annotated_start,annotated_end,annotated) " \
                        "VALUES ('%s',%d,%d,%d,%d,%d);"
    char *transcript_exon_insert = "INSERT INTO transcript_exon " \
                        "(gene_id,transcript_id,start,end) " \
                        "VALUES ('%s','%s',%d,%d);"
    char *junc_insert = "INSERT INTO junction " \
                        "(gene_id,start,end,has_reads,annotated,is_simplified,is_constitutive,has_flag) " \
                        "VALUES ('%s',%d,%d,0,%d,%d,%d,%d);"
    char *junc_reads_insert = "INSERT INTO junction_reads " \
                              "(reads,experiment_name,junction_gene_id,junction_start,junction_end) " \
                              "VALUES (%d,'%s','%s',%d,%d);" \
                              "UPDATE junction SET has_reads=1 " \
                              "WHERE gene_id='%s' AND start=%d AND end=%d;"
    char *ir_insert = "INSERT INTO intron_retention " \
                      "(gene_id,start,end,has_reads,annotated,is_simplified,is_constitutive,has_flag) " \
                      "VALUES ('%s',%d,%d,0,%d,%d,%d,%d);"
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
    cdef char *zErrMsg
    rc = SQLITE_BUSY

    while rc == SQLITE_BUSY or rc == SQLITE_LOCKED:
        rc = sqlite3_exec(db, sql, callback, <int *> 0, &zErrMsg)

    if rc:
        with gil:
            raise Exception('exec_db: ' + zErrMsg.decode('utf-8') + '\n' + sql.decode('utf-8'))

    if zErrMsg:
        sqlite3_free(zErrMsg)

    return rc

cdef int open_db(string file_name, sqlite3 ** db) nogil:
    cdef int rc

    rc = sqlite3_open(file_name.c_str(), db)
    if rc:
        return rc
    return exec_db(db[0], begin_trans)

cdef int close_db(sqlite3 *db) nogil:
    cdef int rc
    rc = exec_db(db, commit)
    if rc:
        return rc
    return sqlite3_close(db)

cdef int gene(sqlite3 *db, string id, string name, string strand, string chromosome) nogil:
    cdef:
        int rc
        int arg_len
        char *sql
        int rm_chars_len

    arg_len = id.length() + name.length() + strand.length() + chromosome.length()
    rm_chars_len = 4 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(gene_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, gene_insert, id.c_str(), name.c_str(), strand.c_str(), chromosome.c_str())

    rc = exec_db(db, sql)
    free(sql)
    return rc

cdef int gene_overlap(sqlite3 *db, string id, vector[string] overlapping) nogil:
    cdef:
        int rc
        int arg_len
        char *sql
        int rm_chars_len

    arg_len = id.length() + overlapping[0].length()
    rm_chars_len = 2 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(gene_overlap_insert) + arg_len - rm_chars_len + 1) * overlapping.size())
    cdef char *pos = sql

    for i in range(overlapping.size()):
        pos += sprintf(pos, gene_overlap_insert, id.c_str(), overlapping[i].c_str())

    rc = exec_db(db, sql)
    free(sql)
    return rc

cdef int alt_gene(sqlite3 *db, string gene_id, int coordinate, char *sql_string) nogil:
    cdef:
        int rc
        int arg_len
        char *sql
        int rm_chars_len

    arg_len = gene_id.length() + int_len(coordinate)
    rm_chars_len = 2 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(sql_string) + arg_len - rm_chars_len + 1))
    sprintf(sql, sql_string, gene_id.c_str(), coordinate)

    rc = exec_db(db, sql)
    free(sql)
    return rc

cdef int alt_start(sqlite3 *db, string gene_id, int coordinate) nogil:
    return alt_gene(db, gene_id, coordinate, alt_start_insert)

cdef int alt_end(sqlite3 *db, string gene_id, int coordinate) nogil:
    return alt_gene(db, gene_id, coordinate, alt_end_insert)

cdef int exon(sqlite3 *db, string gene_id, int start, int end, int annotated_start, int annotated_end,
              bint annotated) nogil:
    cdef:
        int rc
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + int_len(start) + int_len(end) + int_len(annotated_start) + int_len(
        annotated_end) + int_len(annotated)
    rm_chars_len = 6 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(exon_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, exon_insert, gene_id.c_str(), start, end, annotated_start, annotated_end, annotated)

    rc = exec_db(db, sql)
    free(sql)
    return rc

cdef int transcript_exon(sqlite3 *db, string gene_id, string transcript_id, int start, int end) nogil:

    cdef:
        int rc
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + transcript_id.length() + int_len(start) + int_len(end)
    rm_chars_len = 4 * 2

    sql = <char *> malloc(sizeof(char) * (strlen(transcript_exon_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, transcript_exon_insert, gene_id.c_str(), transcript_id.c_str(), start, end)

    rc = exec_db(db, sql)
    free(sql)
    return rc


cdef int junction(sqlite3 *db, string gene_id, int start, int end, bint annotated, bint is_simplified, bint is_constitutive, bint has_flag) nogil:
    cdef:
        int rc
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + int_len(start) + int_len(end) + int_len(annotated) + int_len(is_simplified) + int_len(is_constitutive) + int_len(has_flag)
    rm_chars_len = 7 * 2  # how many characters replaced by args (%s, %d, etc.)

    sql = <char *> malloc(sizeof(char) * (strlen(junc_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, junc_insert, gene_id.c_str(), start, end, annotated, is_simplified, is_constitutive, has_flag)

    rc = exec_db(db, sql)
    free(sql)
    return rc

cdef int intron_retention(sqlite3 *db, string gene_id, int start, int end, bint annotated, bint is_simplified, bint is_constitutive, bint has_flag) nogil:
    cdef:
        int rc
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = gene_id.length() + int_len(start) + int_len(end) + int_len(annotated) + int_len(is_simplified) + int_len(is_constitutive) + int_len(has_flag)
    rm_chars_len = 7 * 2  # how many characters replaced by args (%s, %d, etc.)

    sql = <char *> malloc(sizeof(char) * (strlen(ir_insert) + arg_len - rm_chars_len + 1))
    sprintf(sql, ir_insert, gene_id.c_str(), start, end, annotated, is_simplified, is_constitutive, has_flag)

    rc = exec_db(db, sql)
    free(sql)
    return rc

cdef int reads_update(sqlite3 *db, int reads, string exp_name, string gene_id, int start, int end,
                      const char *sql_string) nogil:
    if not reads:
        return 0

    cdef:
        int rc
        int arg_len
        int rm_chars_len
        char *sql

    arg_len = int_len(reads) + exp_name.length() + ((gene_id.length() + int_len(start) + int_len(end)) * 2)
    rm_chars_len = 8 * 2
    sql = <char *> malloc(sizeof(char) * (strlen(sql_string) + arg_len - rm_chars_len + 1))
    sprintf(sql, sql_string, reads, exp_name.c_str(), gene_id.c_str(), start, end, gene_id.c_str(), start, end)
    rc = exec_db(db, sql)
    free(sql)
    return rc

cdef int junction_reads(sqlite3 *db, int reads, string exp_name, string junc_gene_id, int junc_start,
                        int junc_end) nogil:
    return reads_update(db, reads, exp_name, junc_gene_id, junc_start, junc_end, junc_reads_insert)

cdef int intron_retention_reads(sqlite3 *db, int reads, string exp_name, string ir_gene_id, int ir_start,
                                int ir_end) nogil:
    return reads_update(db, reads, exp_name, ir_gene_id, ir_start, ir_end, ir_reads_insert)
