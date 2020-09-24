"""
constants.py

Constants used by MAJIQ
"""

import os
import string as pystring  # just importing it as string causes cython errors

VERSION = "2.2"

# file extensions
JUNC_FILE_FORMAT = "sj"
SEQ_FILE_FORMAT = "bam"
SEQ_INDEX_FILE_FORMAT = "bam.bai"
MAJIQ_FILE_FORMAT = "majiq"

# filename constants
GROUP_NAME_SEP = '-'
ALLOWED_GROUP_NAME_CHARS = set(pystring.ascii_letters + pystring.digits + '_')

# integer codes for strandedness
UNSTRANDED = 0
FWD_STRANDED = 1
REV_STRANDED = 2

# number of reads to use to estimate read length at first for io_bam
ESTIMATE_NUM_READS = 100

# special coordinates
EMPTY_COORD = -1
FIRST_LAST_JUNC = -2

# file for git hash
GIT_VERSION_FILE = "git_version"

# majiq het constants
HET_SAMPLING_SEED = 20200401
all_stats = ["TNOM", "INFOSCORE", "WILCOXON"]



try:
    import numpy as np

    EPSILON = np.finfo(np.float64).eps
except ModuleNotFoundError:  # when importing in setup, don't need numpy yet
    EPSILON = 2e-16  # slightly less than float64 epsilon


def get_quantifier_voila_filename(outdir, name, deltapsi=False, het=False):
    if deltapsi:
        return f"{outdir}/{name[0]}{GROUP_NAME_SEP}{name[1]}.deltapsi.voila"
    elif het:
        return f"{outdir}/{name[0]}{GROUP_NAME_SEP}{name[1]}.het.voila"
    else:
        return f"{outdir}/{name}.psi.voila"

def get_prior_matrix_filename(outdir, names):
    return "%s/%s_%s.priormatrix.pkl" % (outdir, names[0], names[1])


def get_build_temp_db_filename(outdir):
    return "%s/db.tb" % outdir


def get_builder_majiq_filename(outdir, name):
    return "%s/%s.majiq" % (outdir, name)


def get_builder_splicegraph_filename(outdir):
    return "%s/splicegraph.sql" % (outdir)


def get_weights_filename(outdir, name):
    return "%s/%s.wght" % (outdir, name)


def get_tmp_psisample_file(outdir, name):
    return "%s/%s.psisamples.tmp" % (outdir, name)


def store_git_version(short_sha):

    direc = os.path.dirname(__file__)
    with open("%s/../data/%s" % (direc, GIT_VERSION_FILE), "w+") as ofp:
        ofp.write("%s\n" % short_sha)


def get_git_version():
    direc = os.path.dirname(__file__)
    try:
        with open("%s/../data/%s" % (direc, GIT_VERSION_FILE), "r") as ofp:
            ver = ofp.readline().strip()
    except FileNotFoundError:
        ver = "<hash_not_found>"

    return ver
