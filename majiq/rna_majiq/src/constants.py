"""
constants.py

Constants used by MAJIQ
"""

import os
import tempfile
import string as pystring  # just importing it as string causes cython errors
from importlib.metadata import version, PackageNotFoundError

try:
    VERSION = version("rna_majiq_meta")
except PackageNotFoundError:
    try:
        VERSION = version("rna_majiq")
    except PackageNotFoundError:
        VERSION = "2.3.0"

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

# majiq het constants
HET_SAMPLING_SEED = 20200401


try:
    import numpy as np

    EPSILON = np.finfo(np.float64).eps
except ModuleNotFoundError:  # when importing in setup, don't need numpy yet
    EPSILON = 2e-16  # slightly less than float64 epsilon

run_tempdir = None
def get_tmp_dir(outdir):
    global run_tempdir
    if not run_tempdir:
        run_tempdir = tempfile.mkdtemp(dir=outdir)
    return run_tempdir

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
    return os.path.join(get_tmp_dir(outdir), "db.tb")


def get_builder_majiq_filename(outdir, name):
    return "%s/%s.majiq" % (outdir, name)


def get_builder_splicegraph_filename(outdir):
    return "%s/splicegraph.sql" % (outdir)


def get_weights_filename(outdir, name):
    return "%s/%s.wght" % (outdir, name)


def get_tmp_psisample_file(outdir, name):
    return os.path.join(get_tmp_dir(outdir), f"{name}.psisamples.tmp")
