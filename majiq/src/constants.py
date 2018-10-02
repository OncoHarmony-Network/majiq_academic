import os
VERSION = '2.0'

JUNC_FILE_FORMAT = 'sjdb'
SEQ_FILE_FORMAT = 'bam'
SEQ_INDEX_FILE_FORMAT = 'bam.bai'

UNSTRANDED = 0
FWD_STRANDED = 1
REV_STRANDED = 2

EX_TYPE = 0
IR_TYPE = 1
J_TYPE = 2

#QUEUE MESSAGE CODES

QUEUE_MESSAGE_END_WORKER = -1
QUEUE_MESSAGE_BUILD_JUNCTION = 0
QUEUE_MESSAGE_BUILD_INTRON = 1
QUEUE_MESSAGE_SPLICEGRAPH = 2
QUEUE_MESSAGE_PSI_RESULT = 3
QUEUE_MESSAGE_DELTAPSI_RESULT = 4
QUEUE_MESSAGE_BOOTSTRAP = 5
QUEUE_MESSAGE_HETER_DELTAPSI = 6


MIN_INTRON_LEN = 1000
MIN_BP_OVERLAP = 8
MIN_JUNC_LENGTH = 2
NUM_INTRON_BINS = 10
MAX_DENOVO_DIFFERENCE = 500
NRANDOM_JUNCTIONS = 5000

EMPTY_COORD = -1
FIRST_LAST_JUNC = -2

SIMPLIFY_ALL = 'all'
SIMPLIFY_DB = 'annotated'
SIMPLIFY_DENOVO = 'denovo'

WEIGTHS_AUTO = 'auto'
WEIGTHS_NONE = 'none'
MINVAL = 1e-300

GIT_VERSION_FILE = 'git_version'


all_stats = ["TNOM", "INFOSCORE", "WILCOXON"]


import numpy as np
EPSILON = np.finfo(np.float64).eps


def get_quantifier_voila_filename(outdir, name, deltapsi=False, het=False):
    if deltapsi:
        return "%s/%s_%s.deltapsi.voila" % (outdir, name[0], name[1])
    elif het:
        return "%s/%s_%s.het.voila" % (outdir, name[0], name[1])
    else:
        return "%s/%s.psi.voila" % (outdir, name)


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
    return "%s/%s.psisamples.tmp" %(outdir, name)


def store_git_version(short_sha):

    direc = os.path.dirname(__file__)
    with open('%s/../data/%s' % (direc, GIT_VERSION_FILE), 'w+') as ofp:
        ofp.write('%s\n' % short_sha)


def get_git_version():
    direc = os.path.dirname(__file__)
    try:
        with open('%s/../data/%s' % (direc, GIT_VERSION_FILE), 'r') as ofp:
            ver = ofp.readline().strip()
    except FileNotFoundError:
        ver = '<hash_not_found>'

    return ver
