VERSION = '1.0.0'

JUNCTIONS_DATASET_NAME = '/junctions'
JUNCTIONS_GC_CONTENT = '/junctions_gc'

#QUEUE MESSAGE CODES

QUEUE_MESSAGE_END_WORKER = -1
QUEUE_MESSAGE_BUILD_LSV = 0
QUEUE_MESSAGE_BUILD_CONST_JUNCTION = 1
QUEUE_MESSAGE_SPLICEGRAPH = 2
QUEUE_MESSAGE_PSI_RESULT = 3
QUEUE_MESSAGE_DELTAPSI_RESULT = 4
QUEUE_MESSAGE_BOOTSTRAP = 5


MIN_INTRON_LEN = 1000
MIN_BP_OVERLAP = 8
MIN_JUNC_LENGTH = 10
NUM_INTRON_BINS = 10
MAX_DENOVO_DIFFERENCE = 500
NRANDOM_JUNCTIONS = 5000


SIMPLIFY_ALL = 'all'
SIMPLIFY_DB = 'annotated'
SIMPLIFY_DENOVO = 'denovo'

WEIGTHS_AUTO = 'auto'
WEIGTHS_NONE = 'none'
MINVAL = 1e-300
import numpy as np
EPSILON = np.finfo(np.float64).eps


def get_quantifier_temp_filename(outdir, name):
    return "%s/tmp.%s.filtered.hdf5" % (outdir, name)


def get_quantifier_voila_filename(outdir, name, deltapsi=False):
    if deltapsi:
        return "%s/%s_%s.deltapsi.voila" % (outdir, name[0], name[1])
    else:
        return "%s/%s.psi.voila" % (outdir, name)


def get_prior_matrix_filename(outdir, names):
    return "%s/%s_%s.priormatrix.pkl" % (outdir, names[0], names[1])


def get_build_temp_db_filename(outdir):
    return "%s/tmp/db.hdf5" % outdir


def get_builder_majiq_filename(outdir, name):
    return "%s/%s.majiq.hdf5" % (outdir, name)


def get_builder_temp_majiq_filename(outdir, name):
    return "%s/tmp/%s.tmp.hdf5" % (outdir, name)


# def get_builder_splicegraph_filename(outdir, name):
#     return "%s/%s.splicegraph.hdf5" % (outdir, name)

def get_builder_splicegraph_filename(outdir):
    return "%s/splicegraph.hdf5" % (outdir)


def get_quantifier_norm_temp_files(outdir, name, replica_num):
    #return "%s/%s_%s.temp.hdf5" % (outdir, name, replica_num)
    return "%s/%s_%s.temp.hdf5" % (outdir, name, replica_num)
