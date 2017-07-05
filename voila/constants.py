import multiprocessing
import os
from collections import namedtuple

import collections

import voila

VERSION = '1.0.1'
FILE_VERSION = '0.1'
ANALYSIS_PSI = 'psi'
ANALYSIS_PSI_GENE = 'psi-gene'
ANALYSIS_DELTAPSI = 'deltapsi'
ANALYSIS_DELTAPSI_GENE = 'deltapsi-gene'
ANALYSIS_HETEROGEN = 'heterogen'
LSV_THUMBNAILS = 'lsv-thumbnails'
SPLICE_GRAPHS = 'splice-graphs'
COND_TABLE = 'cond-table'
TOOLS = 'tools'

# Junction Types
JUNCTION_TYPE_DB_RNASEQ = 0
JUNCTION_TYPE_RNASEQ = 1
JUNCTION_TYPE_DB = 2
JUNCTION_TYPE_DB_OTHER_RNASEQ = 3

# Exon Types
EXON_TYPE_DB_RNASEQ = 0
EXON_TYPE_RNASEQ = 1
EXON_TYPE_DB = 2
EXON_TYPE_DB_OTHER_RNASEQ = 3
EXON_TYPE_MISSING_START = 4
EXON_TYPE_MISSING_END = 5

# Intron Retention Types
NONE_IR_TYPE = 0
IR_TYPE_START = 1
IR_TYPE_END = 2

# Summary constants.
SUFFIX_SPLICEGRAPH = 'splicegraph'
DELIMITER = '\t'
EXTENSION = 'txt'
MAX_GENES = 10 # Max. 10 genes per page, create as many HTMLs as needed
MAX_LSVS_DELTAPSI_INDEX = 10000  # Max. LSVs allowed to create full index.html
MAX_LSVS_PSI_INDEX = 15000
COMBINED_PREFIX = 'ALL_'

# subfolder for full (splicegraph summaries)
SUMMARIES_SUBFOLDER = 'summaries'

# Debugging
DEBUG = 1

# logging
VOILA_LOG_NAME = '3a4f4528-e572-404e-b143-acff61cee9ed'

# HDF5
EXPERIMENT_NAMES = 'experiments_names'

# MAJIQ-SPEL
LSV_TEXT_VERSION = 2

SPLICE_GRAPH_FILE_VERSION = 2
VOILA_FILE_VERSION = 3

FILE_MODE = collections.namedtuple('FILE_MODE', 'read write')('r', 'w')

EXEC_DIR = os.path.dirname(os.path.abspath(voila.__file__))
