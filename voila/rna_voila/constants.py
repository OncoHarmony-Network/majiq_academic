import os
import tempfile
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

import rna_voila

EXEC_DIR = Path(os.path.dirname(os.path.abspath(rna_voila.__file__)))

try:
    VERSION = version("rna_majiq_meta")
except PackageNotFoundError:
    try:
        VERSION = version("rna_voila")
    except PackageNotFoundError:
        VERSION = "2.3.0"

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
MAX_GENES = 1  # Max. 10 genes per page, create as many HTMLs as needed
MAX_LSVS_DELTAPSI_INDEX = 10000  # Max. LSVs allowed to create full het_index.html
MAX_LSVS_HET_INDEX = 10000  # Max. LSVs allowed to create full het_index.html
MAX_LSVS_PSI_INDEX = 15000
COMBINED_PREFIX = 'ALL_'

# subfolder for full (splice graph summaries)
SUMMARIES_SUBFOLDER = 'summaries'

# Debugging
DEBUG = 1

# logging
VOILA_LOG_NAME = '3a4f4528-e572-404e-b143-acff61cee9ed'

# MAJIQ-SPEL
LSV_TEXT_VERSION = 5

SPLICE_GRAPH_FILE_VERSION = 6
VOILA_FILE_VERSION = 6

NA_LSV = 'na'

# CONFIG_FILE = '/tmp/voila.ini'
# CONFIG_FILE = tempfile.NamedTemporaryFile(delete=False).name
if 'VOILA_CONFIG_FILE' not in os.environ or not os.environ['VOILA_CONFIG_FILE']:
    CONFIG_FILE = tempfile.NamedTemporaryFile(delete=False).name
else:
    CONFIG_FILE = os.environ['VOILA_CONFIG_FILE']

MINVAL = 1e-300
