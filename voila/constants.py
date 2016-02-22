# Parser
ANALYSIS_PSI = 'psi'
ANALYSIS_PSI_GENE = 'psi-gene'
ANALYSIS_DELTAPSI = 'deltapsi'
ANALYSIS_DELTAPSI_GENE = 'deltapsi-gene'
LSV_THUMBNAILS = 'lsv-thumbnails'
SPLICE_GRAPHS = 'splicegraphs'
COND_TABLE = 'cond-table'

# LSV Graphics
JUNCTION_TYPE_DB_RNASEQ = 0
JUNCTION_TYPE_RNASEQ = 1
JUNCTION_TYPE_DB = 2
JUNCTION_TYPE_DB_OTHER_RNASEQ = 3

EXON_TYPE_DB_RNASEQ = 0
EXON_TYPE_RNASEQ = 1
EXON_TYPE_DB = 2
EXON_TYPE_DB_OTHER_RNASEQ = 3
EXON_TYPE_MISSING_START = 4
EXON_TYPE_MISSING_END = 5

NONE_IR_TYPE = 0
IR_TYPE_START = 1
IR_TYPE_END = 2

# Summary constants.
SUFFIX_SPLICEGRAPH = 'splicegraph'
DELIMITER = '\t'
EXTENSION = 'txt'
MAX_GENES = 10  # Max. 10 genes per page, create as many HTMLs as needed
MAX_LSVS_DELTAPSI_INDEX = 10000 # Max. LSVs allowed to create full index.html
MAX_LSVS_PSI_INDEX = 15000
COMBINED_PREFIX = 'ALL_'

# URL composite (OS dependent, right now for MacOS)
URL_COMPOSITE = "file://%s#%s"

# subfolder for full (splicegraph summaries)
SUMMARIES_SUBFOLDER = 'summaries'

# Debugging
DEBUG = 1