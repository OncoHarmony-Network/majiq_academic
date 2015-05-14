# Parser
ANALYSIS_PSI = 'psi'
ANALYSIS_PSI_GENE = 'psi-gene'
ANALYSIS_DELTAPSI = 'deltapsi'
ANALYSIS_DELTAPSI_GENE = 'deltapsi-gene'
LSV_THUMBNAILS = 'lsv-thumbnails'
SPLICE_GRAPHS = 'splicegraphs'

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

IR_TYPE_START = 1
IR_TYPE_END = 2

# File extensions, delimiters, etc.
SUFFIX_SPLICEGRAPH = 'splicegraph'
DELIMITER = '\t'
EXTENSION = 'txt'

# URL composite (OS dependent, right now for MacOS)
URL_COMPOSITE = "file://%s#%s"

# subfolder for full (splicegraph summaries)
DELTAPSI_SUBFOLDERS = 'summaries'

# Debugging
DEBUG = 1