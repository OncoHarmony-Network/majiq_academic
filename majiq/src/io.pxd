
from libcpp.vector cimport vector
from majiq.src.internals.grimoire cimport  Gene

cdef int read_gff(str filename, vector[Gene *]& glist, object logging) except -1
cdef int dump_lsv_coverage(str filename, dict cov_dict, list type_list, list junc_info, str exp_name)