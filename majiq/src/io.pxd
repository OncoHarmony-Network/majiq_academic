
from libcpp.vector cimport vector
from majiq.src.internals.grimoire cimport  Gene
from libcpp.string cimport string
from libcpp.map cimport map


cdef int read_gff(str filename, map[string, Gene*]& all_genes, vector[string]& gid_vec, object logging) except -1
cdef int dump_lsv_coverage(str filename, dict cov_dict, list type_list, list junc_info, str exp_name)