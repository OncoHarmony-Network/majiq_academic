
from libcpp.vector cimport vector
from majiq.src.internals.grimoire cimport  Gene
from libcpp.string cimport string
from libcpp.map cimport map
from majiq.src.internals.psi cimport psi_distr_t
cimport numpy as np
import numpy as np
from majiq.src.internals.qLSV cimport hetLSV, qLSV

cdef int read_gff(str filename, map[string, Gene*]& all_genes, vector[string]& gid_vec, object logging) except -1
cdef void get_coverage_mat(map[string, vector[psi_distr_t]]& result, map[string, int] lsv_map, list file_list,
                            str weight_fname, int nthreads)
cdef void get_coverage_mat_lsv(map[string, qLSV*]& result, list file_list, str weight_fname, int nthreads)
cdef int dump_lsv_coverage(str filename, dict cov_dict, list type_list, list junc_info, str exp_name)
cdef int dump_lsv_coverage_mat(str filename, list cov_list, list type_list, list junc_info, str exp_name)
cdef int dump_lsv_coverage_mat(str filename, list cov_list, list type_list, list junc_info, str exp_name)
cdef int dump_hettmp_file(str fname, np.ndarray[np.float32_t, ndim=2, mode="c"] osamps)
# cdef map[string, vector[psi_distr_t]] get_coverage_lsv(list list_of_lsv_id, list file_list, str weight_fname)