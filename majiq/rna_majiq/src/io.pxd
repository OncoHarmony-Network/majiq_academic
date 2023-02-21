from libcpp.vector cimport vector
from rna_majiq.src.internals.grimoire cimport  Gene
from libcpp.string cimport string
from libcpp.map cimport map
cimport numpy as np
import numpy as np
from rna_majiq.src.internals.qLSV cimport qLSV

cdef int read_gff(str filename, map[string, Gene*]& all_genes, vector[string]& gid_vec, bint simpl, bint enable_anot_ir,
                  object logging, bint ext3prime, bint ext5prime) except -1
cdef void get_coverage_mat_lsv(map[string, qLSV*]& result, list file_list, int nthreads, bint fltr,
                               int minreads, int minnonzero)
cdef int dump_lsv_coverage_mat(str filename, list cov_list, list type_list, list junc_info, str exp_name)
cdef int dump_hettmp_file(str fname, np.ndarray[np.float32_t, ndim=2, mode="c"] osamps)
