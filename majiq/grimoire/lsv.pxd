from majiq.grimoire.exon cimport Exon
import collections
cimport numpy as np

cdef class LSV:
    cdef list junctions
    cdef Exon exon
    cdef str type
    cdef str gene_id
    cdef str chromosome
    cdef str strand
    cdef str id

    #def get_visual_lsv(self)
    cdef tuple sample_lsvs(LSV self, np.ndarray junc_mtrx, float fitfunc_r, object majiq_config)
    cdef str set_type(LSV self, list jlist, Exon ref_exon, str gene_strand, bint ss)
    cdef int add_lsv(LSV self, np.ndarray junc_mtrx, list type_dict, dict values, list junc_info, float fitfunc_r,
                     object majiq_config) except -1
    #cdef int to_hdf5(LSV self, hdf5grp, int lsv_idx)

quant_lsv = collections.namedtuple('quant_lsv', 'id type coverage')