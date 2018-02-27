from majiq.src.constants import *
from majiq.src.sample import sample_from_junctions
from majiq.grimoire.exon cimport Exon
from majiq.grimoire.junction cimport Junction
from majiq.src.normalize import mark_stacks
cimport numpy as np
import numpy as np
from voila.constants import *
import collections

quant_lsv = collections.namedtuple('quant_lsv', 'id coverage')

ctypedef np.float64_t DTYPE_t

class InvalidLSV(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


cdef class LSV:
    def __init__(self, gene_id, gene_chromosome, gene_strand, ex, ss):

        self.junctions = []
        if ss:
            jncs = [xx for xx in ex.ob if not (xx.donor is None or xx.acceptor is None)]
        else:
            jncs = [xx for xx in ex.ib if not (xx.donor is None or xx.acceptor is None)]

        if len(jncs) < 2:
            raise InvalidLSV("not enougth junctions")

        self.exon = ex
        self.type = self.set_type(jncs, ex, gene_strand, ss)
        if len(self.junctions) < 2:
            raise InvalidLSV("not enougth junctions")

        if len(self.type) >= 250:
            self.type = NA_LSV

        self.gene_id = gene_id
        self.chromosome = gene_chromosome
        self.strand = gene_strand
        stt = ex.start if ex.start > 0 else 'nan'
        end =  ex.end if ex.end > 0 else 'nan'
        self.id = "%s:%s:%s-%s" % (gene_id, self.type[0], stt, end)
        if self.type.endswith('i|i'):
            raise InvalidLSV("incorrect LSV, too many ir junctions")


    cdef int add_lsv(LSV self, np.ndarray junc_mtrx, list type_dict, dict values, list junc_info, float fitfunc_r,
                     object majiq_config) except -1:

        cdef list val
        cdef Junction xx

        type_dict.append((self.id, self.type))
        [junc_info.append((self.id, xx.start, xx.end, junc_mtrx[xx.index].sum(),
                           np.count_nonzero(junc_mtrx[xx.index]))) for xx  in self.junctions]
        val = [junc_mtrx[xx.index] for xx  in self.junctions]
        values[self.id] = sample_junctions(np.array(val), fitfunc_r, majiq_config)

    cdef tuple sample_lsvs(LSV self, np.ndarray junc_mtrx, float fitfunc_r, object majiq_config):
        cdef Junction xx
        cdef list ex_index
        cdef np.ndarray s_lsv, lsv_trs

        ex_index = [xx.index for xx in self.junctions]
        cover = junc_mtrx[ex_index]

        mark_stacks(cover, fitfunc_r, majiq_config.pvalue_limit)

        s_lsv = sample_from_junctions(junction_list=cover,
                                      m=majiq_config.m,
                                      k=majiq_config.k,
                                      fitted_one_over_r=fitfunc_r)

        lsv_trs = np.array([cover.sum(axis=1), np.count_nonzero(cover, axis=1)]).T
        return s_lsv, lsv_trs

    cdef str set_type(LSV self, list jlist, Exon ref_exon, str gene_strand, bint ss):

        cdef list sp_list = []
        cdef str ref_exon_id = "%s-%s" %(ref_exon.start, ref_exon.end)
        cdef Junction junc, xx_junc
        cdef Exon ex
        cdef int coord, ref_coord
        cdef tuple xx

        cdef str ext_type, prev_ex, exid
        cdef int excount, jidx
        cdef list ss_list
        cdef list ref_ss
        cdef set ref_ss_set

        for junc in jlist:
            if ss:
                ex = junc.acceptor
                coord = junc.end
                ref_coord = junc.start
            else:
                ex = junc.donor
                coord = junc.start
                ref_coord = junc.end
            sp_list.append((ex.intron, coord, ref_coord, ex, junc))

        if gene_strand == '-':
            sp_list.sort(key=lambda xx: (-xx[0], xx[1], xx[2]), reverse=True)
        else:
            sp_list.sort(key=lambda xx: (xx[0], xx[1], xx[2]), reverse=False)

        ext_type = 's' if (ss and gene_strand != '-') or (not ss and gene_strand == '-') else 't'
        prev_ex = "%s-%s" % (sp_list[0][3].start, sp_list[0][3].end)
        excount = 1

        ref_ss_set = set([junc.start for junc in ref_exon.ob]) if ss else set([junc.end for junc in ref_exon.ib])
        ref_ss = sorted(list(ref_ss_set))

        for intron, coord, ref_coord, ex, junc in sp_list:
            exid = "%s-%s" % (ex.start, ex.end)
            jidx = ref_ss.index(ref_coord)
            jidx += 1

            if exid == ref_exon_id:
                continue

            if intron:
                ext_type += "|i"
                self.junctions.append(junc)
                continue

            if ex is None:
                ext_type += "|%se0" % jidx
                continue

            if ss:
                ss_list = sorted([xx_junc.end for xx_junc in ex.ib])
            else:
                ss_list = sorted([xx_junc.start for xx_junc in ex.ob])

            if prev_ex != exid:
                prev_ex = exid
                excount += 1
            ext_type += "|%se%s.%so%s" % (jidx, excount, ss_list.index(coord)+1, len(ss_list))

            self.junctions.append(junc)
        return ext_type


cdef int _detect_lsvs(list list_exons, np.ndarray junc_mtrx, float fitfunc_r, str gid, str gchrom, str gstrand,
                 object majiq_config, dict lsv_dict, list type_tlb, list junc_info, object logger) except -1:

    cdef int count = 0
    cdef np.ndarray sum_trx = junc_mtrx.sum(axis=1)
    cdef np.ndarray pos_trx = np.count_nonzero(junc_mtrx, axis=1)
    cdef list lsv_list = [[], []]
    cdef Exon ex
    cdef int ii, lsv_idx
    cdef set jjset
    cdef list jjlist, ex_index
    cdef np.ndarray ex_mtrx_s, ex_mtrx_p, b, c
    cdef LSV ss, st

    for ex in list_exons:
        for ii, jjset in enumerate([ex.ib, ex.ob]):
            jjlist = [xx for xx in jjset if not (xx.donor is None or xx.acceptor is None)]

            if len(jjlist) < 2:
                continue
            ex_index = [xx.index for xx in jjlist]
            ex_mtrx_s = sum_trx[ex_index]
            ex_mtrx_p = pos_trx[ex_index]

            if np.any(np.logical_and(ex_mtrx_s > majiq_config.minpos, ex_mtrx_p > majiq_config.minreads)):
                try:
                    lsv_list[ii].append(LSV(gid, gchrom, gstrand, ex, ss=(ii == 1)))
                except InvalidLSV:
                    continue

    for ss in lsv_list[0]:
        for st in lsv_list[1]:
            if set(ss.junctions).issubset(set(st.junctions)) and not set(ss.junctions).issuperset(set(st.junctions)):
                break
        else:
            ss.add_lsv(junc_mtrx, type_tlb, lsv_dict, junc_info, fitfunc_r, majiq_config)
            count += 1

    for st in lsv_list[1]:
        for ss in lsv_list[0]:
            if set(st.junctions).issubset(set(ss.junctions)):
                break
        else:
            st.add_lsv(junc_mtrx, type_tlb, lsv_dict, junc_info, fitfunc_r, majiq_config)
            count += 1

    return count

###API

cpdef detect_lsvs(object dict_of_genes, dict dict_junctions, dict list_exons, np.ndarray junc_mtrx, float fitfunc_r,
                  object majiq_config, dict lsv_dict, list lsv_type_list, list junc_info, object logger):

    cdef str gne_id
    cdef int gne_idx
    cdef dict gene_obj

    for gne_idx, (gne_id, gene_obj) in enumerate(dict_of_genes.items()):
        # if gene_obj['nreads'] > 0:
        _detect_lsvs(list_exons[gne_id], junc_mtrx, fitfunc_r, gne_id, gene_obj['chromosome'],
                     gene_obj['strand'], majiq_config, lsv_dict, lsv_type_list, junc_info, logger)
        for jj in dict_junctions[gne_id].values():
            jj.reset()
        gene_obj['nreads'] = 0

cpdef tuple sample_junctions_old(np.ndarray junc_mtrx, float fitfunc_r, object majiq_config):
    cdef Junction xx
    cdef list ex_index
    cdef np.ndarray s_lsv, lsv_trs

    lsv_trs = np.array([junc_mtrx.sum(axis=1), np.count_nonzero(junc_mtrx, axis=1)]).T
    mark_stacks(junc_mtrx, fitfunc_r, majiq_config.pvalue_limit)
    s_lsv = sample_from_junctions(junction_list=junc_mtrx,
                                  m=majiq_config.m,
                                  k=majiq_config.k,
                                  fitted_one_over_r=fitfunc_r)

    return s_lsv

cpdef np.ndarray[DTYPE_t, ndim=2] sample_junctions(np.ndarray junc_mtrx, float fitfunc_r, object majiq_config):
    cdef Junction xx
    cdef list ex_index
    cdef np.ndarray s_lsv, lsv_trs

    mark_stacks(junc_mtrx, fitfunc_r, majiq_config.pvalue_limit)
    s_lsv = sample_from_junctions(junction_list=junc_mtrx,
                                  m=majiq_config.m,
                                  k=majiq_config.k,
                                  fitted_one_over_r=fitfunc_r)



    return s_lsv
