#!/usr/bin/python
import gc
import numpy as np

from majiq.grimoire.exon import ExonTx, collapse_list_exons
from majiq.grimoire.lsv import LSV
from majiq.grimoire.junction import Junction
from majiq.src import config as majiq_config


class Gene:
    __eq__ = lambda self, other: (self.chromosome == other.chromosome and self.strand == other.strand
                                  and self.start < other.end and self.end > other.start
                                  and self.strand == other.strand)
    __ne__ = lambda self, other: (self.chromosome != other.chromosome or self.strand != other.strand
                                  or self.start >= other.end or self.end <= other.start)
    __lt__ = lambda self, other: (self.chromosome < other.chromosome
                                  or (self.chromosome == other.chromosome
                                      and (self.end < other.start
                                           or (self.end > other.start and self.start < other.end
                                               and self.strand == '+' and other.strand == '-'))))
    __gt__ = lambda self, other: (self.chromosome > other.chromosome
                                  or (self.chromosome == other.chromosome
                                      and (self.start > other.end
                                           or (self.end > other.start and self.start < other.END
                                               and self.strand == '-' and other.strand == '+'))))

    def __init__(self, gene_id, gene_name, chrom, strand, start, end):
        self.id = gene_id
        self.name = gene_name
        self.chromosome = chrom
        self.strand = strand
        self.transcript_tlb = {}
        self.exons = []
        self.start = start
        self.end = end
        self.readNum = np.zeros(shape=majiq_config.num_experiments, dtype=np.int)
        self.temp_txex_list = []
        self.ir_list = []
        self.ir_definition = []
        self.lsv_list = []
        self.antis_gene = []

    def __hash__(self):
        return hash((self.id, self.chromosome, self.strand, self.start, self.end))

    def get_id(self):
        return self.id

    def get_name(self):
        return self.name

    def get_strand(self):
        return self.strand

    def get_read_count(self):
        return self.readNum

    # def get_rpkm(self):
    # return self.RPKM

    def get_exon_list(self):
        return self.exons

    def get_chromosome(self):
        return self.chromosome

    def get_coordinates(self):
        return self.start, self.end

    def get_transcript(self, trans_id):
        return self.transcript_tlb[trans_id]

    def get_exon_in_coord(self, coord):
        res = None
        for ex in self.exons:
            cc = ex.get_coordinates()
            if cc[0] <= coord <= cc[1]:
                res = ex
                break
        return res

    def get_exon_by_id(self, ex_id):
        # return self.exons[ex_id-1]
        for ex in self.exons:
            if ex.get_id() == ex_id:
                res = ex
                break
        else:
            res = None

        return res

    def get_overlapped_genes(self):
        return self.antis_gene

    def get_ir_definition(self):
        return self.ir_definition

    ''' Set functions '''

    def add_transcript(self, tcrpt):
        if tcrpt.txstart < self.start:
            self.start = tcrpt.txstart
        if tcrpt.txend > self.end:
            self.end = tcrpt.txend

        self.transcript_tlb[tcrpt.get_id()] = tcrpt
        return

    # def add_intron_retention(self, lsv_ir):
    #     self.ir_list.append(lsv_ir)

    def add_read_count(self, read_num, exp_idx):
        self.readNum[exp_idx] += read_num
        return

    def add_exon(self, exon):
        self.exons.append(exon)
        return

    def add_ir_definition(self, start, end):
        self.ir_definition.append((start, end))

    def exist_antisense_gene(self, list_of_genes):
        # strnd = '-'
        # if self.strand == '-':
        #     strnd = '+'
        for strnd in ('+', '-'):
            for gg in list_of_genes[strnd]:
                if gg.get_id() == self.id:
                    continue
                if self.overlaps(gg):
                    # coords = gg.get_coordinates()
                    # if self.start < coords[1] and self.end > coords[0]:
                    self.set_antisense_gene(gg.get_id())
                    gg.set_antisense_gene(self.id)

    def set_antisense_gene(self, gn_id):
        self.antis_gene.append(gn_id)

    def overlaps(self, gne):
        if self.start < gne.end and self.end > gne.start:
            res = True
        else:
            res = False
        return res

    def check_antisense_junctions(self, jstart, jend):
        res = False
        for anti_g in self.antis_gene:
            # if not self.antis_gene is None:
            gg = majiq_config.gene_tlb[anti_g]
            cc = gg.get_coordinates()
            if jend < cc[0] or jstart > cc[1]:
                continue
            j_list = gg.get_all_junctions()
            for jj in j_list:
                if not jj.is_annotated():
                    continue
                (j_st, j_ed) = jj.get_coordinates()
                if j_st > jstart or (j_st == jstart and j_ed > jend):
                    break
                elif jstart == j_st and jend == j_ed:
                    res = True
                    break
            if not res:
                for ex in gg.get_exon_list():
                    coords = ex.get_coordinates()
                    # if jend > coords[0]:
                    #     break
                    coords = [coords[0] - majiq_config.get_max_denovo_difference(),
                              coords[1] + majiq_config.get_max_denovo_difference()]

                    if coords[0] <= jend <= coords[1] or coords[0] <= jstart <= coords[1]:
                        res = True
                        break
            else:
                break
        return res

    # def is_gene_in_list(self, list_of_genes, name):
    #     res = None
    #     for ll in list_of_genes:
    #         if self.chromosome == ll.chromosome and self.strand == ll.strand \
    #                 and self.start < ll.end and self.end > ll.start:
    #
    #             res = ll
    #             ll.start = min(ll.start, self.start)
    #             ll.end = max(ll.end, self.end)
    #             break
    #
    #     return res

    def exist_exon(self, start, end):
        """
         .. function: exist_exon (self, start, end):
            Check if the pair (start, end) are in a known exon in the gene. If not return None. We assume
            that the given exon is for the same chromosome and strand than the gene.

            :param start: Start position of the input exon.
            :param end: End Position of te input exon.
            :rtype: exon instance or None
        """

        res = None

        for ee in self.exons:
            ex_str, ex_end = ee.get_coordinates()
            if start < ex_end and end > ex_str:
                res = ee
                break

        return res

    def exist_junction(self, start, end):
        if start is None or end is None:
            return
        res = None
        for txcpt in self.transcript_tlb.values():
            res = txcpt.in_junction_list(start, end)
            if not res is None:
                break
        return res

    def get_all_junctions(self):
        lst = set()
        for ex in self.get_exon_list():
            for ex_rna in ex.exonRead_list:
                if len(ex_rna.p5_junc) > 0:
                    lst = lst.union(set(ex_rna.p5_junc))

        for tt in self.transcript_tlb.values():
            for jj in tt.get_junction_list():
                if not jj is None and not jj in lst:
                    lst.add(jj)
        s_junc = list(lst)
        return sorted([xx for xx in s_junc if not xx is None])

    def prepare_exons(self):
        self.exons.sort()

    def get_all_ss(self):
        ss = set()
        for ex in self.exons:
            tx_list = ex.get_annotated_exon()
            for txex in tx_list:
                for junc in txex.get_3prime_junc():
                    coord = junc.get_ss_3p()
                    if not coord is None:
                        ss.add((coord, '3prime', junc))
                for junc in txex.get_5prime_junc():
                    coord = junc.get_ss_5p()
                    if not coord is None:
                        ss.add((coord, '5prime', junc))

        return sorted(ss)

    def collapse_exons(self):

        self.temp_txex_list.sort()
        list_ex = collapse_list_exons(self.temp_txex_list, self)
        self.exons.extend(list_ex)
        self.prepare_exons()
        self.remove_temp_attributes()

    def check_exons(self):

        s_exons = set()

        for ex in self.exons:
            s_exons.add(ex.get_coordinates())

        if len(s_exons) != len(self.exons):
            for xx in self.exons:
                print xx, xx.get_coordinates()

            raise RuntimeError

    def new_lsv_definition(self, exon, jlist, lsv_type):

        coords = exon.get_coordinates()
        ret = None
        lsv_id = "%s:%d-%d:%s" % (self.get_id(), coords[0], coords[1], lsv_type)
        for lsv in self.lsv_list:
            if lsv.id == lsv_id:
                ret = lsv
                break
        else:
            ret = LSV(exon, lsv_id, jlist, lsv_type)
            self.lsv_list.append(ret)
        return ret

    def remove_temp_attributes(self):
        del self.temp_txex_list

    def new_annotated_exon(self, start, end, transcript, bl=True, intron=False):
        for txex in self.temp_txex_list:
            if txex.start == start and txex.end == end:
                res = txex
                break
        else:
            res = ExonTx(start, end, transcript, intron=intron)
            if bl:
                self.temp_txex_list.append(res)
        return res

    def new_annotated_junctions(self, start, end, trcpt):
        junc = self.exist_junction(start, end)
        if junc is None:
            junc = Junction(start, end, None, None, self, annotated=True)
        junc.add_transcript(trcpt)

        return junc

    def get_all_introns(self):
        lintrons = []
        if len(self.exons) > 1:
            for ex_idx, ex in enumerate(self.exons[:-1]):
                lintrons.append((ex, self.exons[ex_idx+1]))
        return lintrons

#     def get_rnaseq_mat(self, rand10k, use_annot=True):
#
#         ss3_l = []
#         ss5_l = []
#         tlb = {}
#         exidx = 0
#         ss_3p_vars = [0]*20
#         ss_5p_vars = [0]*20
#         ss_both_var = 0
#         exon_list = []
#         ex_list = self.get_exon_list()
#         for ex in ex_list:
# #            if ex.id is None:
# #                continue
#             l3 = len(set(ex.ss_3p_list))
#             l5 = len(set(ex.ss_5p_list))
#             if l3 == 0 or l5 == 0:
#                 continue
#
#             local_3p, local_5p = ex.ss_variant_counts()
#             if local_3p > 1 and local_5p > 1:
#                 ss_both_var += 1
#
#             ss_3p_vars[local_3p] += 1
#             ss_5p_vars[local_5p] += 1
#
#             st3 = len(ss3_l)
#             st5 = len(ss5_l)
#             ss3_l += sorted([ss3 for ss3 in set(ex.ss_3p_list)])
#             ss5_l += sorted([ss5 for ss5 in set(ex.ss_5p_list)])
#             tlb[exidx] = [range(st3, len(ss3_l)), range(st5, len(ss5_l))]
#             exon_list.append(ex)
#             exidx += 1
#
#         mat = np.zeros(shape=(len(ss5_l), len(ss3_l)), dtype='int')
#         # jmat = np.empty(shape=(len(ss5_l), len(ss3_l)), dtype='object')
#         # jmat.fill(None)
#
#         junc_list = self.get_all_junctions()
#         for junc in junc_list:
#             st, end = junc.get_coordinates()
#             if not st in ss5_l or not end in ss3_l:
#                 continue
#             x = ss5_l.index(st)
#             y = ss3_l.index(end)
#
#             read_num = junc.get_coverage().sum()
#             if read_num > 0:
#                 count_mat = read_num
#             elif junc.is_annotated() and use_annot:
#                 count_mat = -1
#             else:
#                 count_mat = 0
#             mat[x, y] = count_mat
#             # jmat[x, y] = junc
#
#             for exp_idx in range(majiq_config.num_experiments):
#                 if reliable_in_data(junc, exp_idx):
#                     rand10k[exp_idx].add(junc)
#
#         return mat, exon_list, tlb, [ss_3p_vars, ss_5p_vars, ss_both_var]


class Transcript(object):

    def __init__(self, name, gene, txstart, txend):
        self.id = name
        self.gene_id = gene.get_id()
        self.exon_list = []
        #self.junction_list = []
        self.txstart = txstart
        self.txend = txend
        # self.cdsStart = None
        # self.cdsStop = None

    def get_gene(self):

        return majiq_config.gene_tlb[self.gene_id]

    def get_exon_list(self):
        return self.exon_list

    def get_id(self):
        return self.id

    def get_junction_list(self):
        res = set()
        for ex in self.exon_list:
            res.update(set(ex.get_5prime_junc()))
            res.update(set(ex.get_3prime_junc()))

        return sorted(res)

    def prepare_exon_list(self):
        self.exon_list.sort()
        return self.exon_list

    def in_junction_list(self, start, end):
        res = None
        for ex in self.exon_list:
            for jj in ex.get_5prime_junc():
                jjstart, jjend = jj.get_coordinates()
                if jjstart == start and jjend == end:
                    res = jj
                    break
            else:
                continue
            break
        return res

    def add_exon(self, exon):
        self.exon_list.append(exon)
        return

    # def add_junction(self, junc):
    #     if junc not in self.junction_list:
    #         self.junction_list.append(junc)

#     def sort_in_list(self):
#         # strand = self.get_gene().get_strand()
#         # if strand == '+':
#         #     isneg = False
#         # else:
#         #     isneg = True
# #        self.exon_list.sort(reverse=isneg)
#         self.exon_list.sort()


def recreate_gene_tlb(gene_list):

    for gn in gene_list:
        majiq_config.gene_tlb[gn.get_id()] = gn


def clear_gene_tlb():
    majiq_config.gene_tlb.clear()
    gc.collect()