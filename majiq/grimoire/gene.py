#!/usr/bin/python
import numpy as np

from majiq.grimoire.exon import ExonTx, collapse_list_exons
from majiq.grimoire.lsv import LSV, InvalidLSV
from majiq.grimoire.junction import Junction
from majiq.src import config
from majiq.src.analize import reliable_in_data



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
                                           or (self.end > other.start and self.start < other.end
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
        self.readNum = np.zeros(shape=config.num_experiments, dtype=np.int)
        self.temp_txex_list = []
        self.ir_list = []
        self.ir_definition = []
        self.lsv_list = []
        # self.RPKM = np.zeros(shape=mglobals.num_experiments, dtype=np.float)
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

    # def get_transcript_AS_candidates(self):
    #     return self.transAScandidates

    # def get_transcript_CONST_candidates(self):
    #     return self.transCONSTcandidates

    ''' Set functions '''

    # def add_transcript_AS_candidates(self,list_candidates):
    #     self.transAScandidates += list_candidates
    #     return
    #
    # def add_transcript_CONST_candidates(self,list_candidates):
    #     self.transCONSTcandidates += list_candidates
    #     return

    def add_transcript(self, tcrpt):
        if tcrpt.txstart < self.start:
            self.start = tcrpt.txstart
        if tcrpt.txend > self.end:
            self.end = tcrpt.txend

        self.transcript_tlb[tcrpt.get_id()] = tcrpt
        return

    def add_intron_retention(self, lsv_ir):
        self.ir_list.append(lsv_ir)

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
            gg = config.gene_tlb[anti_g]
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
                    if coords[0] <= jend <= coords[1] or coords[0] <= jstart <= coords[1]:
                        res = True
                        break
            else:
                break
        return res

    def is_gene_in_list(self, list_of_genes, name):
        res = None
        for ll in list_of_genes:
            if self.chromosome == ll.chromosome and self.strand == ll.strand \
                    and self.start < ll.end and self.end > ll.start:
                res = ll

                ll.start = min(ll.start, self.start)
                ll.end = max(ll.end, self.end)
                break
        return res

        # def calculate_rpkm(self, experiment_index, total_reads):

    # """
    #          .. function: calculate_RPKM( self, experiment_index, total_Reads )
    #
    #             This function calculates the RPKM as follows
    #             rpk = #Reads in Gene / # of kilobases of the gene exons (some may don't have any read)
    #             rpkm = rpk / (#total reads/1000000)
    #
    #             :param experiment_index: Index of the experiment from the origrinal experiment list
    #             :param total_reads: Total Number of reads of the gene.
    #             :rtype: RPKM value for this gene
    #         """
    #         if len(self.exons) == 0:
    #             return 0
    #         total_kb = float(0)
    #         for ex in self.exons:
    #             start, end = ex.get_coordinates()
    #             #            print "EXON ::",ex.id,end, start
    #             total_kb += float(end-start)
    #
    # #        print self.readNum, experiment_index
    #         rpk = float(self.readNum[experiment_index]) / float(total_kb/1000)
    #         mreads = float(total_reads)/float(1000000)
    #         rpkm = float(rpk)/mreads
    #         #print "Strand",self.strand,"::",total_kb, self.readNum, rpk, mreads, rpkm
    #         self.RPKM[experiment_index] = rpkm
    #         return rpkm

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
#            print "EX GEN:",ee.start, ee.end
#            print "New EX:",start, end
            if start < ee.end and end > ee.start:
                res = ee
                break
#                fnd +=1
#            elif (end < ee.start):
#                if fnd >1 : # overlapping more than one exon.... intron retention
#                    res = None
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

    def get_annotated_junctions(self):
        lst = set()
        for tt in self.transcript_tlb.values():
            for jj in tt.get_junction_list():
                if not jj is None and not jj in lst:
                    lst.add(jj)
        s_junc = list(lst)
        return sorted(s_junc)

    def get_all_rna_junctions(self):
        lst = []
        for ex in self.get_exon_list():
            for ex_rna in ex.exonRead_list:
                if len(ex_rna.p5_junc) > 0:
                    lst.union(set(ex_rna.p5_junc))
        lst.sort()
        return lst

    def prepare_exons(self):
#        self.exons.sort(reverse = isneg)
        self.exons.sort()

    def get_all_ss(self, anot_only=False):

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


        #assert len(s_exons) == len(self.exons), "Exist duplicates in exons in Gene %s, " % self.id
        if len(s_exons) != len(self.exons):
            for xx in self.exons:
                print xx, xx.get_coordinates()

            raise RuntimeError

    def new_lsv_definition(self, exon, jlist, lsv_type, logger=None):

        coords = exon.get_coordinates()
        ret = None
        lsv_id = "%s:%d-%d:%s" % (self.get_id(), coords[0], coords[1], lsv_type)
        for lsv in self.lsv_list:
            if lsv.id == lsv_id:
                ret = lsv
                break
        else:
            try:
                ret = LSV(exon, lsv_id, jlist, lsv_type)
                self.lsv_list.append(ret)
            except InvalidLSV as e:
                if logger:
                    logger.info("Attempt to create LSV with wrong type or not enought junction coverage %s" %
                                exon.get_id())
                raise InvalidLSV(e.msg)

        # for jj in jlist:
        #     if jj.is_virtual() and logger:
        #         logger.info("WE FOUND INTRON RETENTION in exon %s" % exon.get_coordinates())

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

    def get_rnaseq_mat(self, rand10k, use_annot=True):

        ss3_l = []
        ss5_l = []
        tlb = {}
        exidx = 0
        ss_3p_vars = [0]*20
        ss_5p_vars = [0]*20
        ss_both_var = 0
        exon_list = []
        ex_list = self.get_exon_list()
        for ex in ex_list:
#            if ex.id is None:
#                continue
            l3 = len(set(ex.ss_3p_list))
            l5 = len(set(ex.ss_5p_list))
            if l3 == 0 or l5 == 0:
                continue

            local_3p, local_5p = ex.ss_variant_counts()
            if local_3p > 1 and local_5p > 1:
                ss_both_var += 1

            ss_3p_vars[local_3p] += 1
            ss_5p_vars[local_5p] += 1

            st3 = len(ss3_l)
            st5 = len(ss5_l)
            ss3_l += sorted([ss3 for ss3 in set(ex.ss_3p_list)])
            ss5_l += sorted([ss5 for ss5 in set(ex.ss_5p_list)])
            tlb[exidx] = [range(st3, len(ss3_l)), range(st5, len(ss5_l))]
            exon_list.append(ex)
            exidx += 1

        mat = np.zeros(shape=(len(ss5_l), len(ss3_l)), dtype='int')
        # jmat = np.empty(shape=(len(ss5_l), len(ss3_l)), dtype='object')
        # jmat.fill(None)

        junc_list = self.get_all_junctions()
        for junc in junc_list:
            st, end = junc.get_coordinates()
            if not st in ss5_l or not end in ss3_l:
                continue
            x = ss5_l.index(st)
            y = ss3_l.index(end)

            read_num = junc.get_coverage().sum()
            if read_num > 0:
                count_mat = read_num
            elif junc.is_annotated() and use_annot:
                count_mat = -1
            else:
                count_mat = 0
            mat[x, y] = count_mat
            # jmat[x, y] = junc

            for exp_idx in range(config.num_experiments):
                if reliable_in_data(junc, exp_idx):
                    rand10k[exp_idx].add(junc)

        return mat, exon_list, tlb, [ss_3p_vars, ss_5p_vars, ss_both_var]


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

        return config.gene_tlb[self.gene_id]

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

    def sort_in_list(self):
        # strand = self.get_gene().get_strand()
        # if strand == '+':
        #     isneg = False
        # else:
        #     isneg = True
#        self.exon_list.sort(reverse=isneg)
        self.exon_list.sort()
