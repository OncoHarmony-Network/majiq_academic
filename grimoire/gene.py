#!/usr/bin/python
import numpy as np
import grimoire.utils.utils as utils
import mglobals
from grimoire.exon import Exon, ExonTx, collapse_list_exons, print_list_exons
from grimoire.lsv import LSV


class Gene:
    __eq__ = lambda self, other: self.chromosome == other.chromosome and self.strand == other.strand and self.start < other.end and self.end > other.start
    __ne__ = lambda self, other: self.chromosome != other.chromosome or self.strand != other.strand or self.start >= other.end or self.end <= other.start
    __lt__ = lambda self, other: ( self.chromosome < other.chromosome  or ( self.chromosome == other.chromosome and  (self.end < other.start or (self.end > other.start and self.start < other.end and self.strand == '+' and other.strand == '-'))))
#    __le__ = lambda self, other: ( self.chromosome < other.chromosome  or ( self.chromosome == other.chromosome and  (self.end < other.start or (self.start == other.start and self.end == other.end and )))
    __gt__ = lambda self, other: ( self.chromosome > other.chromosome  or ( self.chromosome == other.chromosome and (self.start > other.end or (self.end > other.start and self.start < other.end and self.strand == '-' and other.strand == '+'))))
#    __ge__ = lambda self, other: ( self.chromosome > other.chromosome  or ( self.chromosome == other.chromosome and  (self.start > other.end and self.end <= other.start)))
#    __lt__ = lambda self, other: ( self.chromosome < other.chromosome  or ( self.chromosome == other.chromosome and  ((self.strand == '+' and other.strand == '-') or (self.strand == other.strand and self.end < other.start))))
#    __le__ = lambda self, other: ( self.chromosome < other.chromosome  or ( self.chromosome == other.chromosome and  ((self.strand == '+' and other.strand == '-') or (self.strand == other.strand and self.start <= other.end and self.end >= other.start))))
#    __gt__ = lambda self, other: ( self.chromosome > other.chromosome  or ( self.chromosome == other.chromosome and  ((self.strand == '-' and other.strand == '+') or (self.strand == other.strand and self.start > other.end))))
#    __ge__ = lambda self, other: ( self.chromosome > other.chromosome  or ( self.chromosome == other.chromosome and  ((self.strand == '-' and other.strand == '+') or (self.strand == other.strand and self.start > other.end and self.end <= other.start))))

    def __init__(self, gene_id, chrom, strand, start, end):
        self.id = gene_id
        self.chromosome = chrom
        self.strand = strand
        self.transcript_list = []
        self.exons = []
        self.start = start
        self.end = end
        self.otherNames = [gene_id]
        self.exonNum = 0
        self.readNum = np.zeros(shape=(mglobals.num_experiments), dtype=np.int)
        self.temp_txex_list = []
        self.ir_list = []
        self.lsv_list = []

    def __hash__(self):
        return hash((self.id,self.chromosome, self.strand, self.start, self.end))

    def get_id (self):
        return self.id

    def get_strand(self):
        return self.strand

    def get_RNAread_list(self):
#        print "GET RNAlist",len(self.RNAread_list[0])
        return self.RNAread_list

    def get_read_count(self):
        return self.readNum

    def get_RPKM(self):
        return self.RPKM

    def get_chromosome(self):
        return self.chromosome

    def get_coordinates(self):
        return self.start, self.end

    def get_transcript_AS_candidates(self):
        return self.transAScandidates

    def get_transcript_CONST_candidates(self):
        return self.transCONSTcandidates

    ''' Set functions '''

    # def add_transcript_AS_candidates(self,list_candidates):
    #     self.transAScandidates += list_candidates
    #     return
    #
    # def add_transcript_CONST_candidates(self,list_candidates):
    #     self.transCONSTcandidates += list_candidates
    #     return

    def add_transcript(self, tcrpt ):
        if tcrpt.txstart < self.start :
            self.start = tcrpt.txstart
        if tcrpt.txend > self.end :
            self.end = tcrpt.txend

        self.transcript_list.append(tcrpt)
        return
    
    def add_intron_retention(self, lsv_IR):
        self.ir_list.append(lsv_IR)


    def add_read_count(self, readNum,exp_idx):
        self.readNum[exp_idx] += readNum
        return

    def add_read(self, read, exp_idx ):
        self.readNum[exp_idx] += read.get_read_count()
#        (self.RNAread_list[exp_idx]).append(read)
        return

    def add_exon(self, exon):
        self.exons.append(exon)
        return

    def in_transcript_list(self,tcrpt_name):
        res = False
        for ff in self.transcript_list :
            if ff.get_id() == tcrpt_name:
                res = True
                break
        return res

    def is_gene_in_list(self, list, name):
        res = None
        for ll in list:
            if ( self.chromosome == ll.chromosome and self.strand == ll.strand and self.start < ll.end and self.end > ll.start):
                res = ll
                if not name in ll.otherNames:
                    ll.otherNames.append(name)
                ll.start = min(ll.start,self.start)
                ll.end = max(ll.end,self.end)
                break
        return res


    def calculate_RPKM( self, experiment_index, total_Reads ) :
        '''
         .. function: calculate_RPKM( self, experiment_index, total_Reads )

            This function calculates the RPKM as follows
            rpk = #Reads in Gene / # of kilobases of the gene exons (some may don't have any read)
            rpkm = rpk / (#total reads/1000000)

            :param experiment_index: Index of the experiment from the origrinal experiment list
            :param total_Reads: Total Number of reads of the gene.
            :rtype: RPKM value for this gene
        '''
        if len(self.exons) == 0 : return 0
        total_kb = float(0)
        for ex in self.exons:
            start, end = ex.get_coordinates()
#            print "EXON ::",ex.id,end, start
            total_kb += float(end-start)

#        print self.readNum, experiment_index
        rpk = float(self.readNum[experiment_index]) / float(total_kb/1000)
        mreads = float(total_Reads)/ float(1000000)
        rpkm = float(rpk)/mreads
        #print "Strand",self.strand,"::",total_kb, self.readNum, rpk, mreads, rpkm
        self.RPKM[experiment_index] = rpkm
        return rpkm

    def exist_exon (self, start, end):
        '''
         .. function: exist_exon (self, start, end):
            Check if the pair (start, end) are in a known exon in the gene. If not return None. We assume 
            that the given exon is for the same chromosome and strand than the gene.

            :param start: Start position of the input exon.
            :param end: End Position of te input exon.
            :rtype: exon instance or None
        '''

        res = None
        for ee in self.exons:
#            print "EX GEN:",ee.start, ee.end
#            print "New EX:",start, end
            if (start < ee.end and end > ee.start) :
                res = ee
                break
#                fnd +=1
#            elif (end < ee.start):
#                if fnd >1 : # overlapping more than one exon.... intron retention
#                    res = None
        return res

    def exist_junction(self,start,end):
        if start == None or end == None:
            return
        res = None
        for ff in self.transcript_list :
            res = ff.in_junction_list(start,end)
            if not res is None:
                break
        return res

    def get_all_junctions(self):
        lst = set()
        for ex in self.get_exon_list():
            for ex_rna in ex.exonRead_list:
                if len(ex_rna.p5_junc) > 0:
                    lst = lst.union(set(ex_rna.p5_junc))

        for tt in self.transcript_list:
            for jj in tt.get_junction_list():
                if not jj is None and not jj in lst :
                    lst.add(jj)
        s_junc = list(lst)
        return sorted([xx for xx in s_junc if not xx is None])

    def get_annotated_junctions(self):
        lst = set()
        for tt in self.transcript_list:
            for jj in tt.get_junction_list():
                if not jj is None and not jj in lst :
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

    def prepare_exons( self ) :
#        self.exons.sort(reverse = isneg)
        self.exons.sort()
        for idx,exs in enumerate(self.exons):
            exs.id = idx+1
        self.exonNum = len(self.exons)
        return

    def get_exon_list(self):
        return self.exons

    def get_all_ss(self, anot_only = False ):

        ss = set()
        for ex in self.exons:
            tx_list = ex.get_annotated_exon()
            for txex in tx_list:
                for junc in txex.get_3prime_junc():
                    ss.add((junc.get_ss_3p(),'3prime',junc))
                for junc in txex.get_5prime_junc():
                    ss.add((junc.get_ss_5p(),'5prime',junc))

        print sorted(ss)
        return sorted(ss)

    def collapse_exons(self):

        self.temp_txex_list.sort()
#        print "GEN",self.temp_txex_list
#        print_list_exons(self.temp_txex_list)
        list_ex = collapse_list_exons(self.temp_txex_list, self)
        self.exons.extend(list_ex)
        self.prepare_exons()
        self.remove_temp_attributes()

    def check_exons(self):

        s_exons = set()

        for ex in self.exons:
            s_exons.add(ex.get_coordinates())
            
        assert len(s_exons) == len(self.exons), "Exist duplicates in exons in Gene %s"%(self.id)

    def new_lsv_definition(self, exon, jlist, lsv_type):

        coords = exon.get_coordinates()
        ret = None
        id = "%s:%d-%d:%s" % (self.get_id(), coords[0], coords[1], lsv_type)
        for lsv in self.lsv_list:
            if lsv.id == id:
                ret = lsv
                break
        else:
            try:
                ret = LSV(exon, id, jlist, lsv_type)
                self.lsv_list.append(ret)
            except ValueError:
                print "Attempt to create LSV with wrong type or not enought junction coverage"

        return ret

    def remove_temp_attributes(self):
        del self.temp_txex_list

    def new_annotated_exon(self, start, end, transcript, bl=True):
        for txex in self.temp_txex_list:
            if txex.start == start and txex.end == end:
                res = txex
                break
        else:
            res = ExonTx(start, end, transcript, None)
            if bl:
                self.temp_txex_list.append(res)
        return res

    # mark for delete
    def get_transcript_mat(self, ASvsConst):

        if (len(self.transcript_list) == 1 and ASvsConst == 'AS') or len(self.exons) <3:
            return None

        mat = np.ndarray(shape=(len(self.transcript_list), len(self.exons)), dtype='bool')
        for idx_t, tpt in enumerate(self.transcript_list):
            for g_ex in self.exons:
                if set(tpt.exon_list).intersection(set(g_ex.exonTx_list)):
                    mat[idx_t, g_ex.id-1] = 1
                else:
                    mat[idx_t, g_ex.id-1] = 0

        return mat


    def get_rnaseq_mat(self, rand10k, lsv=False):

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
            if ex.id is None: continue
            l3 = len(set(ex.ss_3p_list))
            l5 = len(set(ex.ss_5p_list))
            if l3 == 0 or l5 == 0: continue

            local_3p, local_5p = ex.ss_variant_counts()
            if local_3p > 1 and local_5p >1 : ss_both_var += 1

            ss_3p_vars[local_3p] += 1
            ss_5p_vars[local_5p] += 1

            st3 = len(ss3_l)
            st5 = len(ss5_l)
            ss3_l += sorted([ss3 for ss3 in set(ex.ss_3p_list)])
            ss5_l += sorted([ss5 for ss5 in set(ex.ss_5p_list)])
            tlb[exidx] = [range(st3,len(ss3_l)),range(st5,len(ss5_l))]
            exon_list.append(ex)
            exidx += 1

        mat  = np.empty(shape=(len(ss5_l),len(ss3_l)),dtype='int')
        jmat = np.empty(shape=(len(ss5_l),len(ss3_l)),dtype='object')
        mat.fill(0)
        jmat.fill(None)

        junc_list = self.get_all_junctions()
        for junc in junc_list:
            st,end = junc.get_coordinates()
            if not st in ss5_l or not end in ss3_l: continue
            x = ss5_l.index(st)
            y = ss3_l.index(end)

            if junc.readN.sum() >0:
                count_mat = junc.readN.sum()
            elif junc.is_annotated() and lsv :
                count_mat = -1
            else:
                count_mat = 0
            mat [ x, y ] = count_mat
            jmat[ x, y ] = junc

            for exp_idx in range(mglobals.num_experiments):
                if junc.get_readN(exp_idx) >= 10 : #and in_DB:
                    rand10k[exp_idx].add(junc)

        if not lsv :
            return mat, jmat, tlb, [ss_3p_vars, ss_5p_vars,ss_both_var]
        else:
            return mat, exon_list, tlb, [ss_3p_vars, ss_5p_vars,ss_both_var]


class Transcript :

    def __init__( self, name, gene,txstart,txend ):
        self.id = name
        self.gene = gene
        self.exon_list = []
        self.junction_list = []
        self.txstart = txstart
        self.txend = txend
        self.cdsStart = None
        self.cdsStop = None

    def add_exon(self, exon):
        self.exon_list.append(exon)
        return

    def get_gene(self):
        return self.gene

    def get_id(self):
        return self.id

    def get_junction_list(self):
        return self.junction_list

    def in_junction_list(self,start,end):
        res = None 
#        if self.gene.strand == '-':
#            tmp = start
#            start = end
#            end = tmp
        for jj in self.junction_list:
            if jj.start == start and jj.end == end:
#                jj.txN += 1
                res = jj
                break
        return res

    def add_junction (self, junc):
        if junc not in self.junction_list:
            self.junction_list.append(junc)

    def sort_in_list(self):
        strand = self.gene.get_strand()
        if strand == '+':
            isneg=False
        else:
            isneg = True
        self.junction_list.sort()
#        self.exon_list.sort(reverse=isneg)
        self.exon_list.sort()

