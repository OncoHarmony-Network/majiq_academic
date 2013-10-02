#!/usr/bin/python
import numpy as np

import my_globals

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

    def __init__(self,gene_id,chrom,strand,start,end):
        self.id = gene_id
        self.chromosome = chrom
        self.strand = strand
        self.transcript_list = []
        self.exons = []
        self.start = start
        self.end = end
        self.otherNames = [gene_id]
#        self.RNAread_list = np.zeros(shape=(globals.num_experiments),dtype=np.dtype('object'))
#        self.RNAread_list.fill([])
        self.exonNum = 0
        self.readNum = np.zeros(shape=(my_globals.num_experiments),dtype=np.int)
        self.RPKM = np.zeros(shape=(my_globals.num_experiments),dtype=np.float)

        self.transAScandidates = []
        self.transCONSTcandidates = []

    def __hash__(self):
        return hash((self.id,self.chromosome, self.strand, self.start, self.end))

    def get_id (self):
        return self.id

    def get_strand(self):
        return self.strand

    def get_RNAread_list(self):
        return self.RNAread_list

    def get_read_count(self):
        return self.readNum

    def get_RPKM(self):
        return self.RPKM

    def get_chromosome(self):
        return self.chromosome

    def get_coordinates( self ):
        return (self.start,self.end)

    def get_transcript_AS_candidates(self):
        return (self.transAScandidates)

    def get_transcript_CONST_candidates(self):
        return (self.transCONSTcandidates)

    ''' Set functions '''

    def add_transcript_AS_candidates(self,list_candidates):
        self.transAScandidates += list_candidates
        return

    def add_transcript_CONST_candidates(self,list_candidates):
        self.transCONSTcandidates += list_candidates
        return

    def add_transcript(self, tcrpt ):
        if tcrpt.txstart < self.start :
            self.start = tcrpt.txstart
        if tcrpt.txend > self.end :
            self.end = tcrpt.txend

        (self.transcript_list).append(tcrpt)
        return
    
    def add_read_count(self, readNum,exp_idx):
        self.readNum[exp_idx] += readNum
        return

    def add_read(self, read, exp_idx ):
        self.readNum[exp_idx] += read.get_read_count()
#        (self.RNAread_list[exp_idx]).append(read)
        return

    def add_exon(self, exon ):
        (self.exons).append(exon)
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
                ll.end = min(ll.end,self.end)
                break
        return res

    def calculate_RPKM( self, experiment_index, total_Reads ) :
        '''
        This function calculates the RPKM as follows
        rpk = #Reads in Gene / # of kilobases of the gene exons (some may don't have any read)
        rpkm = rpk / (#total reads/1000000)
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
        res = None
        for ee in self.exons:
            if (start < ee.end and end > ee.start) :
                res = ee
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
            if res is None:
                break
        return res

    def get_all_junctions(self):
        lst = set()
        for ex in self.get_exon_list():
            for ex_rna in ex.exonRead_list:
                if not ex_rna.p5_junc is None:
                    lst.add(ex_rna.p5_junc)
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
                if not ex_rna.p5_junc is None:
                    lst.append(ex_rna.p5_junc)
        lst.sort()
        return lst

    def prepare_exons( self ) :
#        self.exons.sort(reverse = isneg)
        self.exons.sort()
        for idx,exs in enumerate(self.exons):
            exs.id = idx+1
        self.exonNum = len(self.exons)
        return

    def get_exon_list( self ):
        return self.exons

    def get_all_ss( self ):

        ss = {}
        for ex in self.exons:
            ss3_l = [(ss3,'3prime') for ss3 in set(ex.ss_3p_list)]
            ss5_l = [(ss5,'5prime') for ss5 in set(ex.ss_5p_list)]
            name = ex.get_id()
            if name in ss:
                print "ERROR id is already in "
            ss[name] = sorted(list(ss3_l)+list(ss5_l))
        return ss


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

    def get_id(self):
        return self.id

    def get_junction_list(self):
        return self.junction_list

    def in_junction_list(self,start,end):
        res = None 
        if self.gene.strand == '-':
            tmp = start
            start = end
            end = tmp
        for jj in self.junction_list:
            if jj.start == start and jj.end == end:
                jj.txN += 1
                res = jj
                break
        return res

    def add_junction (self, junc):
        self.junction_list.append(junc)

    def _sort_in_list(self,strand):
        if strand == '+':
            isneg=False
        else:
            isneg = True
        self.junction_list.sort()
#        self.exon_list.sort(reverse=isneg)
        self.exon_list.sort()

