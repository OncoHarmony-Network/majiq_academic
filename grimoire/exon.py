#!/usr/bin/python

import numpy as np
import mglobals

class Exon:


#    __eq__ = lambda self, other: self.start == other.start  and self.end == other.end
#    __ne__ = lambda self, other: self.start != other.start  or  self.end != other.end
#    __lt__ = lambda self, other: self.start < other.start   or (self.start == other.start and self.end < other.end)
#    __le__ = lambda self, other: self.start <= other.start  or (self.start == other.start and self.end <= other.end)
#    __gt__ = lambda self, other: self.start > other.start   or (self.start == other.start and self.end > other.end)
#    __ge__ = lambda self, other: self.start >= other.start  or (self.start == other.start and self.end >= other.end)
    
    __eq__ = lambda self, other: self.start < other.end and  self.end > other.start
    __ne__ = lambda self, other: self.start >= other.end or   self.end <= other.start
    __lt__ = lambda self, other: self.end <= other.start 
   # or (self.end>=other.start and (self.start<other.start or (self.start==other.start and self.end<other.end)))
    __le__ = lambda self, other: self.end < other.start or (self.start < other.end and self.end > other.start)
    __gt__ = lambda self, other: self.start >= other.end 
#or (self.end>=other.start and (self.start>other.start or (self.start==other.start and self.end>other.end)))
    __ge__ = lambda self, other: self.start >= other.end  or (self.start < other.end and self.end > other.start)

    def __init__ ( self, start, end, gene, strand ):
#        print "Creating Exon",start,end, self
        self.start = start
        self.end = end
        self.gene = gene
        self.exonTx_list = []
        self.exonRead_list = []
        self.ss_3p_list = []
        self.ss_5p_list = []
        self.id = None
        self.strand = strand
        self.gc_content = None
        self.coverage = np.zeros(shape=(mglobals.num_experiments))

    def __hash__(self):
        return hash(self.id) ^ hash(self.gene.id)

    def get_id(self):
        return self.id

    def get_strand( self ):
        return self.strand

    def get_coordinates( self ):
        '''
         .. function:: get_coordinates(self)
            Get the exon start and end.
            :rtype: tuple of 2 integers, (exon start, exon end)
        '''
        return (self.start, self.end)

    def get_gene(self):
        return self.gene

    def get_exon_definition(self, transcript):
        res = None
        for ii in self.exonTx_list:
            #print transcript, ii.transcript
            if transcript in ii.transcript:
                res = ii
                break
        if res ==  None:
            print "error"
        return res

    def get_rna_ss(self):
        ss3 = set()
        ss5 = set()

        for ii in self.exonRead_list:
            ss3.add(ii.start)
            ss5.add(ii.end)
        ss3_l = sorted(list(ss3))
        ss5_l = sorted(list(ss5))
        return (ss3_l,ss5_l)

    def add_new_read(self, start, end, readSeq, s3p_junc, s5p_junc):

        if start > end : 
            print " INCORRECT exon definition",start, end
            exit()
        if start < self.start :
            self.start = start
        if end > self.end :
            self.end = end
        found = False
        for ii in self.exonRead_list:
            if start == ii.start and end == ii.end:
                res = ii
                self.exonRead_list.append(res)
                found = True
                break
        if not found :
#            if self.strand == '+':
            ssf = False
            for idx1, i1 in enumerate(self.ss_3p_list):
                if i1 != start : continue
                if self.ss_5p_list[idx1]== end:
                    ssf = True
                    break
            else:
                self.ss_3p_list.append(start)
                self.ss_5p_list.append(end)
            res = ExonRead( start, end, self, s3p_junc, s5p_junc, readSeq )
            self.exonRead_list.append(res)
        return res

    def add_new_definition ( self,start,end,trncpt ):
        if start < self.start :
            self.start = start
        if end > self.end :
            self.end = end
        found = False
        for ii in self.exonTx_list:
            if start == ii.start and end == ii.end:
                res = ii
                ii.add_transcript(trncpt)
                found = True
                break
        if start == 74289871: print "NOUUUUNN", found
        if not found :
#            if self.strand == '+':
            self.ss_3p_list.append(start)
            self.ss_5p_list.append(end)
#            else:
#                self.ss_5p_list.append(start)
#                self.ss_3p_list.append(end)
            res = ExonTx(start,end,trncpt,self)
            self.exonTx_list.append(res)
        return res

    def get_coverage( self,exp_idx ):
        return self.coverage[exp_idx]

    def get_gc_content( self ):
        return self.gc_content

    def update_coverage( self, exp_idx, num ):
        self.coverage[exp_idx] += num

    def set_gc_content(self, sequence):
#        if len(self.exonTx_list) != 0 and len(self.exonRead_list) != 0 :
        cs = sequence.count('c') + sequence.count('C')
        gs = sequence.count('g') + sequence.count('G')
        if len(sequence) == 0 : return
        self.gc_content = float( cs + gs) / float(len(sequence))
            

    def print_triplet_coord(self, fp):
        gene = self.gene
        chr = gene.chromosome
        strand = gene.get_strand()
        startA = self.start
        endA = self.end

        if strand == '-':
            idA = gene.exonNum - (self.id - 1)
            vidC1 = self.id
            vidC2 = self.id - 2
        else:
            idA = self.id
            vidC1 = self.id - 2
            vidC2 = self.id

        startC1,endC1 = gene.exons[vidC1].get_coordinates()
        startC2,endC2 = gene.exons[vidC2].get_coordinates()

        name = "%s.%s"%(gene.id,idA)
        
        fp.write("%s\t%d\t%d\t%s_C1\t0\t%s\n"%(chr,startC1,endC1,name,strand))
        fp.write("%s\t%d\t%d\t%s_A\t0\t%s\n"%(chr,startA,endA,name,strand))
        fp.write("%s\t%d\t%d\t%s_C2\t0\t%s\n"%(chr,startC2,endC2,name,strand))

    def bed_format(self):
        str = ""
        for eRead in self.exonRead_list:
            str += "%s\n"%eRead.bed_format()

        return str

class ExonRead(object):
    
    def __init__ ( self, start, end, exon,pre_junc, post_junc,rnaSeq=None):
        self.start = start
        self.end = end
        self.exon = exon
        self.RNASeq = rnaSeq
        self.p3_junc = pre_junc
        self.p5_junc = post_junc

    def get_coordinates( self ):
        return ( self.start, self.end )

    def get_3p_exon(self):
        return self.p5_junc.get_acceptor()

    def get_5p_junc(self):
        return self.p5_junc

    def bed_format(self):
        chr = self.exon.get_gene().get_chromosome() 
        str = "%s\t%s\t%s\t.\t0\t.\t%s\t%s"%(chr,self.start,self.end,self.start,self.end)
        return str

class ExonTx(object):
    def __init__ ( self, start, end, trnscpt, exon):
        self.start = start
        self.end = end
        self.transcript = [trnscpt]
        #self.p3_junc = pre_junc
        #self.p5_junc = post_junc
        self.exon = exon

    def get_coordinates( self ):
        return ( self.start, self.end )

    def add_transcript( self, trans ):
        self.transcript.append(trans)
        return



def __half_exon(type,junc,readRNA):

    gene = junc.get_gene()
    if type == '3prime':
        coord = junc.get_ss_3p()
    else:
        coord = junc.get_ss_5p()

    for ex in gene.get_exon_list():
        (ex_start,ex_end) = ex.get_coordinates()
        if ex_start <= coord and ex_end >= coord:
            if type == '3prime':
                start = coord
                end = ex_end
                junc.add_acceptor(ex)
                to = junc
                frm = None
            else:
                start = ex_start
                end = coord
                junc.add_donor(ex)
                to = None
                frm = junc
            #print "half",type,"::",ex_start, ex_end, junc.start, junc.end, end 
            if end - start < 10 : continue
            ex.add_new_read( start, end, readRNA, to, frm )

            break
    return 0

def new_exon_definition(start, end, readRNA, s3prime_junc, s5prime_junc, gene):

    if  end - start < 10 : return 0

    ex = gene.exist_exon(start,end)
    newExons = 0
    if ex is None :
        newExons = 1
        #print "CREATE1::",start, end
        ex = Exon(start,end,gene,gene.get_strand())
        gene.add_exon(ex)
    #print "ADD1::",start, end
    ex.add_new_read( start, end, readRNA, s3prime_junc, s5prime_junc )
    s3prime_junc.add_acceptor(ex)
    s5prime_junc.add_donor(ex)

    return newExons


def detect_exons(gene_list, junction_list, readRNA):

    newExons = 0
    for strand, j_list in junction_list.items():
        opened = 0
        opened_exon = []
        last_5prime = None
        first_3prime = None
        for (coord,type, jj) in j_list :
            jj_gene = jj.get_gene()
            if type == '5prime':
                if opened >0 :
                    prev_gene = opened_exon[-1].get_gene()
                    if jj_gene != prev_gene :
                        ''' New Gene will start with 5prime SS. 
                        But there are non-paired 3'ss elements.'''
                        newExons +=  __half_exon('3prime',opened_exon[-1],readRNA)
                        opened = 0
                        opened_exon = []
                        last_5prime = None
                        newExons += __half_exon('5prime',jj,readRNA)
                    else:
                        start = opened_exon[-1].get_ss_3p()
                        end = coord
                        newExons += new_exon_definition(start,end,readRNA, opened_exon[-1],jj, jj_gene)
                        last_5prime = jj
                        opened_exon.pop()
                        opened -=1
                elif opened == 0 :
                    newExons += __half_exon('5prime',jj,readRNA)
                #end elif opened
            else:
                if opened >0 :
                    if not last_5prime is None:
                        end = last_5prime.get_ss_5p()
                        for ss in opened_exon:
                            if ss.get_gene() != last_5prime.get_gene(): continue
                            start = ss.get_ss_3p()
                            newExons += new_exon_definition(start,end,readRNA, ss, last_5prime, ss.get_gene())
                    else:
                        for ss in opened_exon:
                            newExons += __half_exon('3prime', ss, readRNA)
                    #end else ...
                    opened = 0
                    opened_exon = []
                else:
                    first_3prime = jj
                last_5prime = None
                opened_exon.append(jj)
                opened +=1

    print "FOUND new %d exons"%newExons
    return 

