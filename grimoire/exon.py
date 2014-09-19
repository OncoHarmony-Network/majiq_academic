#!/usr/bin/python

import numpy as np
import mglobals
from grimoire.lsv import LSV_IR
from grimoire.junction import Junction


class Exon:

    __eq__ = lambda self, other: self.start < other.end and self.end > other.start
    __ne__ = lambda self, other: self.start >= other.end or self.end <= other.start
    __lt__ = lambda self, other: self.end <= other.start 
    __le__ = lambda self, other: self.end < other.start or (self.start < other.end and self.end > other.start)
    __gt__ = lambda self, other: self.start >= other.end 
    __ge__ = lambda self, other: self.start >= other.end or (self.start < other.end and self.end > other.start)

    def __init__(self, start, end, gene, strand, annot=False):

        self.start = start
        self.end = end
        self.gene_name = gene.get_id()
        self.exonTx_list = []
        self.exonRead_list = []
        self.ss_3p_list = []
        self.ss_5p_list = []
        self.id = None
        self.strand = strand
        self.gc_content = None
        self.coverage = np.zeros(shape=mglobals.num_experiments)
        self.score = None
        self.pcr_name = None
        self.pcr_candidate = None
        self.ir = False
        self.db_coord = (start, end)
        self.annotated = annot

    def __hash__(self):
        return hash(self.id) ^ hash(self.gene_name)

    def get_id(self):
        return self.id

    def get_strand(self):
        return self.strand

    def get_coordinates(self):
        """
         .. function:: get_coordinates(self)
            Get the exon start and end.
            :rtype: tuple of 2 integers, (exon start, exon end)
        """
        return self.start, self.end

    def get_gene(self):
        return mglobals.gene_tlb[self.gene_name]

    def get_annotated_exon(self):
        return self.exonTx_list

    def get_exon_definition(self, transcript):
        res = None
        for ii in self.exonTx_list:
            #print transcript, ii.transcript
            if transcript in ii.transcript:
                res = ii
                break
        if res is None:
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
        return ss3_l, ss5_l

    def set_ir(self, ir):
        self.ir |= ir

    def set_pcr_score(self, pcr_name, score, candidate):
        self.pcr_name = pcr_name
        self.score = score
        self.pcr_candidate = candidate

    def get_pcr_score(self):
        return self.score

    def get_pcr_candidate(self):
        return self.pcr_candidate

    def get_pcr_name(self):
        return self.pcr_name

    def add_new_read(self, start, end, read_seq, s3p_junc, s5p_junc):

        #assert start < end , " INCORRECT exon definition %s - %s "%(start, end)
        if start >= end:
            return None
        if start < self.start:
            self.start = start
        if end > self.end:
            self.end = end
        found = False
        for ii in self.exonRead_list:
            if start == ii.start and end == ii.end:
                res = ii
                #self.exonRead_list.append(res)
                found = True
                break
        if not found:
            ssf = False
            for idx1, i1 in enumerate(self.ss_3p_list):
                if i1 != start:
                    continue
                if self.ss_5p_list[idx1] == end:
                    ssf = True
                    break
            else:
                self.ss_3p_list.append(start)
                self.ss_5p_list.append(end)
            res = ExonRead(start, end, self, s3p_junc, s5p_junc, read_seq)
            self.exonRead_list.append(res)
        return res

    def add_new_definition(self, start, end, trncpt):
        if start < self.start:
            self.start = start
        if end > self.end:
            self.end = end
        found = False
        for ii in self.exonTx_list:
            if start == ii.start and end == ii.end:
                res = ii
                # ii.add_transcript(trncpt)
                found = True
                break
        if not found:
#            if self.strand == '+':
            self.ss_3p_list.append(start)
            self.ss_5p_list.append(end)
#            else:
#                self.ss_5p_list.append(start)
#                self.ss_3p_list.append(end)
            res = ExonTx(start, end, trncpt, self)
            self.exonTx_list.append(res)
        return res

    def get_coverage(self, exp_idx):
        return self.coverage[exp_idx]

    def get_gc_content(self):
        return self.gc_content

    def update_coverage( self, exp_idx, num ):
        self.coverage[exp_idx] += num

    def set_gc_content(self, sequence):
#        if len(self.exonTx_list) != 0 and len(self.exonRead_list) != 0 :
        cs = sequence.count('c') + sequence.count('C')
        gs = sequence.count('g') + sequence.count('G')
        if len(sequence) == 0:
            return
        self.gc_content = float(cs + gs) / float(len(sequence))

    def get_junctions(self, jtype):
        if jtype != '3prime' and jtype != '5prime':
            raise RuntimeError('Incorrect splicesite type %s' % jtype)
        jlist = set()
        for exon_list in (self.exonTx_list, self.exonRead_list):
            for ex in exon_list:
                if jtype == '3prime':
                    for junc in ex.p3_junc:
                        jlist.add(junc)
                else:
                    for junc in ex.p5_junc:
                        jlist.add(junc)
        return jlist

    def print_triplet_coord(self, fp):
        gene = mglobals.gene_tlb[self.gene_name]
        chrom = gene.chromosome
        strand = gene.get_strand()
        start_a = self.start
        end_a = self.end

        if strand == '-':
            id_a = gene.exonNum - (self.id - 1)
            vid_c1 = self.id
            vid_c2 = self.id - 2
        else:
            id_a = self.id
            vid_c1 = self.id - 2
            vid_c2 = self.id

        start_c1, end_c1 = gene.exons[vid_c1].get_coordinates()
        start_c2, end_c2 = gene.exons[vid_c2].get_coordinates()

        name = "%s.%s" % (gene.id, id_a)
        
        fp.write("%s\t%d\t%d\t%s_C1\t0\t%s\n" % (chrom, start_c1, end_c1, name, strand))
        fp.write("%s\t%d\t%d\t%s_A\t0\t%s\n" % (chrom, start_a, end_a, name, strand))
        fp.write("%s\t%d\t%d\t%s_C2\t0\t%s\n" % (chrom, start_c2, end_c2, name, strand))

    def bed_format(self):
        bed_str = ""
        for eRead in self.exonRead_list:
            bed_str += "%s\n" % eRead.bed_format()

        return bed_str

    def ss_variant_counts(self, minreads=5):
        local_3p = 0
        local_5p = 0

        temp_set = set()
        for ss3p in self.ss_3p_list:
            for exread in self.exonRead_list:
                if ss3p != exread.start:
                    continue
                sum_reads = 0
                for jj in exread.p3_junc:
                    if jj is None:
                        continue
                    sum_reads += jj.readN.sum()
                if sum_reads >= minreads:
                    temp_set.add(ss3p)
        local_3p = len(temp_set)

        temp_set = set()
        for ss5p in self.ss_5p_list:
            for exread in self.exonRead_list:
                if ss5p != exread.end:
                    continue
                sum_reads = 0
                for jj in exread.p3_junc:
                    if jj is None:
                        continue
                    sum_reads += jj.readN.sum()
                if sum_reads >= minreads:
                    temp_set.add(ss5p)
        local_5p = len(temp_set)

        if local_3p > 19:
            local_3p = 19
        if local_5p > 19:
            local_5p = 19

        return local_3p, local_5p


class ExonRead(object):
    
    def __init__(self, start, end, exon, pre_junc, post_junc, rna_seq=None):
        self.start = start
        self.end = end
#        self.exon = exon
        self.RNASeq = rna_seq
        self.p3_junc = [pre_junc]
        self.p5_junc = [post_junc]

    def get_coordinates(self):
        return self.start, self.end

    def get_5p_junc(self):
        return self.p5_junc

    def bed_format(self):
        chrom = self.exon.get_gene.get_chromosome()
        strng = "%s\t%s\t%s\t.\t0\t.\t%s\t%s" % (chrom, self.start, self.end, self.start, self.end)
        return strng


class ExonTx(object):

    __eq__ = lambda self, other: self.start == other.start and self.end == other.end
    __ne__ = lambda self, other: self.start != other.start or self.end != other.end
    __lt__ = lambda self, other: self.start < other.start or (self.start == other.start and self.end < other.end)
    __le__ = lambda self, other: self.start <= other.start or (self.start == other.start and self.end <= other.end)
    __gt__ = lambda self, other: self.start > other.start or (self.start == other.start and self.end > other.end)
    __ge__ = lambda self, other: self.start >= other.start or (self.start == other.start and self.end >= other.end)

    def __init__(self, start, end, trnscpt, exon):
        self.start = start
        self.end = end
        #self.transcript = [trnscpt]
        self.transcript_name = [trnscpt.get_id()]
        self.gene_name = trnscpt.get_gene().get_id()
        self.p3_junc = []
        self.p5_junc = []
#        self.exon = exon
        self.ir = False

    def get_coordinates(self):
        return self.start, self.end

    def add_transcript(self, trans):
        self.transcript.append(trans)

    def add_5prime_junc(self, p5_junc):
        if p5_junc not in self.p5_junc:
            self.p5_junc.append(p5_junc)

    def add_3prime_junc(self, p3_junc):
        if p3_junc not in self.p3_junc:
            self.p3_junc.append(p3_junc)

    def get_transcript(self):
        res = []
        for tx_name in self.transcript_name:
            res.append(mglobals.gene_tlb[self.gene_name].get_transcript(tx_name))
        return res

    def get_5prime_junc(self):
        return self.p5_junc

    def get_3prime_junc(self):
        return self.p3_junc

    def overlaps(self, start, end):
        res = False
        if self.start < end and self.end > start:
            res = True
        return res

    def split_exon(self, intron_coords, gn):

        res = []
        exb1 = False
        exb2 = False
        if self.end - intron_coords[1]+1 > 5:
            txex1 = gn.new_annotated_exon(intron_coords[1]+1, self.end, self.get_transcript()[0], bl=False)
            txex1.p5_junc.extend(self.p5_junc)
            res.append(txex1)
            txex1.junction_consistency()
            exb1 = True
        if intron_coords[0] - 1 - self.start > 5:
            txex2 = gn.new_annotated_exon(self.start, intron_coords[0]-1, self.get_transcript()[0], bl=False)
            txex2.p3_junc.extend(self.p3_junc)
            exb2 = True
            res.append(txex2)
            txex2.junction_consistency()

        exb = exb1 & exb2

        if exb:
            junc = gn.exist_junction(txex2.end, txex1.start)
            if junc is None:
                junc = Junction(txex2.end, txex1.start, None, None, gn, annotated=True)
#                junc.add_donor(txex2)
#                junc.add_acceptor(txex1)
            txex2.p5_junc.append(junc)
            txex1.p3_junc.append(junc)

        for trn in self.get_transcript():
            if exb:
                trn.add_junction(junc)
            # if exb1:
            #     txex1.add_transcript(trn)
            # if exb2:
            #     txex2.add_transcript(trn)

        del self
        return res

    def junction_consistency(self):

        #print "EXONTX:", self.start, self.end

        j5_list = []
        for j5 in self.p5_junc:
            jcoord = j5.get_coordinates()
            #print 'P5',j5.get_gene().get_id(), j5.get_coordinates(), self.start, self.end
            if self.start <= jcoord[1] <= self.end:
                #j5.add_donor(self)
                j5_list.append(j5)
                
        j3_list = []
        for j3 in self.p3_junc:
            jcoord = j3.get_coordinates()
            #print 'P3::',j3.get_gene().get_id(), j3.get_coordinates(), self.start, self.end
            if self.start <= jcoord[1] <= self.end:
                #j3.add_acceptor(self)
                j3_list.append(j3)

        self.p3_junc = j3_list[:]
        self.p5_junc = j5_list[:]

        return

    def collapse(self, list_exontx, gne):
    
        all_5prime = [xx.end for xx in list_exontx]
        all_3prime = [xx.start for xx in list_exontx]
        all_5prime = sorted(set(all_5prime))
        all_3prime = sorted(set(all_3prime))
        exlist = []

        if max(all_3prime) > min(all_5prime):
            ''' Intron retention '''
            introns = []
            last_p5 = 0
            jdx = 0
            in_found = False
            for idx, p5 in enumerate(all_5prime):
                while jdx < len(all_3prime):
                    p3 = all_3prime[jdx]
                    if p3 < p5:
                        if in_found:
                            introns.append((last_p5+1, max(p3-1, last_p5+1)))
                            in_found = False
                        jdx += 1
                    else:
                        last_p5 = p5
                        in_found = True
                        break

            for idx, txex in enumerate(list_exontx):
                for intr in introns:
                    if not txex.overlaps(intr[0], intr[1]):
                        if intr[0] > txex.end or (intr[0] <= txex.end < intr[1]):
                            txex.ir = True
                        continue
                    ''' intron retention'''
                    LSV_IR(txex.start, txex.end, [], gne)
                    dummy = txex.split_exon(intr, gne)
                    list_exontx.remove(txex)
                    for dm in dummy:
                        if dm not in list_exontx:
                            list_exontx.append(dm)
                    break

            list_exontx.sort()
            exlist.extend(collapse_list_exons(list_exontx, gne))

        else:
            ex = Exon(min(all_3prime), max(all_5prime), gne, gne.get_strand(), annot=True)

            for txex in list_exontx:
                ex.set_ir(txex.ir)
                ex.ss_3p_list.append(txex.start)
                ex.ss_5p_list.append(txex.end)
                ex.exonTx_list.append(txex)
                txex.exon = ex
                for p3_junc in txex.p3_junc:
                    p3_junc.add_acceptor(ex)
                for p5_junc in txex.p5_junc:
                    p5_junc.add_donor(ex)
            exlist.append(ex)
        return exlist


def print_list_exons(list_ex, msg=""):
    #list_ex.sort()
    print "%%%%%%%%%%%%LIST_EXONS %s" % msg
    for ex in list_ex:
        print "\t\t", ex.get_coordinates(), ex
    print "%%%%%%%%%%%%%%%%%%%%%%"

num_it = 0


def collapse_list_exons(listexons, gne):
    global num_it
    num_it += 1
    #print "[%s] INIT COLLAPSE_LIST EXONS "%(num_it)
    #print_list_exons(listexons,"[%s] IN INIT"%num_it)
    overlp = []
    exlist = []
    start = 0
    end = 0
    for idx, ex in enumerate(listexons):
        if ex.overlaps(start, end):
            if start > ex.start:
                start = ex.start
            if end < ex.end:
                end = ex.end
            overlp.append(ex)
#            continue
        else:
            if len(overlp) > 0:
                exlist.extend(ex.collapse(overlp, gne))

            overlp = [ex]
            start, end = ex.get_coordinates()
        if idx == len(listexons)-1:
            exlist.extend(ex.collapse(overlp, gne))
    num_it -= 1
    return exlist


def __half_exon(ss_type, junc, read_rna):
    gene = junc.get_gene
    if ss_type == '3prime':
        coord = junc.get_ss_3p()
    else:
        coord = junc.get_ss_5p()

    for ex in gene.get_exon_list():
        (ex_start, ex_end) = ex.get_coordinates()
        if ex_start <= coord <= ex_end:
            if ss_type == '3prime':
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
#            if end - start < 10 : continue
            res = ex.add_new_read(start, end, read_rna, to, frm)
            if res:
                ex.ss_3p_list.append(start)
                ex.ss_5p_list.append(end)

            break
    return 0


def new_exon_definition(start, end, read_rna, s3prime_junc, s5prime_junc, gene):

    if end - start < 5:
        return 0
#    print "NEW DEFINITION::",start, end

    ex = gene.exist_exon(start, end)
    new_exons = 0
    if ex is None:
        new_exons = 1
        ex = Exon(start, end, gene, gene.get_strand())
        gene.add_exon(ex)
#    else:
#        print "EXON FOUND", ex, ex.get_coordinates(), ex.annotated
    ex.add_new_read(start, end, read_rna, s3prime_junc, s5prime_junc)
    s3prime_junc.add_acceptor(ex)
    s5prime_junc.add_donor(ex)

    return new_exons


def detect_exons(gene, junction_list, read_rna):

    new_exons = 0
    opened = 0
    opened_exon = []
    last_5prime = None
    first_3prime = None

#    for jj in  gene.get_annotated_junctions():
#        if not (jj.get_ss_5p(),'5prime',jj) in junction_list:
#            junction_list.append((jj.get_ss_5p(),'5prime',jj))
#        if not (jj.get_ss_3p(),'3prime',jj) in junction_list:
#            junction_list.append((jj.get_ss_3p(),'3prime',jj))

    junction_list.extend(gene.get_all_ss())

    junction_list.sort()
#    print "JUNC EXTENDED", junction_list
#    print "DETECT EXONS::",gene.get_id()
    for (coord, jtype, jj) in junction_list:
#        print "---NEW-------------------------------------------------------------"
#        print coord, jtype, jj, jj.coverage.sum(), jj.annotated, jj.is_annotated()
        if jj.coverage.sum() < mglobals.MINREADS and not jj.is_annotated():
            continue
#        print coord, jtype, jj, jj.coverage.sum(), jj.annotated
#        print "LIST",opened_exon
#        print "LASTS",first_3prime, last_5prime
        jj_gene = jj.get_gene()
        if jtype == '5prime':
#            print "CHECK 1",coord,jj.get_ss_5p()
            if opened > 0:
                start = opened_exon[-1].get_ss_3p()
                end = coord
                new_exons += new_exon_definition(start, end, read_rna, opened_exon[-1], jj, jj_gene)
                pp = opened_exon.pop()
                opened -= 1
            elif opened == 0:
                if first_3prime is None:
                    new_exons += __half_exon('5prime', jj, read_rna)
                else:
                    new_exons += new_exon_definition(first_3prime.get_ss_3p(),
                                                     coord, read_rna, first_3prime,
                                                     jj, jj_gene)
            last_5prime = jj
            #end elif opened
        else:
            if opened > 0:
                if not last_5prime is None:
                    end = last_5prime.get_ss_5p()
                    for ss in opened_exon:
                        if ss.get_gene != last_5prime.get_gene:
                            continue
                        start = ss.get_ss_3p()
                        new_exons += new_exon_definition(start, end, read_rna, ss, last_5prime, ss.get_gene)
                    last_5prime = None
                    opened = 0
                    opened_exon = []
                    first_3prime = jj
            else:
#                print "CHECK 2.2",coord,jj.get_ss_3p()
                last_5prime = None
                first_3prime = jj
            #end else ...
            opened_exon.append(jj)
            opened += 1

    for ss in opened_exon:
        new_exons += __half_exon('3prime', ss, read_rna)

    for (coord, jtype, jj) in junction_list:
        if jj.coverage.sum() < mglobals.MINREADS and not jj.is_annotated():
            junction_list.remove((coord, jtype, jj))
            del jj
            continue
        if jj.donor is None and jj.acceptor is None:
#            print "JUNCTIONS MISSING EXONS",jj.donor, jj.acceptor, jj.readN.sum(), jj.start, jj.end
            junction_list.remove((coord, jtype, jj))
            del jj

    print "FOUND new %d exons" % new_exons
    return 

