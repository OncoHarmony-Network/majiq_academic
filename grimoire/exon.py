#!/usr/bin/python

import os

import numpy as np

from grimoire.junction import Junction
from src import config



# FLAGS
ANNOTATED = 0b0001
INTRON = 0b0010
MISS_START = 0b0100
MISS_END = 0b1000


class Exon:
    __eq__ = lambda self, other: self.start < other.end and self.end > other.start
    __ne__ = lambda self, other: self.start >= other.end or self.end <= other.start
    __lt__ = lambda self, other: self.end <= other.start
    __le__ = lambda self, other: self.end < other.start or (self.start < other.end and self.end > other.start)
    __gt__ = lambda self, other: self.start >= other.end
    __ge__ = lambda self, other: self.start >= other.end or (self.start < other.end and self.end > other.start)

    def __init__(self, start, end, gene, strand, annot=False, isintron=False):

        self.flag = 0b0000
        if start == EMPTY_COORD:
            start = end - 10
            self.flag |= 0b0100
        if end == EMPTY_COORD:
            end = start + 10
            self.flag |= 0b1000

        self.start = start
        self.end = end
        self.gene_name = gene.get_id()
        self.exonTx_list = []
        self.exonRead_list = []
        self.ss_3p_list = []
        self.ss_5p_list = []
        self.id = "%s:%d-%d" % (self.gene_name, start, end)
        # self.strand = strand
        self.gc_content = 0
        self.coverage = np.zeros(shape=config.num_experiments)
        self.score = None
        self.pcr_name = None
        self.pcr_candidate = None
        self.db_coord = (start, end)

        self.flag |= 0b0001 if annot else 0b0000
        self.flag |= 0b0010 if isintron else 0b0000

    def __hash__(self):
        return hash(self.id) ^ hash(self.gene_name)

    def get_id(self):
        return self.id

    # def get_strand(self):
    # return self.strand

    def get_coordinates(self):
        """
         .. function:: get_coordinates(self)
            Get the exon start and end.
            :rtype: tuple of 2 integers, (exon start, exon end)
        """
        return self.start, self.end

    def get_gene(self):
        return config.gene_tlb[self.gene_name]

    def get_strand(self):
        return self.get_gene().get_strand()

    def is_annotated(self):
        return self.flag & ANNOTATED == ANNOTATED

    def is_intron(self):
        return self.flag & INTRON == INTRON

    def is_miss_start(self):
        return self.flag & MISS_START == MISS_START

    def is_miss_end(self):
        return self.flag & MISS_END == MISS_END

    def get_ir(self):

        return self.is_intron()

    def get_annotated_exon(self):
        return self.exonTx_list

    def get_rna_ss(self):
        ss3 = set()
        ss5 = set()

        for ii in self.exonRead_list:
            ss3.add(ii.start)
            ss5.add(ii.end)
        ss3_l = sorted(list(ss3))
        ss5_l = sorted(list(ss5))
        return ss3_l, ss5_l

    def get_length(self):
        if self.end is None or self.start is None:
            ln = 0
        else:
            ln = self.end - self.start
        return ln

    def set_ir(self, ir):
        self.flag |= 0b0011 if ir else 0b0000

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
        self.start = min(self.start, start)
        self.end = max(self.end, end)

        for ii in self.exonRead_list:
            if start == ii.start and end == ii.end:
                res = ii
                res.add_5prime_junc(s5p_junc)
                res.add_3prime_junc(s3p_junc)

                #self.exonRead_list.append(res)
                break
        else:
            for idx1, i1 in enumerate(self.ss_3p_list):
                if i1 != start:
                    continue
                if self.ss_5p_list[idx1] == end:
                    break
            else:
                self.ss_3p_list.append(start)
                self.ss_5p_list.append(end)
            res = ExonRead(start, end, s3p_junc, s5p_junc, read_seq)
            self.exonRead_list.append(res)
        return res

    def add_new_definition(self, start, end, trncpt):
        self.start = min(self.start, start)
        self.end = max(self.end, end)
        for ii in self.exonTx_list:
            if start == ii.start and end == ii.end:
                res = ii
                break
        else:
            #            if self.strand == '+':
            self.ss_3p_list.append(start)
            self.ss_5p_list.append(end)
            #            else:
            #                self.ss_5p_list.append(start)
            #                self.ss_3p_list.append(end)
            res = ExonTx(start, end, trncpt)
            self.exonTx_list.append(res)
        return res

    def get_coverage(self, exp_idx):
        return self.coverage[exp_idx]

    def get_gc_content(self):
        return self.gc_content

    def get_total_read_num(self, exp_idx):
        ex_reads = self.coverage[exp_idx].sum()
        junc3 = self.get_junctions('3prime')
        for j3 in junc3:
            ex_reads += j3.get_read_num(exp_idx)

        junc5 = self.get_junctions('5prime')
        for j5 in junc5:
            ex_reads += j5.get_read_num(exp_idx)
        return ex_reads

    def update_coverage(self, exp_idx, num):
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
                        if junc is None:
                            continue
                        jlist.add(junc)
                else:
                    for junc in ex.p5_junc:
                        if junc is None:
                            continue
                        jlist.add(junc)
        return jlist

    # def print_triplet_coord(self, fp):
    #     gene = mglobals.gene_tlb[self.gene_name]
    #     chrom = gene.chromosome
    #     strand = gene.get_strand()
    #     start_a = self.start
    #     end_a = self.end
    #
    #     if strand == '-':
    #         id_a = gene.exonNum - (self.id - 1)
    #         vid_c1 = self.id
    #         vid_c2 = self.id - 2
    #     else:
    #         id_a = self.id
    #         vid_c1 = self.id - 2
    #         vid_c2 = self.id
    #
    #     start_c1, end_c1 = gene.exons[vid_c1].get_coordinates()
    #     start_c2, end_c2 = gene.exons[vid_c2].get_coordinates()
    #
    #     name = "%s.%s" % (gene.id, id_a)
    #
    #     fp.write("%s\t%d\t%d\t%s_C1\t0\t%s\n" % (chrom, start_c1, end_c1, name, strand))
    #     fp.write("%s\t%d\t%d\t%s_A\t0\t%s\n" % (chrom, start_a, end_a, name, strand))
    #     fp.write("%s\t%d\t%d\t%s_C2\t0\t%s\n" % (chrom, start_c2, end_c2, name, strand))

    def bed_format(self):
        bed_str = ""
        for eRead in self.exonRead_list:
            bed_str += "%s\n" % eRead.bed_format()

        return bed_str

    def ss_variant_counts(self, minreads=5):

        temp_set = set()
        for ss3p in self.ss_3p_list:
            for exread in self.exonRead_list:
                if ss3p != exread.start:
                    continue
                sum_reads = 0
                for jj in exread.p3_junc:
                    if jj is None:
                        continue
                    sum_reads += jj.coverage.sum()
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
                    sum_reads += jj.coverage.sum()
                if sum_reads >= minreads:
                    temp_set.add(ss5p)
        local_5p = len(temp_set)

        if local_3p > 19:
            local_3p = 19
        if local_5p > 19:
            local_5p = 19

        return local_3p, local_5p


class ExonRead(object):
    def __init__(self, start, end, pre_junc, post_junc, rna_seq=None):
        self.start = start
        self.end = end
        self.RNASeq = rna_seq
        self.p3_junc = []
        if not pre_junc is None:
            self.p3_junc.append(pre_junc)
        self.p5_junc = []
        if not post_junc is None:
            self.p5_junc.append(post_junc)

    def get_coordinates(self):
        return self.start, self.end

    def get_5p_junc(self):
        return self.p5_junc

    def add_5prime_junc(self, junc):

        if not junc is None and not junc in self.p5_junc:
            self.p5_junc.append(junc)

    def add_3prime_junc(self, junc):
        if not junc is None and not junc in self.p3_junc:
            self.p3_junc.append(junc)

    def bed_format(self):
        chrom = self.exon.get_gene().get_chromosome()
        strng = "%s\t%s\t%s\t.\t0\t.\t%s\t%s" % (chrom, self.start, self.end, self.start, self.end)
        return strng


class ExonTx(object):
    __eq__ = lambda self, other: self.start == other.start and self.end == other.end
    __ne__ = lambda self, other: self.start != other.start or self.end != other.end
    __lt__ = lambda self, other: self.start < other.start or (self.start == other.start and self.end < other.end)
    __le__ = lambda self, other: self.start <= other.start or (self.start == other.start and self.end <= other.end)
    __gt__ = lambda self, other: self.start > other.start or (self.start == other.start and self.end > other.end)
    __ge__ = lambda self, other: self.start >= other.start or (self.start == other.start and self.end >= other.end)

    def __init__(self, start, end, trnscpt, intron=False):
        self.start = start
        self.end = end
        self.transcript_name = [trnscpt.get_id()]
        self.gene_name = trnscpt.get_gene().get_id()
        self.p3_junc = []
        self.p5_junc = []
        self.ir = intron

    def get_coordinates(self):
        return self.start, self.end

    def add_transcript(self, trans):
        self.transcript_name.append(trans.get_id())

    def add_5prime_junc(self, p5_junc):
        if p5_junc not in self.p5_junc:
            self.p5_junc.append(p5_junc)

    def add_3prime_junc(self, p3_junc):
        if p3_junc not in self.p3_junc:
            self.p3_junc.append(p3_junc)

    def get_transcript(self):
        res = []
        for tx_name in self.transcript_name:
            res.append(config.gene_tlb[self.gene_name].get_transcript(tx_name))
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
        if self.end - intron_coords[1] + 1 > 5:
            txex1 = gn.new_annotated_exon(intron_coords[1] + 1, self.end, self.get_transcript()[0], bl=False)
            txex1.p5_junc.extend(self.p5_junc)
            res.append(txex1)
            txex1.junction_consistency()
            exb1 = True
        if intron_coords[0] - 1 - self.start > 5:
            txex2 = gn.new_annotated_exon(self.start, intron_coords[0] - 1, self.get_transcript()[0], bl=False)
            txex2.p3_junc.extend(self.p3_junc)
            exb2 = True
            res.append(txex2)
            txex2.junction_consistency()

        exb = exb1 & exb2
        if exb:
            # #creating intron exon
            #         intron = gn.new_annotated_exon(intron_coords[0], intron_coords[1], self.get_transcript()[0],
            #                                        bl=False, intron=True)
            #         res.append(intron)
            #         junc1 = gn.exist_junction(txex2.end, intron_coords[0])
            #         if junc1 is None:
            #             junc1 = Junction(txex2.end, intron_coords[0], None, None, gn, annotated=True)
            # #           junc1.add_donor(txex2)
            # #           junc1.add_acceptor(intron)`
            #         txex2.p5_junc.append(junc1)
            #         intron.p3_junc.append(junc1)
            #
            #         junc2 = gn.exist_junction(intron_coords, txex1.start)
            #         if junc2 is None:
            #             junc2 = Junction(intron_coords, txex1.start, None, None, gn, annotated=True)
            # #           junc.add_donor(intron)
            # #           junc.add_acceptor(txex1)
            #         intron.p5_junc.append(junc2)
            #         txex1.p3_junc.append(junc2)

            junc = gn.exist_junction(txex2.end, txex1.start)
            if junc is None:
                junc = Junction(txex2.end, txex1.start, None, None, gn, annotated=True)
            #                junc.add_donor(txex2)
            #                junc.add_acceptor(txex1)
            txex2.p5_junc.append(junc)
            txex1.p3_junc.append(junc)

        del self
        return res

    def junction_consistency(self):

        j5_list = []
        for j5 in self.p5_junc:
            jcoord = j5.get_coordinates()
            if self.start <= jcoord[0] <= self.end:
                # j5.add_donor(self)
                j5_list.append(j5)

        j3_list = []
        for j3 in self.p3_junc:
            jcoord = j3.get_coordinates()
            if self.start <= jcoord[1] <= self.end:
                # j3.add_acceptor(self)
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
                            introns.append((last_p5 + 1, max(p3 - 1, last_p5 + 1)))
                            in_found = False
                        jdx += 1
                    else:
                        last_p5 = p5
                        in_found = True
                        break

            for idx, txex in enumerate(list_exontx):
                for intr in introns:
                    if not txex.overlaps(intr[0], intr[1]):
                        # if intr[0] > txex.end or (intr[0] <= txex.end < intr[1]):
                        # txex.ir = True
                        continue
                    ''' intron retention'''
                    gne.add_ir_definition(intr[0], intr[1])
                    # LSV_IR(txex.start, txex.end, [], gne)
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
                # txex.exon = ex
                for p3_junc in txex.p3_junc:
                    p3_junc.add_acceptor(ex)
                for p5_junc in txex.p5_junc:
                    p5_junc.add_donor(ex)
            exlist.append(ex)
        return exlist


def print_list_exons(list_ex, msg=""):
    # list_ex.sort()
    print "%%%%%%%%%%%%LIST_EXONS %s" % msg
    for ex in list_ex:
        print "\t\t", ex.get_coordinates(), ex
    print "%%%%%%%%%%%%%%%%%%%%%%"


num_it = 0


def collapse_list_exons(listexons, gne):
    global num_it
    num_it += 1
    # print "[%s] INIT COLLAPSE_LIST EXONS "%(num_it)
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
        if idx == len(listexons) - 1:
            exlist.extend(ex.collapse(overlp, gne))
    num_it -= 1
    return exlist


def __half_exon(ss_type, junc, read_rna):
    gene = junc.get_gene()
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
                # print "half",type,"::",ex_start, ex_end, junc.start, junc.end, end
            # if end - start < 10 : continue
            res = ex.add_new_read(start, end, read_rna, to, frm)
            if res:
                ex.ss_3p_list.append(start)
                ex.ss_5p_list.append(end)

            break
    return 0


EMPTY_COORD = -1


def new_exon_definition(start, end, read_rna, s3prime_junc, s5prime_junc, gene, isintron=False):
    if end - start < 5:
        return 0

    ex = gene.exist_exon(start, end)
    new_exons = 0
    half = False

    if ex is None:
        new_exons = 1
        in_db = False
        for xx in gene.get_ir_definition():
            if start <= xx[1] and end >= xx[0]:
                in_db = True
                break
        ex = Exon(start, end, gene, gene.get_strand(), annot=in_db, isintron=isintron)
        gene.add_exon(ex)
    else:
        coords = ex.get_coordinates()
        if start != EMPTY_COORD and start < (coords[0] - config.get_max_denovo_difference()):
            if gene.exist_exon(start, start + 10) is None:
                new_exons += 1
                ex1 = Exon(start, EMPTY_COORD, gene, gene.get_strand(), isintron)
                cc = ex1.get_coordinates()
                s3prime_junc.add_acceptor(ex1)
                gene.add_exon(ex1)
                ex1.add_new_read(cc[0], cc[1], read_rna, s3prime_junc, None)
            half = True

        if end != EMPTY_COORD and end > (coords[1] + config.get_max_denovo_difference()):
            if gene.exist_exon(end - 10, end) is None:
                new_exons += 1
                ex2 = Exon(EMPTY_COORD, end, gene, gene.get_strand(), isintron)
                cc = ex2.get_coordinates()
                s5prime_junc.add_donor(ex2)
                gene.add_exon(ex2)
                ex2.add_new_read(cc[0], cc[1], read_rna, None, s5prime_junc)
            half = True

    if not half:
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

    junction_list.extend(gene.get_all_ss())
    junction_list.sort()
    for (coord, jtype, jj) in junction_list:

        if not jj.is_reliable() and not jj.is_annotated():
            continue

        jj_gene = jj.get_gene()
        if jtype == '5prime':
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
            # end elif opened
        else:
            if opened > 0:
                if not last_5prime is None:
                    end = last_5prime.get_ss_5p()
                    for ss in opened_exon:
                        if ss.get_gene() != last_5prime.get_gene():
                            continue
                        start = ss.get_ss_3p()
                        new_exons += new_exon_definition(start, end, read_rna, ss, last_5prime, ss.get_gene())
                    last_5prime = None
                    opened = 0
                    opened_exon = []
                    first_3prime = jj
            else:
                last_5prime = None
                first_3prime = jj
            # end else ...
            opened_exon.append(jj)
            opened += 1

    for ss in opened_exon:
        new_exons += __half_exon('3prime', ss, read_rna)

    for (coord, jtype, jj) in junction_list:
        if not jj.is_reliable() and not jj.is_annotated():
            junction_list.remove((coord, jtype, jj))
            del jj
            continue
        if jj.get_donor() is None and jj.get_acceptor() is None:
            junction_list.remove((coord, jtype, jj))
            del jj

    # print "FOUND new %d exons" % new_exons
    return


def set_exons_gc_content(chrom, exon_list):
    fastadir_path = "%s/" % config.genome_path

    # print "Loading chromosome... %s"%chrom
    chrom_path = fastadir_path + chrom + ".fa"
    if not os.path.exists(chrom_path):
        return
    chrom_file = open(chrom_path)
    loaded_chrom = []
    for chrom_line in chrom_file:
        if not chrom_line.startswith(">"):
            loaded_chrom.append(chrom_line.strip("\n"))
    loaded_chrom = ''.join(loaded_chrom)
    chrom_file.close()
    # print exon_list
    for exon in exon_list:
        strt, end = exon.get_coordinates()
        if end - strt < 5:
            continue
        sequence = loaded_chrom[strt:end]
        #reverse the sequence if the strand is reverse
        sequence = sequence.lower()
        if exon.get_strand() == "-":
            new_seq = []
            for char in sequence[::-1]:
                if char == 'g':
                    new_seq.append('c')
                elif char == 'c':
                    new_seq.append('g')
                elif char == 'a':
                    new_seq.append('t')
                elif char == 't':
                    new_seq.append('a')
                else:
                    new_seq.append(char)
            sequence = ''.join(new_seq)
        exon.set_gc_content(sequence)
    del loaded_chrom
