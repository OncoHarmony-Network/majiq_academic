from gene import Gene, Transcript
from exon import Exon, detect_exons
from junction import  Junction
from rna_reads import RNA_read, parse_read
from RNAexperiment import RNAexperiment
import grimoire.utils.utils as utils
import pysam
import gc
import numpy as np
import mglobals
from sys import getrefcount
import pdb



def __cross_junctions(read):
    '''
      This part will parse the jI from STAR
      # Column 17: jI:B:I,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1- based)
      # jI:B:I,-1 means that the read doesn't cross any junction
    '''
    jlist = []
    cross = False
    try: 
        list_junc = read.opt('jI')
        if len(list_junc) > 1 and list_junc[0] != -1:
            for idx in range(0,len(list_junc),2):
                junc_start = int(list_junc[idx])-1
                junc_end = int(list_junc[idx+1])+1
                jlist.append((junc_start,junc_end))
            #end for idx ...
            cross = True
        #end else ...
    except KeyError:
#    if len(jlist) != 0: print "STAR ::",jlist
#    print"THIS IS NOT a WELL defined STAR output"
        off = 0
        for op, num in read.cigar:
            if op in [0,5,6,7,8] : 
                off += num
            elif op in [1,5]:
                off += 0
            elif op == 2:
                off+= num
            elif op == 3:
                jlist.append((read.pos+off,read.pos+off+num+1))
                off += num
#    if len(jlist) !=0 : print "NOSTAR:", jlist, read.cigar
        
    return (cross,jlist)


def __is_unique(read):
    unique = True
    try:
        if read.opt('NH') >1 :
            unique = False
    except KeyError:
        if read.flag & 0x100 == 1: 
            unique = False
    return True
    #return unique

def __get_num_reads(read):
    try:
        nreads = int(read.opt('HI'))
    except KeyError:
        nreads = 1
#    return 1
    return nreads

def count_mapped_reads( filename, exp_idx):
    stats = pysam.flagstat(filename)
    mapped_reads = int(stats[2].split()[0])
    mglobals.num_mapped_reads[exp_idx] = mapped_reads


def is_neg_strand (read):
    res = False
    if read.flag & 0x10 == 0x10 : 
#        print "FLAG",read.flag
        res = True

    return res 


import pdb
def read_sam_or_bam(filenames, gene_list, readlen, chrom):
#def read_sam_or_bam(filename, gene_list, readlen, chrom, exp_index):

    counter = [0] * 6
    samfile = [ pysam.Samfile( xx, "rb" ) for xx in filenames ]
    temp_ex = []
    NUnum = 0
    skip_gene = 0
    non_skip = 0

    for strand in ('+','-') :
        for gne in gene_list[strand]:
            junctions = []
            strt,end = gne.get_coordinates()
            j_list = gne.get_all_junctions()
            ex_list = gne.get_exon_list()

            for exp_index in range(len(filenames)):
                try:
                    read_iter = samfile[exp_index].fetch( chrom, strt, end)
                except ValueError:
                    continue
                for read in read_iter:
                    strand_read = '+' if not is_neg_strand(read) else '-'
    #                strand_read = '+' if not read.is_reverse else '-'
    #                print "STRAND",strand_read, strand
                    if strand_read != strand: continue
                    unique = __is_unique(read)
                    if not unique : 
                        NUnum  += 1
                        continue
                    nreads = __get_num_reads(read)
                    gne.add_read_count(nreads, exp_index)
                    is_cross, junc_list = __cross_junctions(read)
                    r_start = read.pos
                    if r_start < strt or r_start > end : continue
    #                print "READ",strt,end, r_start, read
    
                    for ex_idx in range(len(ex_list)):
                        ex_start, ex_end = ex_list[ex_idx].get_coordinates()
                        if r_start >= ex_start and r_start <= ex_end :
                            ex_list[ex_idx].update_coverage(exp_index, nreads)
                            temp_ex.append(ex_list[ex_idx])
                            break
                    #else:
                    #    break
    
                    if not is_cross : continue
                    nc = read.seq.count('C') + read.seq.count('c') 
                    ng = read.seq.count('g') + read.seq.count('G') 
                    gc_content = float( nc + ng ) / float(len(read.seq))
                    for (junc_start,junc_end) in junc_list:
                        if junc_start - r_start > readlen : 
                            r_start = junc_start - (readlen - 16) -1
                        elif junc_start - r_start >= readlen - 8 or junc_start -r_start <= 8: continue
                        found = False
                        if junc_end - junc_start < 10 :
                            counter[0] +=1
                            continue

                        for jj in j_list:
                            ( j_st,j_ed ) = jj.get_coordinates()
                            if   j_st > junc_start  or (j_st == junc_start and j_ed > junc_end): break
                            elif j_st < junc_start  or (j_st == junc_start and j_ed < junc_end): continue
                            elif junc_start == j_st and junc_end == j_ed:
                                ''' update junction and add to list'''
                                counter[3] +=1
                                jj.update_junction_read(exp_index,nreads,r_start, gc_content, unique)
                                if not (junc_start,'5prime',jj) in junctions:
                                    junctions.append((junc_start,'5prime',jj))
                                    junctions.append((junc_end,'3prime',jj))
                                found = True
                                break
                            #end elif junc_start == ...
                        #end for jj in j_list

                        if not found :
                            ''' update junction and add to list'''
                            junc = None
                            for (coord,t,junc) in junctions:
                                if (junc.start == junc_start and junc.end == junc_end):
                                    junc.update_junction_read(exp_index, nreads, r_start, gc_content, unique)
                                    if not (junc_start,'5prime',junc) in junctions:
                                        junctions.append((junc_start,'5prime',junc))
                                        junctions.append((junc_end,'3prime',junc))
                                    break
                                #end if (j.start) == ...
                            #end for (coord,t,j) ...
                            if junc is None:
                                '''mark a new junction '''
                                counter[4] += 1
                                junc = Junction( junc_start, junc_end, None, None, gne, readN=nreads)
                                junc.update_junction_read(exp_index, nreads, r_start, gc_content, unique)
                                junctions.append((junc_start,'5prime',junc))
                                junctions.append((junc_end,'3prime',junc))
                        #end if not found ...
                    #end for junc ...
    #            print "JJJunctions", junctions
            if len(junctions) > 0 :
                detect_exons(gne, junctions, None)
            gne.prepare_exons()
            re = gne.calculate_RPKM(exp_index, mglobals.num_mapped_reads[exp_index])
            if re == 0 :
                #print "SIN READS??? :",gene.id
                skip_gene +=1
            else:
                non_skip +=1
    print "INVALID JUNC",    counter[0]
    print "READ WRONG GENE", counter[1]
    print "READ IN GENE",    counter[2]
    print "READ FOUND JUNC", counter[3]
    print "READ NEW JUNC",   counter[4]
    print "READ ALL JUNC",   counter[5]
    print "Non Unique", NUnum
            
    print "Skipped genes without exons",skip_gene
    print " Non skipped",non_skip
    return 



def read_transcript_ucsc( filename ):
    print "Reading transcript ucsc format"
    BUFFER = int(10E6) #10 megabyte buffer
    file = open(filename, 'r')
    text = file.readlines(BUFFER)
    all_genes = {}
    while text != []:
        print "New Buffer"
        for t in text:
            t = t.strip()
            if t.startswith('#') or t == '':
                continue
            tab = t.split('\t')
            chrom = tab[0]
            ''' BED format file is 0 based, so the first base of the chromosome is 0. That means that the first real base of an exon is start+1'''
            start = int(tab[1])+1
            end = int(tab[2])
            strand = tab[5]
            gene_id = tab[3]
            transcript_id = tab[3]
            off_start = tab[11].split(',')
            off_len = tab[10].split(',')
            nblocks = int(tab[9])

            ex_start = [0]*nblocks
            ex_end = [0]*nblocks
            for ii in range(nblocks):
                ex_start[ii] = int(off_start[ii])+start
                ex_end[ii] = ex_start[ii] + int(off_len[ii])-1

            if not chrom in all_genes:
                all_genes[chrom] = {'+':[],'-':[]}
            gn = Gene(gene_id,chrom,strand,start,end)
            gene = gn.is_gene_in_list(all_genes[chrom][strand], gene_id)
            if not gene is None:
                del gn
                gn = gene
            else:
                all_genes[chrom][strand].append(gn)

            trcpt = Transcript(transcript_id, gn,start,end )
            gn.add_transcript(trcpt)
            pre_end = None
            pre_ex = None
            pre_txex = None
            for ii in xrange(nblocks):
                start = ex_start[ii]
                end = ex_end[ii]

                txex = gn.new_annotated_exon(start, end, trcpt)
                trcpt.add_exon(txex)
                junc = gn.exist_junction(pre_end,start)
                if junc is None:
                    junc = Junction(pre_end,start, None, None,gn,annotated=True)
                trcpt.add_junction(junc)
                txex.add_3prime_junc(junc)
                if not pre_end is None:
                    pre_txex.add_5prime_junc(junc)
                pre_end = end
                pre_txex = txex
            #END for ii in range(nblocks)
            if not pre_end is None: 
                junc = Junction( pre_end,None, None, None, gn,annotated=True )
                trcpt.add_junction(junc)
                pre_txex.add_5prime_junc(junc)
            trcpt._sort_in_list(strand)
        #end for t in text
        gc.collect()
        text = file.readlines(BUFFER)
    #END WHILE
    file.close()

    n_genes = 0
    for chr in all_genes.keys():
        temp_ex = []
        for strand, gg in all_genes[chr].items():
            n_genes += len(gg)
            all_genes[chr][strand] = sorted(gg)
            for gene in all_genes[chr][strand]:
#                print "COLLAPSING GENE", gene.get_id()
                gene.collapse_exons()
                temp_ex.extend(gene.get_exon_list())
#                print "NUM IR:", len(gene.ir_list)
        print "Calculating gc_content.........",
        utils.set_exons_gc_content(chr, temp_ex)
        print "Done."

    print "NUM_GENES",n_genes
    return all_genes


def read_bed_pcr( filename , list_genes):
    
    input_f = open(filename,'r')
    readlines = input_f.readlines()
    alt_exon = []
    pre_chrom = ''
    gene_list = {}
    lnum = 0
    while lnum < len(readlines):
        more = True
        event = {}
        while more :
            rl = readlines[lnum]
            tab = rl.strip().split()
            t = tab[3].split('_')
            event['name'] = '_'.join(t[:-1])
            reg = t[-1]
            if reg == 'C2' : more = False
            event['chrom'] = tab[0]
            event['strand'] = tab[5]
            event[reg] = [int(tab[1]), int(tab[2])]
            score = float(tab[4].split('|')[0])
            if reg in ['A','A2']: alt_exon = event[reg]
            lnum +=1

        chrom = event['chrom']
        strand = event['strand']
        print event
        exon_start = event['C1'][0]
        exon_end = event['C1'][1]

        if chrom != pre_chrom :
            try:
                gene_list = list_genes[chrom]
            except KeyError:
                continue

            pre_chrom = chrom
            idx = {'+':0, '-':0}
        name = event['name']
        while idx[strand]< len(gene_list[strand]):
            gn = gene_list[strand][idx[strand]]
            (g_start, g_end) = gn.get_coordinates()
            if exon_end < g_start:  break
            elif exon_start > g_end :
                idx[strand] +=1
                continue
            ex = gn.exist_exon(exon_start,exon_end)
            if ex is None : break
            ex.set_pcr_score(name, score, alt_exon)

            break

