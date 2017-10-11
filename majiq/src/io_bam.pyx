from majiq.src.constants import *
import numpy as np
from majiq.grimoire.exon cimport Intron
from majiq.grimoire.junction cimport Junction
from majiq.grimoire.junction import Junction
from majiq.src.io import dump_junctions
from majiq.src.multiproc import QueueMessage

from majiq.src.config import Config
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
# from libcpp cimport bool
import pysam
import cython

# READING BAM FILES

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef __cross_junctions(AlignedSegment read):
    """
      This part will parse the jI from STAR
      # Column 17: jI:B:I,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1- based)
      # jI:B:I,-1 means that the read doesn't cross any junction
    """

    cdef list jlist = []
    cdef bint cross = False
    cdef int idx
    cdef unsigned int op
    cdef unsigned int num
    cdef unsigned int off
    cdef int junc_end, junc_start

    # if len(jlist) != 0: print "STAR ::",jlist
    # print"THIS IS NOT a WELL defined STAR output"
    off = 0
    for op, num in read.cigar:
        if op in [0, 5, 6, 7, 8]:
            off += num
        elif op in [1, 5]:
            off += 0
        elif op == 2:
            off += num
        elif op == 3:
            jlist.append((read.pos + off, read.pos + off + num + 1))
            off += num
            cross = True
            # if len(jlist) !=0 : print "NOSTAR:", jlist, read.cigar

    return cross, jlist, off


cpdef int find_introns(str filename, dict list_introns, float intron_threshold, queue, str gname) except -1:

    cdef AlignmentFile samfl
    cdef int nchunks, chunk_len
    cdef bint b_included
    cdef list junc_list
    cdef int intron_len, ii, i_st, i_nd, ub, lb, num_bins = 10
    cdef float val
    cdef str chrom, strand

    samfl = open_rnaseq(filename)


    for gne_id, chrom, strand, i_st, i_nd in list_introns.keys():
        #print(chrom, i_st, i_nd)
        intron_len = (i_nd - i_st)
        nchunks = 1 if intron_len <= MIN_INTRON_LEN else num_bins
        chunk_len = int(intron_len / nchunks)
        b_included = True
        for ii in range(nchunks):
            lb = i_st + ii*chunk_len
            ub = i_st + (ii+1)*chunk_len -1
            ub = min(ub, i_nd)
            try:
                val = samfl.count(contig=chrom, start=lb, stop=ub, until_eof=True, read_callback=__is_unique,
                                reference=None, end=None)
            except ValueError as e:
                #
                # logging.error('\t[%s]There are no reads in %s:%d-%d' % (info_msg, gne['chromosome'], gne['start'], gne['end']))
                continue

            val /= (ub-lb)
            b_included = b_included and (val>=intron_threshold)

        if b_included :
            qm = QueueMessage(QUEUE_MESSAGE_BUILD_INTRON, (gne_id, i_st, i_nd, gname), 0)
            queue.put(qm, block=True)
            #list_introns[(gne_id, chrom, strand, i_st, i_nd)] += 1

    samfl.close()
    return 0


cdef inline dict __read_STAR_junc_file(str filename, set in_jj, bint stranded):

    cdef dict out_dd = {}
    cdef dict strand_list
    cdef object fp
    cdef str ln
    cdef list tab
    cdef str stnd
    cdef tuple jid

    if stranded:
        strand_list = {'1':['+'], '2':['-'], '0':['+', '-'] }
    else:
        strand_list = {'1':['.'], '2':['.'], '0':['.'] }

    with open(filename, 'r') as fp:
        for ln in fp.readlines():
            tab = ln.strip().split()

            for stnd in strand_list[tab[3]]:
                jid = (tab[0], stnd, tab[1], tab[2])
                if jid not in in_jj:
                    try:
                        out_dd[tab[0]].add((tab[1], tab[2], stnd))
                    except KeyError:
                        out_dd[tab[0]] = set()
                        out_dd[tab[0]].add((tab[1], tab[2], stnd))
    return out_dd



cdef dict __read_juncs_from_bam(str filename, set in_jj, bint stranded):

    cdef AlignedSegment read
    cdef AlignmentFile samfl
    cdef int readlen
    cdef bint is_cross
    cdef dict out_dd = {}
    cdef list junc_list
    cdef int junc_start, junc_end, rend
    cdef tuple jid, jid2
    cdef str strand

    samfl = open_rnaseq(filename)
    read_iter = samfl.fetch(until_eof=True)


    for read in read_iter:
        is_cross, junc_list, rend = __cross_junctions(read)
        if not __is_unique(read) or not is_cross:
            continue
        chrom = samfl.getrname(read.reference_id)
        readlen = len(read.seq)
        for junc_start, junc_end in junc_list:

            if (junc_start - read.pos >= readlen - MIN_BP_OVERLAP) or (junc_start - read.pos <= MIN_BP_OVERLAP) or \
                (junc_end - junc_start < MIN_JUNC_LENGTH):
                    continue
            if stranded:
                strand = '-' if read.is_reverse() else '+'
                jid = (chrom, strand, junc_start, junc_end)
                jid2 = (junc_start, junc_end, strand)
            else:
                jid = (chrom, '.', junc_start, junc_end)
                jid2 = (junc_start, junc_end, '.')

            if chrom not in out_dd:
                out_dd[chrom] = set()

            if jid not in in_jj and jid2 not in out_dd[chrom]:
                out_dd[chrom].add(jid2)


    return out_dd


def read_juncs(str fname, bint is_junc_file, dict dict_exons, dict dict_genes, dict junctions, bint stranded, queue,
               str gname):

    cdef str ln, chrom, gid
    cdef list tab
    cdef Junction jj
    cdef tuple junc
    cdef int gidx=0, ngenes=0
    cdef dict new_junctions
    cdef set jj_set
    cdef set set_junctions

    if stranded:
        set_junctions = set([(dict_genes[gene_id]['chromosome'], dict_genes[gene_id]['strand'], yy.start, yy.end)
                              for gene_id, xx in junctions.items() for yy in xx.values()])
    else:
        set_junctions = set([(dict_genes[gene_id]['chromosome'], '.', yy.start, yy.end)
                             for gene_id, xx in junctions.items() for yy in xx.values()])

    if is_junc_file:
       new_junctions = __read_STAR_junc_file(fname, set_junctions, stranded)
    else:
       new_junctions = __read_juncs_from_bam(fname, set_junctions, stranded)

    for chrom, jj_set in new_junctions.items():
        gne_list = sorted([xx for xx in dict_genes.values() if xx['chromosome']==chrom], key=lambda x: (x['start'], x['end']))
        ngenes = len(gne_list)

        if stranded:
            init_gidx = {'+': 0, '-': 0}
        else:
            init_gidx = {'.': 0}

        found = False
        for junc in sorted(jj_set):
            gidx = init_gidx[junc[2]]
            possible_genes = []
            while gidx < ngenes:
                gobj = gne_list[gidx]
                if stranded and gobj['strand'] != junc[2]:
                    gidx += 1
                    continue

                if gobj['end'] < junc[0]:
                    gidx += 1
                    init_gidx[junc[2]] += 1
                    continue

                if gobj['start'] > junc[1]:
                    if not found and len(possible_genes) > 0:
                        for gobj in possible_genes:
                            junc_obj = Junction(junc[0],  junc[1], gid, -1, annot=False)
                            qm = QueueMessage(QUEUE_MESSAGE_BUILD_JUNCTION, (gid, junc[0], junc[1], gname), 0)
                            queue.put(qm, block=True)
                            junctions[gobj['id']][junc[:-1]] = junc_obj

                    break
                if gobj['start']<=junc[1] and gobj['end']>= junc[0]:
                    gid = gobj['id']
                    start_sp = [jj.start for ex in dict_exons[gid] for jj in ex.ob if jj.start > 0 and jj.end > 0]
                    end_sp = [jj.end for ex in dict_exons[gid] for jj in ex.ib if jj.start > 0 and jj.end > 0]
                    #if junc[0] == 20757728 : print [dict_exons[gid]
                    if junc[0] in start_sp or junc[1] in end_sp:
                        found = True
                        junc_obj = Junction(junc[0],  junc[1], gid, -1, annot=False)
                        qm = QueueMessage(QUEUE_MESSAGE_BUILD_JUNCTION, (gid, junc[0], junc[1], gname), 0)
                        queue.put(qm, block=True)
                        junctions[gobj['id']][junc[:-1]] = junc_obj
                    else:
                        possible_genes.append(gobj)
                    gidx +=1



cdef inline bint __is_unique(AlignedSegment read):
    return not(read.flag & 0x100 == 0x100)


cdef inline int __get_num_reads(AlignedSegment read):
    return 1


cdef inline bint _match_strand(AlignedSegment read, str gene_strand):
    majiq_config = Config()
    res = True
    if majiq_config.strand_specific:
        #TODO: REMOVE
        if (read.flag & 0x10 == 0x10 and gene_strand == b'+') or (read.flag & 0x10 == 0x00 and gene_strand == b'-'):
            res = True
        else:
            res = False
    return res


cpdef AlignmentFile open_rnaseq(str samfile):
    return pysam.Samfile(samfile, "rb")


cpdef long close_rnaseq(AlignmentFile samfl) except -1:
    samfl.close()

cpdef int read_sam_or_bam(object gne, AlignmentFile samfl, list matrx, dict junctions, list exon_list, list intron_list,
                          str info_msg='0', object logging=None) except -1:
    cdef int res
    res = _read_sam_or_bam(gne, samfl, matrx, junctions=junctions, exon_list=exon_list, intron_list=intron_list,
                           info_msg=info_msg, logging=logging)
    return res


# cpdef long rnaseq_intron_retention(dict gne, list list_exons, AlignmentFile samfl, list matrx, list out_junctions,
#                                    object logging=None) except -1:
#
#     cdef int res
#     res = _rnaseq_intron_retention(gne, list_exons, samfl, matrx, out_junctions, logging=None)
#     return res


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef int __junction_read(AlignedSegment read, list junc_list, int max_readlen, int effective_len,
                         list matrx, dict junctions) except -1:

    cdef float nc = read.seq.count('C') + read.seq.count('c')
    cdef float ng = read.seq.count('g') + read.seq.count('G')
    cdef float gc_content = float(nc + ng) / float(len(read.seq))
    cdef int readlen = len(read.seq)
    cdef int nreads = __get_num_reads(read)
    cdef long r_start = read.pos
    cdef Junction junc

    for (junc_start, junc_end) in junc_list:

        if junc_start - r_start > readlen:
            r_start_offset = junc_list[0][0] - r_start
            r_start = junc_start - r_start_offset

        if (junc_start - r_start >= readlen - MIN_BP_OVERLAP) or (junc_start - r_start <= MIN_BP_OVERLAP) or \
                (junc_end - junc_start < MIN_JUNC_LENGTH):
            continue

        left_ind = max_readlen - (junc_start - r_start) - MIN_BP_OVERLAP + 1
        if (junc_start, junc_end) in junctions:
            ''' update junction and add to list'''
            try:
                junc = junctions[(junc_start, junc_end)]
                junc.update_junction_read(nreads)
                if junc.index == -1:
                    junc.index = len(matrx)
                    matrx.append([0] * effective_len)

                matrx[junc.index][left_ind] += nreads
            except KeyError:

                continue


cdef int __intronic_read(AlignedSegment read, Intron intron, str gne_id, list junc_list, int max_readlen,
                         int effective_len, list matrx, dict junctions) except -1:

    cdef int nreads = __get_num_reads(read)
    cdef long r_start = read.pos
    cdef int indx, readlen = len(read.seq)
    cdef int rel_start, intron_idx, left_ind = max_readlen - (intron.start - r_start) - MIN_BP_OVERLAP
    cdef int offset = readlen - MIN_BP_OVERLAP

    if intron.start - r_start > readlen:
        r_start = intron.start - (readlen - MIN_BP_OVERLAP*2) - 1

    if r_start < intron.start - MIN_BP_OVERLAP:
        if intron.junc1 is None:
            try:
                intron.junc1 = junctions[(intron.start - 1, intron.start)]
            except KeyError:
                intron.junc1 = Junction(intron.start - 1, intron.start, gne_id, cov_idx=len(matrx), intron=True)
            intron.junc1_cov = [0] * effective_len
        intron.junc1.update_junction_read(nreads)

    elif (intron.end - offset ) < r_start < (intron.end+1):
        if intron.junc2 is None:
            try:
                intron.junc2 = junctions[(intron.end, intron.end+1)]
            except KeyError:
                intron.junc2 = Junction(intron.end, intron.end+1, gne_id, cov_idx=len(matrx), intron=True)
            intron.junc2_cov = [0] * effective_len
        intron.junc2.update_junction_read(nreads)

    else:
        # section 3
        if r_start > intron.start - 1:
            intron_idx = r_start - intron.start
            rel_start = int(intron_idx / intron.chunk_len)
            indx = -1 if rel_start >= intron.nchunks else rel_start
            intron.parts[indx] += nreads



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(True)  # turn on negative index wrapping for entire function

cdef int _read_sam_or_bam(object gne, AlignmentFile samfl, list matrx, dict junctions, list exon_list, list intron_list,
                          str info_msg='0', object logging=None) except -1:

    cdef unsigned int r_start, junc_start, junc_end, readlen, nc, ng
    cdef AlignedSegment read
    cdef bint unique, found
    cdef float gc_content
    cdef bint bb
    #cdef dict junctions = {} #{xx.get_coordinates(): xx for xx in gne.get_all_junctions(filter=False)}
    cdef int counter = 0
    cdef int tot_reads = 0
    cdef int effective_len = 0
    cdef list junc_list
    cdef Intron intron

    try:

        majiq_config = Config()
        effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1
        read_iter = samfl.fetch(gne['chromosome'], gne['start'], gne['end'], multiple_iterators=False)

        for read in read_iter:
            is_cross, junc_list, end_r = __cross_junctions(read)
            unique = __is_unique(read)
            if not _match_strand(read, gene_strand=gne['strand']) or read.pos < gne['start'] or not unique:
                continue

            tot_reads += 1
            if is_cross:
              __junction_read(read, junc_list, majiq_config.readLen, effective_len, matrx, junctions)
            for intron in intron_list:

                if read.pos >= intron.end:
                    continue
                elif end_r <= intron.start:
                    break
                __intronic_read(read, intron, gne['id'], junc_list, majiq_config.readLen, effective_len,
                                matrx, junctions)
        for intron in intron_list:
            if intron.junc1 is None or intron.junc2 is None:
                del intron
                continue

            if intron.junc1.nreads >= majiq_config.min_denovo and intron.junc2.nreads >= majiq_config.min_denovo:
                junctions[(intron.junc1.start, intron.junc1.end)] = intron.junc1
                junctions[(intron.junc2.start, intron.junc2.end)] = intron.junc2

                matrx[intron.junc1.index].append(intron.junc1_cov)
                matrx[intron.junc2.index].append(intron.junc2_cov)
                ex = intron.to_exon()
                exon_list.append(ex)
                del intron

        return tot_reads
    except ValueError as e:
        print(e)
        logging.debug('\t[%s]There are no reads in %s:%d-%d' % (info_msg, gne['chromosome'], gne['start'], gne['end']))
        return 0
