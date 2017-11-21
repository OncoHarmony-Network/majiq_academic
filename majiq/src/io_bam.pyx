from majiq.src.constants import *
import numpy as np
cimport numpy as np
from majiq.grimoire.exon cimport Intron
from majiq.grimoire.junction cimport Junction
from majiq.grimoire.junction import Junction
from majiq.src.io import dump_junctions
from majiq.src.multiproc import QueueMessage

from majiq.src.config import Config
from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from pysam.libcalignedsegment cimport PileupColumn, PileupRead
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

    return cross, jlist, read.reference_end


cdef inline bint __is_unique(AlignedSegment read):
    return not(read.flag & 0x100 == 0x100)


cdef inline int __get_num_reads(AlignedSegment read):
    return 1


cdef inline bint _match_strand(AlignedSegment read, str gene_strand):
    majiq_config = Config()
    res = True
    #print(read.is_reverse, read.flag, read.flag & 0x10, read.flag & 0x10 == 0x10, gene_strand, gene_strand == b'+',  gene_strand == '+')
    if majiq_config.strand_specific:
        #TODO: REMOVE
        if (read.flag & 0x10 == 0x10 and gene_strand == '+') or (read.flag & 0x10 == 0x00 and gene_strand == '-'):
            res = True
        else:
            res = False
    return res


cdef inline int _get_left_index(junc_start, read_start):
    return junc_start - (read_start + MIN_BP_OVERLAP)


cdef inline bint __valid_intron_read(AlignedSegment read):
    cdef bint is_cross
    cdef list junc_list
    cdef bint bo = __is_unique(read)

    is_cross, junc_list, rend = __cross_junctions(read)

    return bo and not is_cross and _match_strand(read, gstrand)


cdef str gstrand

cpdef int find_introns(str filename, dict list_introns, float intron_threshold, queue, str gname) except -1:

    cdef AlignedSegment read
    cdef AlignmentFile samfl
    cdef int readlen
    cdef bint is_cross
    cdef dict detected_introns = {}
    cdef list junc_list, chrom_list
    cdef int i_str, i_end, nchunks, chk_len, offs, rend
    cdef tuple introns
    cdef str strand, gne_id
    cdef np.ndarray cover

    samfl = open_rnaseq(filename)
    read_iter = samfl.fetch(until_eof=True)

    detected_introns = {}

    for read in read_iter:
        chrom = samfl.getrname(read.reference_id)
        is_cross, junc_list, rend = __cross_junctions(read)
        if not __is_unique(read) or is_cross or chrom not in list_introns.keys() :
            continue

        overlap = len(read.seq) - MIN_BP_OVERLAP
        for gne_id, strand, i_str, i_end, nchunks, chk_len in list_introns[chrom]:

            if read.pos <= (i_str - overlap) or not _match_strand(read, gene_strand=strand):
                break
            elif read.pos >= i_end:
                continue

            offs = max(0, int((read.pos - i_str) / chk_len))
            try:
                detected_introns[(gne_id, i_str, i_end, chk_len)][offs] += 1
            except KeyError:
                detected_introns[(gne_id, i_str, i_end, chk_len)] = np.zeros(shape=nchunks, dtype=np.float)
                detected_introns[(gne_id, i_str, i_end, chk_len)][offs] += 1

    for introns, cover in detected_introns.items():
        cover = cover/chk_len
        if np.any(cover < intron_threshold):
            continue

        qm = QueueMessage(QUEUE_MESSAGE_BUILD_INTRON, (introns[0], introns[1], introns[2], gname), 0)
        queue.put(qm, block=True)

    samfl.close()
    return 0



# cpdef int find_introns2(str filename, dict list_introns, float intron_threshold, queue, str gname) except -1:
#
#     cdef AlignmentFile samfl
#     cdef int nchunks, chunk_len
#     cdef bint b_included
#     cdef int intron_len, i_st, i_nd, val, ibin, num_bins = 10
#     cdef str chrom, strand
#     cdef np.ndarray intron_bins
#     cdef PileupColumn pile
#     cdef PileupRead xx
#     global gstrand
#
#     samfl = open_rnaseq(filename)
#
#
#     for gne_id, chrom, strand, i_st, i_nd in list_introns.keys():
#
#         intron_len = (i_nd - i_st)
#         nchunks = 1 if intron_len <= MIN_INTRON_LEN else num_bins
#
#         chunk_len = int(intron_len / nchunks)+1
#         b_included = False
#
#         for ii in range(nchunks):
#             lb = i_st + ii*chunk_len
#             ub = i_st + (ii+1)*chunk_len -1
#             ub = min(ub, i_nd)
#             try:
#                 gstrand = strand
#                 val = samfl.count(contig=chrom, start=lb, stop=ub, until_eof=True,
#                                   read_callback=__valid_intron_read, reference=None, end=None)
#             except ValueError as e:
#                 b_included = False
#                 break
#
#             val /= (ub-lb)
#             b_included = (val>=intron_threshold)
#             if not b_included:
#                 break
#
#         if b_included :
#             qm = QueueMessage(QUEUE_MESSAGE_BUILD_INTRON, (gne_id, i_st, i_nd, gname), 0)
#             queue.put(qm, block=True)
#             #list_introns[(gne_id, chrom, strand, i_st, i_nd)] += 1
#
#     samfl.close()
#     return 0


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
    cdef list junc_list, chrom_list
    cdef int junc_start, junc_end, rend
    cdef tuple jid, jid2
    cdef str strand

    chrom_list = [xx[0] for xx in in_jj]

    samfl = open_rnaseq(filename)
    read_iter = samfl.fetch(until_eof=True)


    for read in read_iter:
        is_cross, junc_list, rend = __cross_junctions(read)
        chrom = samfl.getrname(read.reference_id)
        if not __is_unique(read) or not is_cross or chrom not in chrom_list:
            continue

        readlen = len(read.seq)
        for junc_start, junc_end in junc_list:

            if (junc_start - read.pos >= readlen - MIN_BP_OVERLAP) or (junc_start - read.pos <= MIN_BP_OVERLAP) or \
                (junc_end - junc_start < MIN_JUNC_LENGTH):
                    continue
            if stranded:
                strand = '-' if read.is_reverse else '+'
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
                    if junc[0] in start_sp or junc[1] in end_sp:
                        found = True
                        qm = QueueMessage(QUEUE_MESSAGE_BUILD_JUNCTION, (gid, junc[0], junc[1], gname), 0)
                        queue.put(qm, block=True)
                    else:
                        possible_genes.append(gobj)
                    gidx +=1


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

        left_ind = _get_left_index(junc_start, r_start)
        if (junc_start, junc_end) in junctions:
            ''' update junction and add to list'''
            try:
                junc = junctions[(junc_start, junc_end)]
                junc.update_junction_read(nreads)
                if junc.index == 0:
                    junc.index = len(matrx)
                    matrx.append([0] * effective_len)

                matrx[junc.index][left_ind] += nreads
            except KeyError:
                # print (junctions.values()[0].gene_id,junc_start, junc_end)
                continue

cdef int __intronic_read(AlignedSegment read, junc_start, junc_end, list ref_pos, dict junctions,
                         int effective_len, list matrx) except -1:

    cdef long r_start = read.pos
    cdef int nreads = __get_num_reads(read)

    if junc_start in ref_pos[MIN_BP_OVERLAP:-MIN_BP_OVERLAP] and junc_end in ref_pos[MIN_BP_OVERLAP:-MIN_BP_OVERLAP]:
        try:
            junc = junctions[(junc_start, junc_end)]
            junc.update_junction_read(nreads)
            if junc.index == 0:
                junc.index = len(matrx)
                matrx.append([0] * effective_len)
            left_ind = _get_left_index(junc_start, r_start)
            try:
                matrx[junc.index][left_ind] += nreads
            except IndexError:
                print('INDEX ERROR', len(matrx), junc.index, left_ind)

        except KeyError:
            return 0


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
    cdef list junc_list, ref_pos
    cdef Intron intron

    try:

        majiq_config = Config()
        effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1
        read_iter = samfl.fetch(gne['chromosome'], gne['start'], gne['end'], multiple_iterators=False)

        for read in read_iter:
            is_cross, junc_list, end_r = __cross_junctions(read)
            unique = __is_unique(read)
            #print(read, _match_strand(read, gene_strand=gne['strand']), read.pos < gne['start'], unique)
            if not _match_strand(read, gene_strand=gne['strand']) or read.pos < gne['start'] or not unique:
                continue


            tot_reads += 1
            if is_cross:
              __junction_read(read, junc_list, majiq_config.readLen, effective_len, matrx, junctions)
            else:
                ref_pos = read.get_reference_positions()
                for intron in intron_list:
                    if read.pos >= intron.end:
                        break
                    elif end_r <= intron.start:
                        continue
                    #junc1
                    __intronic_read(read, intron.start-1, intron.start, ref_pos, junctions,
                                    effective_len, matrx)
                    #junc2
                    __intronic_read(read, intron.end, intron.end+1, ref_pos, junctions,
                                    effective_len, matrx)

        return tot_reads
    except ValueError as e:
        logging.debug('\t[%s]There are no reads in %s:%d-%d' % (info_msg, gne['chromosome'], gne['start'], gne['end']))
        return 0


## API
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


