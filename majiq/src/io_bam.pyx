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


cpdef int  find_introns(str filename, dict list_introns, float intron_threshold, int group_indx) except -1:


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
        chunk_len = int(intron_len/ num_bins)
        b_included = True
        for ii in range(nchunks):
            lb = i_st + ii*chunk_len
            ub = i_st + (ii+1)*chunk_len -1
            ub = min(ub, i_nd)
            val = samfl.count(contig=chrom, start=lb, stop=ub, until_eof=True, read_callback=__is_unique,
                              reference=None, end=None)

            val /= (ub-lb)
            b_included = b_included and (val>=intron_threshold)

        if b_included :
            list_introns[(gne_id, chrom, strand, i_st, i_nd)] += 1

    samfl.close()
    return 0

cdef inline int __read_STAR_junc_file(int fidx, str filename, set out_jj, dict out_dd) except -1:

    with open(filename, 'r') as fp:
        for ln in fp.readlines():
            tab = ln.strip().split()
            # if (tab[0], tab[1], tab[2]) in in_juncs:
            #     continue
            try:
                out_dd[(tab[0], tab[1], tab[2])] += tab[6]
            except KeyError:
                out_dd[(tab[0], tab[1], tab[2])] = tab[6]
            out_jj.append((tab[0], tab[1], tab[2], tab[3]))


cdef int __read_juncs_from_bam(int fidx, str filename, list out_jj, dict out_dd) except -1:

    cdef AlignedSegment read
    cdef AlignmentFile samfl
    cdef int readlen
    cdef bint is_cross
    cdef list junc_list
    cdef int junc_start, junc_end, rend

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
            # if (chrom, junc_start, junc_end) in out_jj:
            #     continue
            #out_jj.add((read.reference, junc_start, junc_end))
            try:
                out_dd[(chrom, junc_start, junc_end)] += 1
            except KeyError:
                out_dd[(chrom, junc_start, junc_end)] = 1

def read_juncs(db_f, list file_list, set in_juncs, dict junctions, dict all_genes, int min_denovo, object q):

    cdef str ln, fname, chrom='', jjid
    cdef list tab
    cdef list jj_list = []
    cdef dict jj_dd = {}
    cdef tuple jj
    cdef int gidx=0, ngenes=0
    cdef object fp

    for fidx, (bb, fname) in enumerate(file_list):
        if bb:
            __read_STAR_junc_file(fidx, fname, jj_list, jj_dd)
        else:
            __read_juncs_from_bam(fidx, fname, jj_list, jj_dd)
            jj_list.extend(list(jj_dd.keys()))

    jj_list.sort(key=lambda xx: (xx[0], xx[1], xx[2]))
    jj_dd = {kk: float(vv)/len(file_list) for kk, vv in jj_dd.items()}

    for jj in jj_list:
        if jj_dd[(jj[0], jj[1], jj[2])] < min_denovo:
            continue

        if chrom != jj[0]:
            chrom = jj[0]
            ngenes = len(all_genes[chrom])

        while gidx < ngenes :
            gstart, gend, gid = all_genes[chrom][gidx]
            if gend < jj[1]:
                gidx += 1
                continue
            if gstart > jj[2]:
                break

            if gstart<=jj[1] and gend>= jj[2]:
                jjid = '%s/junctions/%s-%s' % (gid, jj[1], jj[2])
                if jjid not in in_juncs:
                    in_juncs.add((jj[0], jj[1], jj[2]))
                    # qm = QueueMessage(QUEUE_MESSAGE_BUILD_JUNCTION, (gid,jj[1],  jj[2]), 0)
                    # q.put(qm, block=True)
                    dump_junctions(db_f, gid, jj[1], jj[2], None, annot=False)
                junctions[gid][(jj[1], jj[2])] = Junction(jj[1],  jj[2], gid, -1, annot=False)
                junctions[gid][(jj[1], jj[2])].update_junction_read(jj_dd[(jj[0], jj[1], jj[2])])

                break
    return 0



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
        logging.error('\t[%s]There are no reads in %s:%d-%d' % (info_msg, gne['chrom'], gne['start'], gne['end']))
        return 0









#
#
#
#
#
#
#
#
#
#
#
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# cdef int _read_sam_or_bam2(object gne, AlignmentFile samfl, list matrx, list out_junctions, str info_msg='0',
#                           object logging=None) except -1:
#
#     cdef unsigned int r_start, junc_start, junc_end, readlen, nc, ng
#     cdef AlignedSegment read
#     cdef bint unique, found
#     cdef float gc_content
#     cdef bint bb
#     cdef dict junctions = {} #{xx.get_coordinates(): xx for xx in gne.get_all_junctions(filter=False)}
#     cdef int counter = 0
#     cdef int tot_reads = 0
#     cdef int effective_len = 0
#     cdef Junction junc
#     cdef list junc_list
#
#
#     try:
#         majiq_config = Config()
#         effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1
#         read_iter = samfl.fetch(gne['chromosome'], gne['start'], gne['end'], multiple_iterators=False)
#
#         for read in read_iter:
#             is_cross, junc_list = __cross_junctions(read)
#             r_start = read.pos
#             unique = __is_unique(read)
#             if not _match_strand(read, gene_strand=gne['strand']) or r_start < gne['start'] or not unique:
#                 continue
#
#             nreads = __get_num_reads(read)
#             tot_reads += nreads
#             nc = read.seq.count('C') + read.seq.count('c')
#             ng = read.seq.count('g') + read.seq.count('G')
#             gc_content = float(nc + ng) / float(len(read.seq))
#             readlen = len(read.seq)
#
#             if not is_cross:
#                 continue
#
#             if majiq_config.gcnorm:
#                 pass
#
#             for (junc_start, junc_end) in junc_list:
#
#                 if junc_start - r_start > readlen:
#                     r_start_offset = junc_list[0][0] - r_start
#                     r_start = junc_start - r_start_offset
#
#                 if (junc_start - r_start >= readlen - MIN_BP_OVERLAP) or (junc_start - r_start <= MIN_BP_OVERLAP) or \
#                         (junc_end - junc_start < MIN_JUNC_LENGTH):
#                     continue
#
#                 left_ind = majiq_config.readLen - (junc.start - r_start) - MIN_BP_OVERLAP + 1
#                 if (junc_start, junc_end) in junctions:
#                     ''' update junction and add to list'''
#                     junc = junctions[(junc_start, junc_end)]
#                     junc.update_junction_read(nreads)
#                     matrx[junc.index][left_ind] += nreads
#
#                 elif not majiq_config.non_denovo:
#
#                     #TODO: fix antisense
#                     bb = False #gne.check_antisense_junctions_hdf5(junc_start, junc_end, h5py_file)
#                     if not bb:
#
#                         counter += 1
#                         junc = Junction(junc_start, junc_end, gne['id'], cov_idx=len(matrx))
#                         junc.update_junction_read(nreads)
#                         matrx.append([0]*effective_len)
#                         matrx[junc.index][left_ind] += nreads
#                         junctions[(junc_start, junc_end)] = junc
#                         out_junctions.append(junc)
#
#         return tot_reads
#     except ValueError as e:
#         print(e)
#         logging.error('\t[%s]There are no reads in %s:%d-%d' % (info_msg, gne['chrom'], gne['start'], gne['end']))
#         return 0
#
#
# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# cdef long _rnaseq_intron_retention(dict gne, list list_exons, list matrx, list out_junctions,
#                                    object logging=None) except -1:
#
#     # cdef unsigned short num_bins = NUM_INTRON_BINS, nchunks
#     # cdef str strand = gne['strand']
#     # cdef str chrom = gne['chromosome']
#     # cdef AlignedSegment read
#     # cdef float gc_content
#     # cdef bint is_cross, unique, intron_body_covered
#     # cdef int nreads, offset, intron_len, strt, end, r_start, intron_start, intron_end, readlen, nc, ng
#     # cdef unsigned int intron_idx, num_positions, chunk_len, xx, yy
#     # cdef Junction junc1, junc2
#     # cdef Exon exon1, exon2
#     # cdef list jvals, junc1_cov, junc2_cov
#
#     cdef int effective_len, ex_idx
#     cdef object majiq_config
#
#     majiq_config = Config()
#     effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1
#     for exp_idx, samfl in majiq_config.sam_list:
#
#     for ex_idx in list_exons[:-1]:
#         exon1 = list_exons[ex_idx]
#         exon2 = list_exons[ex_idx+1]
#         intron_start = exon1.end + 1
#         intron_end = exon2.start - 1
#         intron_len = intron_end - intron_start
#         if intron_len <= 0:
#             continue
#         try:
#             read_iter = samfl.fetch(chrom, intron_start + MIN_BP_OVERLAP, intron_end - MIN_BP_OVERLAP,
#                                     multiple_iterators=False)
#         except ValueError:
#             continue
#
#         nchunks = 1 if intron_len <= MIN_INTRON_LEN else num_bins
#
#         # we want to take just the middle part not the reads that are crossing the junctions
#         # since 8 is the overlapping number of nucleotites we accept, the inner part is the
#         # real intron size - (readlen-8)/*start part*/ - (readlen-8)/*end part*/
#
#         chunk_len = int(intron_len / nchunks)
#
#         # bmap = np.ones(shape=intron_len, dtype=np.bool)
#         index_list = []
#         for ii in range(nchunks):
#             start = ii * chunk_len
#             end = min(intron_len, (ii + 1) * chunk_len)
#             index_list.append((start, end))
#
#         intron_parts = np.zeros(shape=nchunks, dtype=np.float)
#         junc1 = None
#         junc2 = None
#
#         for read in read_iter:
#             is_cross, junc_list = __cross_junctions(read)
#             r_start = read.pos
#             unique = __is_unique(read)
#             if not _match_strand(read, gene_strand=gne['strand']) or r_start < gne['start'] or not unique:
#                 continue
#             nreads = __get_num_reads(read)
#
#             if is_cross:
#                 jvals = [xx for xx, yy in junc_list if not (yy < intron_start or xx > intron_end)]
#                 if len(jvals) > 0:
#                     continue
#
#             nc = read.seq.count('C') + read.seq.count('c')
#             ng = read.seq.count('g') + read.seq.count('G')
#             gc_content = float(nc + ng) / float(len(read.seq))
#             readlen = len(read.seq)
#             offset = readlen - MIN_BP_OVERLAP
#
#             left_ind = majiq_config.readLen - (exon1.end - r_start) - MIN_BP_OVERLAP + 1
#             if intron_start - r_start > readlen:
#                 r_start = intron_start - (readlen - MIN_BP_OVERLAP*2) - 1
#
#             if r_start < exon1.end - MIN_BP_OVERLAP:
#                 if junc1 is None:
#                     junc1 = Junction(exon1.end, intron_start, gne['id'], cov_idx=len(matrx), intron=True)
#                     junc1_cov = [0] * effective_len
#
#                 junc1.update_junction_read(nreads)
#
#             elif (exon2.start - offset - 1) < r_start < exon2.start:
#                 if junc2 is None:
#                     junc2 = Junction(intron_end, exon2.start, gne['id'], cov_idx=len(matrx), intron=True)
#                     junc2_cov = [0] * effective_len
#                 junc2.update_junction_read(nreads)
#
#             else:
#                 # section 3
#                 if r_start <= exon1.end: continue
#                 intron_idx = r_start - (exon1.end + 1)
#                 rel_start = int(intron_idx / chunk_len)
#                 indx = -1 if rel_start >= nchunks else rel_start
#                 intron_parts[indx] += nreads
#
#         if junc1 is None or junc2 is None:
#             continue
#
#         intron_body_covered = True
#         if intron_len > 2 * (majiq_config.readLen - MIN_BP_OVERLAP):
#             for ii in range(nchunks):
#                 num_positions = chunk_len
#                 if intron_parts[ii] == 0:
#                     val = 0
#                 elif num_positions == 0:
#                     continue
#                 else:
#                     val = float(intron_parts[ii]) / num_positions
#                 if val < majiq_config.min_intronic_cov:
#                     intron_body_covered = False
#                     break
#
#         if (junc1.nreads >= majiq_config.min_denovo and
#             junc2.nreads >= majiq_config.min_denovo and
#             intron_body_covered):
#
#             out_junctions.append(junc1)
#             out_junctions.append(junc2)
#
#             matrx[junc1.index].append(junc1_cov)
#             matrx[junc2.index].append(junc2_cov)
#             del junc1
#             del junc2
#
#     return 0
