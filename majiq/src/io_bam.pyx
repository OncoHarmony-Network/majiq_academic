from majiq.src.constants import *
import numpy as np
from majiq.grimoire.exon cimport Exon
from majiq.grimoire.junction cimport Junction
from majiq.grimoire.junction import Junction

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

    try:
        list_junc = read.opt('jI')
        if len(list_junc) > 1 and list_junc[0] != -1:
            for idx in range(0, len(list_junc), 2):
                junc_start = int(list_junc[idx]) - 1
                junc_end = int(list_junc[idx + 1]) + 1
                jlist.append((junc_start, junc_end))
            # end for idx ...
            cross = True
            # end else ...
    except KeyError:
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

    return cross, jlist

def __cross_junctions2(read):
    """
      This part will parse the jI from STAR
      # Column 17: jI:B:I,Start1,End1,Start2,End2,... Start and End of introns for all junctions (1- based)
      # jI:B:I,-1 means that the read doesn't cross any junction
    """

    jlist = []
    cross = False
    try:
        list_junc = read.opt('jI')
        if len(list_junc) > 1 and list_junc[0] != -1:
            for idx in range(0, len(list_junc), 2):
                junc_start = int(list_junc[idx]) - 1
                junc_end = int(list_junc[idx + 1]) + 1
                jlist.append((junc_start, junc_end))
            # end for idx ...
            cross = True
            # end else ...
    except KeyError:
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

    return cross, jlist

cdef inline int read_juncs_from_bam(str filename, set in_juncs, list out_jj, dict out_dd) except -1:
    cdef AlignedSegment read
    cdef AlignmentFile samfl

    samfl = open_rnaseq(filename)
    read_iter = samfl.fetch(until_eof=True)

    jj_set = set()

    for read in read_iter:
        is_cross, junc_list = __cross_junctions(read)
        if not __is_unique(read) or not is_cross:
            continue
        readlen = len(read.seq)
        for junc_start, junc_end in junc_list:
            if (junc_start - read.pos >= readlen - MIN_BP_OVERLAP) or (junc_start - read.pos <= MIN_BP_OVERLAP) or \
                (junc_end - junc_start < MIN_JUNC_LENGTH) or ((read.reference,junc_start, junc_end) in in_juncs):
                    continue

            out_jj.add((read.reference,junc_start, junc_end))
            try:
                out_dd[(read.reference,junc_start, junc_end)] += 1
            except KeyError:
                out_dd[(read.reference,junc_start, junc_end)] = 1



cdef inline bint __is_unique(AlignedSegment read):
    return not(read.flag & 0x100 == 0x100)


cdef inline int __get_num_reads(AlignedSegment read):
    return 1


cdef inline bint _match_strand(AlignedSegment read, str gene_strand):
    majiq_config = Config()
    res = True
    if majiq_config.strand_specific:
        #TODO: REMOVE
        print (gene_strand, b'+')
        if (read.flag & 0x10 == 0x10 and gene_strand == b'+') or (read.flag & 0x10 == 0x00 and gene_strand == b'-'):
            res = True
        else:
            res = False
    return res


cpdef AlignmentFile open_rnaseq(str samfile):
    return pysam.Samfile(samfile, "rb")


cpdef long close_rnaseq(AlignmentFile samfl) except -1:
    samfl.close()

cpdef int read_sam_or_bam(object gne, AlignmentFile samfl, object h5py_file, list out_junctions,
                          str info_msg='0', object logging=None) except -1:
    cdef int res
    res = _read_sam_or_bam(gne, samfl, h5py_file, info_msg=info_msg, out_junctions=out_junctions, logging=logging)
    return res


cpdef long rnaseq_intron_retention(dict gne, list list_exons, AlignmentFile samfl, list matrx, list out_junctions,
                                   object logging=None) except -1:

    cdef int res
    res = _rnaseq_intron_retention(gne, list_exons, samfl, matrx, out_junctions, logging=None)
    return res

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef int _read_sam_or_bam(object gne, AlignmentFile samfl, list matrx, list out_junctions, str info_msg='0',
                          object logging=None) except -1:

    cdef unsigned int r_start, junc_start, junc_end, readlen, nc, ng
    cdef AlignedSegment read
    cdef bint unique, found
    cdef float gc_content
    cdef bint bb
    cdef dict junctions = {} #{xx.get_coordinates(): xx for xx in gne.get_all_junctions(filter=False)}
    cdef int counter = 0
    cdef int tot_reads = 0
    cdef int effective_len = 0
    cdef Junction junc
    cdef list junc_list


    try:
        majiq_config = Config()
        effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1
        read_iter = samfl.fetch(gne['chromosome'], gne['start'], gne['end'], multiple_iterators=False)

        for read in read_iter:
            is_cross, junc_list = __cross_junctions(read)
            r_start = read.pos
            unique = __is_unique(read)
            if not _match_strand(read, gene_strand=gne['strand']) or r_start < gne['start'] or not unique:
                continue

            nreads = __get_num_reads(read)
            tot_reads += nreads

            if not is_cross:
                continue

            if majiq_config.gcnorm:
                pass

            nc = read.seq.count('C') + read.seq.count('c')
            ng = read.seq.count('g') + read.seq.count('G')
            gc_content = float(nc + ng) / float(len(read.seq))
            readlen = len(read.seq)
            for (junc_start, junc_end) in junc_list:

                if junc_start - r_start > readlen:
                    r_start_offset = junc_list[0][0] - r_start
                    r_start = junc_start - r_start_offset

                if (junc_start - r_start >= readlen - MIN_BP_OVERLAP) or (junc_start - r_start <= MIN_BP_OVERLAP) or \
                        (junc_end - junc_start < MIN_JUNC_LENGTH):
                    continue

                left_ind = majiq_config.readLen - (junc.start - r_start) - MIN_BP_OVERLAP + 1
                if (junc_start, junc_end) in junctions:
                    ''' update junction and add to list'''
                    junc = junctions[(junc_start, junc_end)]
                    junc.update_junction_read(nreads)
                    matrx[junc.index][left_ind] += nreads

                elif not majiq_config.non_denovo:

                    #TODO: fix antisense
                    bb = False #gne.check_antisense_junctions_hdf5(junc_start, junc_end, h5py_file)
                    if not bb:

                        counter += 1
                        junc = Junction(junc_start, junc_end, gne['id'], cov_idx=len(matrx))
                        junc.update_junction_read(nreads)
                        matrx.append([0]*effective_len)
                        matrx[junc.index][left_ind] += nreads
                        junctions[(junc_start, junc_end)] = junc
                        out_junctions.append(junc)

        return tot_reads
    except ValueError as e:
        print(e)
        logging.error('\t[%s]There are no reads in %s:%d-%d' % (info_msg, gne['chrom'], gne['start'], gne['end']))
        return 0


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef long _rnaseq_intron_retention(dict gne, list list_exons, list matrx, list out_junctions,
                                   object logging=None) except -1:

    # cdef unsigned short num_bins = NUM_INTRON_BINS, nchunks
    # cdef str strand = gne['strand']
    # cdef str chrom = gne['chromosome']
    # cdef AlignedSegment read
    # cdef float gc_content
    # cdef bint is_cross, unique, intron_body_covered
    # cdef int nreads, offset, intron_len, strt, end, r_start, intron_start, intron_end, readlen, nc, ng
    # cdef unsigned int intron_idx, num_positions, chunk_len, xx, yy
    # cdef Junction junc1, junc2
    # cdef Exon exon1, exon2
    # cdef list jvals, junc1_cov, junc2_cov

    cdef int effective_len, ex_idx
    cdef object majiq_config

    majiq_config = Config()
    effective_len = (majiq_config.readLen - 2*MIN_BP_OVERLAP) + 1
    for exp_idx, samfl in majiq_config.sam_list:

    for ex_idx in list_exons[:-1]:
        exon1 = list_exons[ex_idx]
        exon2 = list_exons[ex_idx+1]
        intron_start = exon1.end + 1
        intron_end = exon2.start - 1
        intron_len = intron_end - intron_start
        if intron_len <= 0:
            continue
        try:
            read_iter = samfl.fetch(chrom, intron_start + MIN_BP_OVERLAP, intron_end - MIN_BP_OVERLAP,
                                    multiple_iterators=False)
        except ValueError:
            continue

        nchunks = 1 if intron_len <= MIN_INTRON_LEN else num_bins

        # we want to take just the middle part not the reads that are crossing the junctions
        # since 8 is the overlapping number of nucleotites we accept, the inner part is the
        # real intron size - (readlen-8)/*start part*/ - (readlen-8)/*end part*/

        chunk_len = int(intron_len / nchunks)

        # bmap = np.ones(shape=intron_len, dtype=np.bool)
        index_list = []
        for ii in range(nchunks):
            start = ii * chunk_len
            end = min(intron_len, (ii + 1) * chunk_len)
            index_list.append((start, end))

        intron_parts = np.zeros(shape=nchunks, dtype=np.float)
        junc1 = None
        junc2 = None

        for read in read_iter:
            is_cross, junc_list = __cross_junctions(read)
            r_start = read.pos
            unique = __is_unique(read)
            if not _match_strand(read, gene_strand=gne['strand']) or r_start < gne['start'] or not unique:
                continue
            nreads = __get_num_reads(read)

            if is_cross:
                jvals = [xx for xx, yy in junc_list if not (yy < intron_start or xx > intron_end)]
                if len(jvals) > 0:
                    continue

            nc = read.seq.count('C') + read.seq.count('c')
            ng = read.seq.count('g') + read.seq.count('G')
            gc_content = float(nc + ng) / float(len(read.seq))
            readlen = len(read.seq)
            offset = readlen - MIN_BP_OVERLAP

            left_ind = majiq_config.readLen - (exon1.end - r_start) - MIN_BP_OVERLAP + 1
            if intron_start - r_start > readlen:
                r_start = intron_start - (readlen - MIN_BP_OVERLAP*2) - 1

            if r_start < exon1.end - MIN_BP_OVERLAP:
                if junc1 is None:
                    junc1 = Junction(exon1.end, intron_start, gne['id'], cov_idx=len(matrx), intron=True)
                    junc1_cov = [0] * effective_len

                junc1.update_junction_read(nreads)

            elif (exon2.start - offset - 1) < r_start < exon2.start:
                if junc2 is None:
                    junc2 = Junction(intron_end, exon2.start, gne['id'], cov_idx=len(matrx), intron=True)
                    junc2_cov = [0] * effective_len
                junc2.update_junction_read(nreads)

            else:
                # section 3
                if r_start <= exon1.end: continue
                intron_idx = r_start - (exon1.end + 1)
                rel_start = int(intron_idx / chunk_len)
                indx = -1 if rel_start >= nchunks else rel_start
                intron_parts[indx] += nreads

        if junc1 is None or junc2 is None:
            continue

        intron_body_covered = True
        if intron_len > 2 * (majiq_config.readLen - MIN_BP_OVERLAP):
            for ii in range(nchunks):
                num_positions = chunk_len
                if intron_parts[ii] == 0:
                    val = 0
                elif num_positions == 0:
                    continue
                else:
                    val = float(intron_parts[ii]) / num_positions
                if val < majiq_config.min_intronic_cov:
                    intron_body_covered = False
                    break

        if (junc1.nreads >= majiq_config.min_denovo and
            junc2.nreads >= majiq_config.min_denovo and
            intron_body_covered):

            out_junctions.append(junc1)
            out_junctions.append(junc2)

            matrx[junc1.index].append(junc1_cov)
            matrx[junc2.index].append(junc2_cov)
            del junc1
            del junc2

    return 0
