from majiq.src.constants import *
import numpy as np
from majiq.grimoire.exon import detect_exons, new_exon_definition
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

cpdef long read_sam_or_bam(object gne, AlignmentFile samfl, object h5py_file,
                          str info_msg='0', object logging=None) except -1:
    _read_sam_or_bam(gne, samfl, h5py_file, info_msg=info_msg, logging=logging)


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef long _read_sam_or_bam2(object gne, AlignmentFile samfl, object h5py_file,
                           str info_msg='0', object logging=None) except -1:

    cdef unsigned int strt, end, r_start, junc_start, junc_end, readlen, nc, ng
    cdef str strand
    cdef AlignedSegment read
    cdef bint unique, found
    cdef float gc_content
    cdef bint bb
    cdef dict junctions = {xx.get_coordinates(): xx for xx in gne.get_all_junctions(filter=False)}
    cdef list ex_list
    cdef str chrom
    cdef int counter = 0

    try:
        strand = gne.get_strand()
        strt, end = gne.get_coordinates()
        # ex_list = gne.get_exon_list()
        chrom = gne.get_chromosome()
        majiq_config = Config()


        read_iter = samfl.fetch(chrom, strt, end, multiple_iterators=False)
        junc_k = samfl.find_introns(read_iter)

        for (junc_start, junc_end) in junc_k.keys():
            if junc_end - junc_start < MIN_JUNC_LENGTH :
                continue

            junc_start_reg = junc_start - (majiq_config.readLen - MIN_BP_OVERLAP)
            read_iter = samfl.fetch(chrom, junc_start_reg, junc_start - MIN_BP_OVERLAP, multiple_iterators=False)

            if (junc_start, junc_end) in junctions:
                ''' update junction and add to list'''
                junc = junctions[(junc_start, junc_end)]
            elif not majiq_config.non_denovo:
                bb = gne.check_antisense_junctions_hdf5(junc_start, junc_end, h5py_file)
                if not bb:
                    counter += 1
                    junc = Junction(junc_start, junc_end, None, None, gne.get_id(), retrieve=True)
                    junctions[(junc_start, junc_end)] = junc

            for read in read_iter:
                unique = __is_unique(read)

                r_start = read.pos
                is_cross, junc_list = __cross_junctions(read)
                if not _match_strand(read, gene_strand=strand) or r_start < strt or not unique or not is_cross:
                    continue

                if (junc_start, junc_end) in junc_list:
                    nreads = __get_num_reads(read)
                    nc = read.seq.count('C') + read.seq.count('c')
                    ng = read.seq.count('g') + read.seq.count('G')
                    gc_content = float(nc + ng) / float(len(read.seq))
                    readlen = len(read.seq)
                    if junc_start - r_start > readlen:
                        r_start_offset = junc_list[0][0] - r_start
                        r_start = junc_start - r_start_offset
                    junc.update_junction_read(nreads, r_start, gc_content, unique)


        if counter > 0:
            detect_exons(gne, junctions, None)
    except ValueError as e:
        print(e)
        logging.error('\t[%s]There are no reads in %s:%d-%d' % (info_msg, chrom, strt, end))


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef long _read_sam_or_bam(object gne, AlignmentFile samfl, object h5py_file,
                           str info_msg='0', object logging=None) except -1:

    cdef unsigned int strt, end, r_start, junc_start, junc_end, readlen, nc, ng
    cdef str strand
    cdef AlignedSegment read
    cdef bint unique, found
    cdef float gc_content
    cdef bint bb
    cdef dict junctions = {xx.get_coordinates(): xx for xx in gne.get_all_junctions(filter=False)}
    cdef list ex_list
    cdef str chrom
    cdef int counter = 0
    cdef int tot_reads = 0

    try:
        strand = gne.get_strand()
        strt, end = gne.get_coordinates()
        ex_list = gne.get_exon_list()
        chrom = gne.get_chromosome()
        majiq_config = Config()


        read_iter = samfl.fetch(chrom, strt, end, multiple_iterators=False)

        for read in read_iter:
            is_cross, junc_list = __cross_junctions(read)
            r_start = read.pos
            unique = __is_unique(read)
            if not _match_strand(read, gene_strand=strand) or r_start < strt or not unique:
                continue

            nreads = __get_num_reads(read)
            tot_reads += nreads

            if majiq_config.gcnorm:
                for ex_idx in range(len(ex_list)):
                    ex_start, ex_end = ex_list[ex_idx].get_coordinates()
                    if ex_start <= r_start <= ex_end:
                        ex_list[ex_idx].update_coverage(nreads)
                        break

            if not is_cross:
                continue

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

                if (junc_start, junc_end) in junctions:
                    ''' update junction and add to list'''
                    jj = junctions[(junc_start, junc_end)]
                    jj.update_junction_read(nreads, r_start, gc_content, unique)

                elif not majiq_config.non_denovo:
                    bb = gne.check_antisense_junctions_hdf5(junc_start, junc_end, h5py_file)
                    if not bb:
                        counter += 1
                        junc = Junction(junc_start, junc_end, None, None, gne.get_id(), retrieve=True)
                        junc.update_junction_read(nreads, r_start, gc_content, unique)
                        junctions[(junc_start, junc_end)] = junc

        gne.add_read_count(tot_reads)
        if counter > 0:
            detect_exons(gne, junctions, None)
    except ValueError as e:
        print(e)
        logging.error('\t[%s]There are no reads in %s:%d-%d' % (info_msg, chrom, strt, end))
    finally:
        gne.prepare_exons()



@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cpdef long rnaseq_intron_retention(object gne, AlignmentFile samfl, int chnk, object logging=None) except -1:

    cdef unsigned short num_bins = NUM_INTRON_BINS, nchunks
    cdef str strand = gne.get_strand()
    cdef AlignedSegment read
    cdef float gc_content
    cdef bint bb, is_cross, unique, intron_body_covered
    cdef int strt, end, r_start, intron_start, intron_end, readlen, nc, ng
    cdef int nreads, offset, ex1_end, ex2_start, intron_len, cov1, cov2
    cdef unsigned int intron_idx, num_positions, chunk_len
    cdef object junc1, junc2, ex, exon1, exon2

    intron_list = gne.get_all_introns()
    chrom = gne.get_chromosome()
    majiq_config = Config()
    for exon1, exon2 in intron_list:
        ex1_end = exon1.get_coordinates()[1]
        ex2_start = exon2.get_coordinates()[0]
        intron_start = ex1_end + 1
        intron_end = ex2_start - 1

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
            if not _match_strand(read, gene_strand=strand):
                continue

            unique = __is_unique(read)
            r_start = read.pos
            nreads = __get_num_reads(read)

            if not unique:
                # intron_idx = r_start - (ex1_end + 1)
                # if not (0 <= intron_idx <= intron_len):
                #     continue
                # bmap[intron_idx] = False
                continue

            if is_cross:
                jvals = [xx for xx, yy in junc_list if not (yy < intron_start or xx > intron_end)]
                if len(jvals) > 0:
                    continue

            nc = read.seq.count('C') + read.seq.count('c')
            ng = read.seq.count('g') + read.seq.count('G')
            gc_content = float(nc + ng) / float(len(read.seq))
            readlen = len(read.seq)
            offset = readlen - MIN_BP_OVERLAP

            if intron_start - r_start > readlen:
                r_start = intron_start - (readlen - MIN_BP_OVERLAP*2) - 1

            if r_start < ex1_end - MIN_BP_OVERLAP:
                if junc1 is None:
                    junc1 = Junction(ex1_end, intron_start, exon1, None, gne.get_id(), retrieve=True)
                junc1.update_junction_read(nreads, r_start, gc_content, unique)

            elif (ex2_start - offset - 1) < r_start < ex2_start:
                if junc2 is None:
                    junc2 = Junction(intron_end, ex2_start, exon2, None, gne.get_id(), retrieve=True)
                junc2.update_junction_read(nreads, r_start, gc_content, unique)

            else:
                # section 3
                if r_start <= ex1_end: continue
                intron_idx = r_start - (ex1_end + 1)
                # print (intron_idx, r_start, ex1_end)
                rel_start = int(intron_idx / chunk_len)
                indx = -1 if rel_start >= nchunks else rel_start
                # if not bmap[intron_idx]:
                #     bmap[intron_idx] = True
                intron_parts[indx] += nreads

        if junc1 is None or junc2 is None:
            continue

        cov1 = junc1.get_coverage().sum()
        cov2 = junc2.get_coverage().sum()

        # intron_parts /= chunk_len

        intron_body_covered = True

        if intron_len > 2 * (majiq_config.readLen - MIN_BP_OVERLAP):
            for ii in range(nchunks):
                #num_positions = np.count_nonzero(bmap[index_list[ii][0]:index_list[ii][1]])
                num_positions = chunk_len
                nii = intron_parts[ii]
                if nii == 0:
                    val = 0
                elif num_positions == 0:
                    continue
                else:
                    val = float(nii) / num_positions
                if val < majiq_config.min_intronic_cov:
                    intron_body_covered = False
                    break

        if cov1 >= majiq_config.min_denovo and cov2 >= majiq_config.min_denovo and intron_body_covered:
            exnum = new_exon_definition(intron_start, intron_end,
                                        junc1, junc2, gne, nondenovo=majiq_config.non_denovo,
                                        isintron=True)
            if exnum == -1:
                continue
            logging.debug("NEW INTRON RETENTION EVENT %s, %d-%d" % (gne.get_name(), intron_start, intron_end))
            junc1.add_donor(exon1)
            for ex in exon1.exonRead_list:
                st, end = ex.get_coordinates()
                if end == junc1.get_coordinates()[0]:
                    ex.add_5prime_junc(junc1)
                    break

            junc2.add_acceptor(exon2)
            for ex in exon2.exonRead_list:
                st, end = ex.get_coordinates()
                if st == junc2.get_coordinates()[1]:
                    ex.add_3prime_junc(junc2)
                    break
    gne.prepare_exons()


