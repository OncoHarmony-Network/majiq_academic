from majiq.src.constants import *
from majiq.grimoire.junction cimport Junction
cimport numpy as np
import numpy as np

cdef class Exon:
    def __init__(self, start, end, annot=False, intronic=False, db_idx=-1):
        self.start = start
        self.end = end
        self.annot = annot

        self.db_coords = (start, end)
        self.ib = set()
        self.ob = set()
        self.intron = intronic
        #self.db_idx = db_idx


cdef class Intron:

    num_bins = 10
    def __init__(self, start, end, annot=False, db_idx=-1):

        self.nchunks = 1 if (end-start) <= MIN_INTRON_LEN else self.num_bins
        self.start = start
        self.end = end
        self.chunk_len = int((end-start) / self.nchunks)
        self.skip = False

        self.junc1 = None
        self.junc2 = None
        self.annot = annot
        self.parts = np.zeros(shape=self.nchunks, dtype=float)
        self.db_idx = db_idx

# def exon_from_intron(gid, start, end, annot, junc1, junc2, list out_list):
#     ex = Exon(start, end, annot=annot)
#     ex.ib.add(Junction(start-1,  start, gid, -1, annot=annot))
#     ex.ob.add(Junction(end,  end+1, gid, -1, annot=annot))


cdef Exon exon_overlap(dict dd, int st, int end):
    cdef tuple kk
    cdef Exon vv

    for kk, vv in dd.items():
        if vv.start != EMPTY_COORD and vv.end != EMPTY_COORD and st < kk[1] and end > kk[0]:
            return vv


cdef int new_exon_definition(int start, int end, dict exon_dict, list out_list, Junction inbound_j, Junction outbound_j,
                             bint in_db=False):
    cdef Exon ex1, ex2
    cdef int new_exons = 0
    if end - start < 1:
        return 0

    ex1 = exon_overlap(exon_dict, start, end)

    if inbound_j.intronic and outbound_j.intronic:
        try:
            ex1 = exon_dict[(inbound_j.end, outbound_j.start)]
            ex2 = ex1

        except KeyError:
            return 0

    elif ex1 is None:

        if (end - start) <= MAX_DENOVO_DIFFERENCE:
            ex1 = Exon(start, end, annot=in_db)
            ex2 = ex1
            exon_dict[(start, end)] = ex1
            out_list.append(ex1)
            new_exons = 1
        else:
            try:
                ex1 = exon_dict[(start, start + 10)]
            except KeyError:
                ex1 = Exon(start, EMPTY_COORD, annot=in_db)
                exon_dict[(start, start + 10)] = ex1
                out_list.append(ex1)
                new_exons += 1

            try:
                ex2 = exon_dict[(end - 10, end)]
            except KeyError:
                ex2 = Exon(EMPTY_COORD, end, annot=in_db)
                exon_dict[(end - 10, end)] = ex2
                out_list.append(ex2)
                new_exons += 1


    else:
        ex2 = ex1
        if start < (ex1.start - MAX_DENOVO_DIFFERENCE):
            try:
                ex1 = exon_dict[(start, start + 10)]
            except KeyError:
                ex1 = Exon(start, EMPTY_COORD, annot=in_db)
                exon_dict[(start, start + 10)] = ex1
                out_list.append(ex1)
                new_exons += 1

        else:
            if start < ex1.start:
                ex1.start = start

        if end > (ex2.end + MAX_DENOVO_DIFFERENCE):
            try:
                ex2 = exon_dict[(end - 10, end)]
            except KeyError:
                ex2 = Exon(EMPTY_COORD, end, annot=in_db)
                exon_dict[(end - 10, end)] = ex2
                out_list.append(ex2)
                new_exons += 1
        else:
            if end > ex2.end:
                ex2.end = end

    if inbound_j is not None:
        ex1.ib.add(inbound_j)
        inbound_j.acceptor = ex1

    if outbound_j is not None:
        ex2.ob.add(outbound_j)
        outbound_j.donor = ex2

    return new_exons


def detect_exons(dict junction_dict, list exon_list):
    cdef int opened, new_exons = 0
    cdef list opened_exon = []
    cdef Junction last_5prime = None
    cdef Junction first_3prime = None
    cdef Exon ex
    cdef dict exon_dict = {(ex.start, ex.end): ex for ex in exon_list}
    cdef list junction_list = []
    cdef tuple kk
    cdef Junction jj, jj2
    cdef long coord
    cdef bint jtype

    for kk, jj in junction_dict.items():
        if kk[0]>0:
            junction_list.append((kk[0], True, jj))
        if kk[1]>0:
            junction_list.append((kk[1], False, jj))
    junction_list.sort(key=lambda jj: (jj[0], -jj[1]))

    for (coord, jtype, jj) in junction_list:
        if not (jj.is_reliable() or jj.annotated):
            continue
        opened = len(opened_exon)
        if jtype:
            if opened > 0:

                # print('ST1', opened_exon[-1].end, coord)
                new_exons += new_exon_definition(opened_exon[-1].end, coord, exon_dict, exon_list, opened_exon[-1], jj, in_db=False)
                opened_exon.pop()
            elif opened == 0:
                if first_3prime is None:
                    # print('ST2', opened_exon[-1].end, coord)
                    new_exons += new_exon_definition(jj.start-10, jj.start, exon_dict, exon_list, None, jj)
                    #new_exons += __half_exon('5prime', jj)
                else:
                    new_exons += new_exon_definition(first_3prime.end, coord, exon_dict, exon_list, first_3prime, jj, in_db=False)
            last_5prime = jj
            # end elif opened
        else:
            if opened > 0:
                if last_5prime is not None:
                    for jj2 in opened_exon:
                        # print('ST4', opened_exon[-1].end, coord)
                        new_exons += new_exon_definition(jj2.end, last_5prime.start, exon_dict, exon_list, jj2, last_5prime, in_db=False)
                    last_5prime = None
                    opened_exon = []
                    first_3prime = jj
            else:
                last_5prime = None
                first_3prime = jj
            opened_exon.append(jj)

    for jj in opened_exon:
        # print('ST5', opened_exon[-1].end, coord)
        new_exons += new_exon_definition(jj.end, jj.end+10, exon_dict, exon_list, jj, None, in_db=False)

    exon_list = sorted(exon_dict.values(), key=lambda ex: ex.start)
    return new_exons


cpdef expand_introns(str gne_id, list list_introns, list list_exons, dict dict_junctions, int default_index=-1):

    cdef list new_list = []
    cdef Exon ex, donor_ex, acceptor_ex, ir_ex
    cdef Junction jj
    cdef Intron intron
    cdef int ex_start, ex_end

    for intron in list_introns:

        for ex in list_exons:
            if ex.intron:
                continue

            ex_start = ex.start if ex.start != EMPTY_COORD else ex.end-1
            ex_end = ex.end if ex.end != EMPTY_COORD else ex.start+1

            if ex_start <= intron.start and ex_end >= intron.end:
                intron.skip = True
                break

            if ex.end != EMPTY_COORD and ex_start<= intron.start <= (ex_end+1) :
                intron.start = ex_end + 1
                donor_ex = ex

            if ex_start > intron.start and ex_end < intron.end:

                if ex.end != EMPTY_COORD and ex.start != EMPTY_COORD :
                    ir_ex = Exon(intron.start, ex.end-1, annot=intron.annot, intronic=True)
                    new_list.append(ir_ex)

                    jj = Junction(intron.start-1, intron.start, gne_id, default_index, intron=True, annot=intron.annot)
                    jj.donor = donor_ex
                    jj.acceptor = ir_ex
                    ir_ex.ib.add(jj)
                    donor_ex.ob.add(jj)
                    dict_junctions[(intron.start-1, intron.start)] = jj

                    jj = Junction(ex.start-1, ex.start, gne_id, default_index, intron=True, annot=intron.annot)
                    jj.donor = ir_ex
                    jj.acceptor = ex
                    ir_ex.ob.add(jj)
                    ex.ib.add(jj)
                    dict_junctions[(ex.start-1, ex.start)] = jj

                    intron.start = ex.end +1
                    donor_ex = ex

                elif ex.end == EMPTY_COORD:
                    intron.end = ex_start - 1
                    acceptor_ex = ex
                elif ex.start == EMPTY_COORD:
                    intron.start = ex_end + 1
                    donor_ex = ex

            if ex.start != EMPTY_COORD and (ex_start-1)<= intron.end <= ex_end:
                intron.end = ex_start - 1
                acceptor_ex = ex

        if intron.skip:
            continue
        ir_ex = Exon(intron.start, intron.end, annot=intron.annot, intronic=True)
        new_list.append(ir_ex)

        jj = Junction(intron.start-1, intron.start, gne_id, default_index, intron=True, annot=intron.annot)

        jj.donor = donor_ex
        jj.acceptor = ir_ex
        ir_ex.ib.add(jj)
        donor_ex.ob.add(jj)
        dict_junctions[(intron.start-1, intron.start)] = jj

        jj = Junction(intron.end, intron.end+1, gne_id, default_index, intron=True, annot=intron.annot)
        jj.donor = ir_ex
        jj.acceptor = acceptor_ex
        ir_ex.ob.add(jj)
        acceptor_ex.ib.add(jj)
        dict_junctions[(intron.end, intron.end+1)] = jj

    list_exons.extend(new_list)