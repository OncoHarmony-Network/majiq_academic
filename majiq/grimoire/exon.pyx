from majiq.src.constants import *
from majiq.grimoire.junction cimport Junction
cimport numpy as np
import numpy as np

cdef class Exon:
    def __init__(self, start, end, annot=False):
        self.start = start
        self.end = end
        self.annot = annot

        self.db_coords = (start, end)
        self.ib = set()
        self.ob = set()
        self.intron = False


cdef class Intron:

    num_bins = 10
    def __init__(self, start, end, annot=False):

        self.nchunks = 1 if (end-start) <= MIN_INTRON_LEN else self.num_bins
        self.start = start
        self.end = end
        self.chunk_len = int((end-start) / self.nchunks)

        self.junc1 = None
        self.junc2 = None
        self.annot = annot
        self.parts = np.zeros(shape=self.nchunks, dtype=float)

    cdef Exon to_exon(self):
        ex = Exon(self.start, self.end, annot=False)
        ex.ib.add(self.junc1)
        ex.ob.add(self.junc2)
        ex.intron = True
        return ex




def new_intron_definition():
    pass


def exon_overlap(dd, st, end):
    for kk, vv in dd.items():
        if st <= kk[1] and end >= kk[0]:
            return vv


def new_exon_definition(start, end, exon_dict, inbound_j, outbound_j, in_db=False):
    if end - start < 1:
        return 0
    ex1 = exon_overlap(exon_dict, start, end)
    count = 0
    if ex1 is None:
        if (end - start) <= MAX_DENOVO_DIFFERENCE:
            ex1 = Exon(start, end, annot=in_db)
            ex2 = ex1
            exon_dict[(start, end)] = ex1
            count = 1
        else:
            ex1 = Exon(start, EMPTY_COORD, annot=in_db)
            exon_dict[(start, start + 2)] = ex1

            ex2 = Exon(EMPTY_COORD, end, annot=in_db)
            exon_dict[(end - 10, end)] = ex2
            count = 2

    else:
        ex2 = ex1
        if start < (ex1.start - MAX_DENOVO_DIFFERENCE):
            ex1 = Exon(start, EMPTY_COORD, annot=in_db)
            exon_dict[(start, start + 10)] = ex1
            count = 1
        else:
            if start < ex1.start:
                ex1.start = start

        if end > (ex2.end + MAX_DENOVO_DIFFERENCE):
            ex2 = Exon(EMPTY_COORD, end, annot=in_db)
            exon_dict[(end - 10, end)] = ex2
            count += 1
        else:
            if end > ex2.end:
                ex2.end = end

    if inbound_j is not None:
        ex1.ib.add(inbound_j)
        inbound_j.acceptor = ex1

    if outbound_j is not None:
        ex2.ob.add(outbound_j)
        outbound_j.donor = ex2

    return count


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
        if not jj.intronic:
            junction_list.append((kk[0], True, jj))
            junction_list.append((kk[1], False, jj))
    junction_list.sort(key=lambda jj: jj[0])

    for (coord, jtype, jj) in junction_list:
        if not (jj.is_reliable() or jj.annotated):
            continue
        opened = len(opened_exon)
        if jtype:
            if opened > 0:
                new_exons += new_exon_definition(opened_exon[-1].end, coord, exon_dict, opened_exon[-1], jj, in_db=False)
                opened_exon.pop()
            elif opened == 0:
                if first_3prime is None:
                    new_exons += new_exon_definition(jj.start-10, jj.start, exon_dict, None, jj)
                    #new_exons += __half_exon('5prime', jj)
                else:
                    new_exons += new_exon_definition(first_3prime.end, coord, exon_dict, first_3prime, jj, in_db=False)
            last_5prime = jj
            # end elif opened
        else:
            if opened > 0:
                if last_5prime is not None:
                    for jj2 in opened_exon:
                        new_exons += new_exon_definition(jj2.end, last_5prime.start, exon_dict, jj2, last_5prime, in_db=False)
                    last_5prime = None
                    opened_exon = []
                    first_3prime = jj
            else:
                last_5prime = None
                first_3prime = jj
            opened_exon.append(jj)

    for jj in opened_exon:
        new_exons += new_exon_definition(jj.end, jj.end+10, exon_dict, jj, None, in_db=False)

    exon_list = sorted(exon_dict.values(), key=lambda ex: ex.start)
    return new_exons