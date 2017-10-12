
cdef class Junction:

    def __init__(self, start, end, gene_id, cov_idx, annot=False, intron=False):
        self.start = start
        self.end = end
        self.gene_id = gene_id
        self.index = cov_idx

        self.nreads = 0
        self.acceptor = None
        self.donor = None
        self.annot = annot
        self.intronic = intron

    def update_junction_read(self, int n):
        self.nreads += n

    def is_reliable(self):
        return True

    def reset(self):
        self.nreads = 0
        self.index = 0



def create_junction(start, end, gene_id, cov_idx, annot=False, intron=False):
    return Junction(start, end, gene_id, cov_idx, annot, intron)
