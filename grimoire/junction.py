import numpy as np
import scipy.sparse
import mglobals


class Junction:

    __eq__ = lambda self, other: self.start == other.start and self.end == other.end
    __ne__ = lambda self, other: self.start != other.start or self.end != other.end
    __lt__ = lambda self, other: self.start < other.start or (self.start == other.start and self.end < other.end)
    __le__ = lambda self, other: self.start <= other.start or (self.start == other.start and self.end <= other.end)
    __gt__ = lambda self, other: self.start > other.start or (self.start == other.start and self.end > other.end)
    __ge__ = lambda self, other: self.start >= other.start or (self.start == other.start and self.end >= other.end)

    def __init__(self, start, end, donor, acceptor, gene, readN=0, annotated=False):
        ''' The start and end in junctions are the last exon in '''
        
        self.start = start
        self.end = end
        if donor is None:
            self.donor_id = -1
        else:
            self.donor_id = donor.get_id()
        if acceptor is None:
            self.acceptor_id = -1
        else:
            self.acceptor_id = acceptor.get_id()
        self.gene_name = gene.get_id()
        self.txN = 1
        self.annotated = annotated
        self.readN = np.zeros(shape=(mglobals.num_experiments), dtype=np.int)
        self.coverage = scipy.sparse.lil_matrix((mglobals.num_experiments, (mglobals.readLen-16)+1), dtype=np.int)
        self.gc_content = np.zeros(shape=((mglobals.readLen-16)+1), dtype=np.float)
        self.id = "%s:%s-%s" % (self.gene_name, start, end)

    def __hash__(self):
        return hash(self.start) ^ hash(self.end) ^ hash(self.gene_name)

    # GETTERs
    def get_id(self):
        return self.id

    def get_ss_5p(self):
        return self.start

    def get_ss_3p(self):
        return self.end

    def get_count(self):
        return self.txN

    def get_gene(self):
        return mglobals.gene_tlb[self.gene_name]

    def get_coordinates(self):
        return self.start, self.end

    def get_donor(self):
        if self.donor_id == -1:
            ex = None
        else:
            ex = self.get_gene().get_exon_by_id(self.donor_id)
        return ex

    def get_acceptor(self):
        if self.acceptor_id == -1:
            ex = None
        else:
            ex = self.get_gene().get_exon_by_id(self.acceptor_id)
        return ex

    def get_gc_content(self):
        return self.gc_content

    def get_readN(self, idx):
        return self.readN[idx]

    def is_annotated(self):
        return self.annotated

    #MODIFIERs

    def add_gc_content_positions(self, pos, gc):
        self.gc_content[pos] = gc
#        self.gc_factor[exp_idx,:] = gc_factor

    def add_donor(self, donor):
        if donor is None:
            self.donor_id = -1
        else:
            self.donor_id = donor.get_id()

    def add_acceptor(self, acceptor):
        if acceptor is None:
            self.acceptor_id = -1
        else:
            self.acceptor_id = acceptor.get_id()

    def update_junction_read(self, exp_idx, read_n, start, gc, unique):
#        print "J3",self, getrefcount(self)

        self.readN[exp_idx] += read_n
        left_ind = mglobals.readLen - (self.start - start) - 8 + 1
        if unique:
            self.coverage[exp_idx, left_ind] += read_n
        else:
            self.coverage[exp_idx, left_ind] = -1
        self.gc_content[left_ind] = gc

    def add_read_number(self, exp_idx, read_n):
        self.readN[exp_idx] += read_n


class MajiqJunc:

    def __init__(self, jnc, exp_idx):
        self.exons = {}
        if jnc is None:
            self.gc_index = scipy.sparse.lil_matrix((1,(mglobals.readLen-16)+1), dtype=np.int)
            self.name = None
            self.id = "None"

            self.exons['chrom'] = None
            self.exons['coord1'] = [0, 0]
            self.exons['coord2'] = [0, 0]
            self.exons['strand'] = None
            self.coverage = scipy.sparse.lil_matrix((1,(mglobals.readLen-16)+1), dtype=np.int)
        else:
            self.name = "%s:%s-%s" % (jnc.get_gene.get_id(), jnc.get_ss_5p(), jnc.get_ss_3p())
            self.id = "%s:%s-%s" % (jnc.get_gene.get_chromosome(), jnc.get_ss_5p(), jnc.get_ss_3p())

            self.exons['chrom'] = jnc.get_gene.get_chromosome()
            self.exons['strand'] = jnc.get_gene.get_strand()
            if jnc.get_donor() is None:
                self.exons['coord1'] = [0, 0]
            else:
                self.exons['coord1'] = jnc.get_donor().get_coordinates

            if jnc.get_acceptor() is None:
                self.exons['coord2'] = [0, 0]
            else:
                self.exons['coord2'] = jnc.get_acceptor().get_coordinates

            self.coverage = jnc.coverage[exp_idx, :]
            self.gc_factor = scipy.sparse.lil_matrix((1,(mglobals.readLen-16)+1), dtype=np.float)
            for jj in range(mglobals.readLen-16+1):
                dummy = jnc.gc_content[jj]
                self.gc_factor[0, jj] = dummy

    def set_gc_factor(self, exp_idx):
        for jj in xrange(self.gc_factor.shape[0]):
            dummy = self.gc_factor[0, jj]
            if dummy == 0:
                gc_f = 0
            else:
                gc_f = mglobals.gc_factor[exp_idx]( dummy )
            self.gc_factor[0, jj] = gc_f
