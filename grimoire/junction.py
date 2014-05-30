import numpy as np
from utils.decorators import *
import scipy

import mglobals

class Junction:

    __eq__ = lambda self, other: self.start == other.start  and self.end == other.end
    __ne__ = lambda self, other: self.start != other.start  or  self.end != other.end
    __lt__ = lambda self, other: self.start < other.start   or (self.start == other.start and self.end < other.end)
    __le__ = lambda self, other: self.start <= other.start  or (self.start == other.start and self.end <= other.end)
    __gt__ = lambda self, other: self.start > other.start   or (self.start == other.start and self.end > other.end)
    __ge__ = lambda self, other: self.start >= other.start  or (self.start == other.start and self.end >= other.end)

    def __init__ (self, start, end,donor, acceptor,gene,readN=0, annotated=False):
        ''' The start and end in junctions are the last exon in '''
        
        readLength = mglobals.readLen

        self.start = start
        self.end = end
        self.donor = donor
        self.acceptor = acceptor
        self.gene = gene
        self.txN = 1
        self.readN           = np.zeros((mglobals.num_experiments),dtype=np.int)
#        self.gc_content = scipy.sparse.lil_matrix((mglobals.num_experiments,(readLength-16)+1),dtype=np.float)
        self.coverage        = scipy.sparse.lil_matrix((mglobals.num_experiments,(readLength-16)+1),dtype=np.int)
        self.gc_content        = np.zeros( shape=((readLength-16)+1), dtype=np.float )
#        self.gc_factor       = scipy.sparse.lil_matrix((mglobals.num_experiments,(readLength-16)+1),dtype=np.float)
        self.annotated = annotated


    def __hash__(self):
        return hash(self.start) ^ hash(self.end) ^ hash(self.gene.id)

    # GETTERs
    def get_ss_5p(self):
        return self.start
    def get_ss_3p(self):
        return self.end
    def get_count(self):
        return self.txN
    def get_gene (self):
        return self.gene
    def get_coordinates( self ):
        return (self.start,self.end)
    def get_donor (self):
        return (self.donor)
    def get_acceptor (self):
        return (self.acceptor)
    def get_gc_content(self):
        return (self.gc_content)
    def get_readN(self, idx):
        return (self.readN[idx])
#    def get_gc_content(self):
#        return(self.gc_index, self.gc_factor)
    def is_annotated(self):
        return(self.annotated)

    #MODIFIERs


    def add_gc_content_positions(self, pos, gc):
        self.gc_content[pos] = gc
#        self.gc_factor[exp_idx,:] = gc_factor

    def add_donor(self, donor):
        self.donor = donor

    def add_acceptor(self, acceptor):
        self.acceptor = acceptor

    def update_junction_read( self, exp_idx, readN, start,gc,unique ) :
        self.readN[exp_idx] += readN
        left_ind = mglobals.readLen - (self.start - start) - 8 +1
        if unique :
            self.coverage[exp_idx,left_ind]+= readN
        else:
            self.coverage[exp_idx,left_ind]= -1
        self.gc_content[left_ind] = gc

    def add_read_number(self,exp_idx,readN):
        self.readN[exp_idx] += readN


class majiq_junc:

    def __init__ (self, jnc, exp_idx ) :
        if jnc is None:
            self.gc_index = scipy.sparse.lil_matrix((1,(mglobals.readLen-16)+1),dtype=np.int)
            self.name = None
            self.id = "None"
            self.exons = {}
            self.exons['chrom']= None
            self.exons['coord1'] = [0,0]
            self.exons['coord2'] = [0,0]
            self.exons['strand'] = None
            self.coverage = scipy.sparse.lil_matrix((1,(mglobals.readLen-16)+1),dtype=np.int)
        else:
            self.name     = "%s:%s-%s"%(jnc.get_gene().get_id(),jnc.get_ss_5p(),jnc.get_ss_3p())
            self.id     = "%s:%s-%s"%(jnc.get_gene().get_chromosome(),jnc.get_ss_5p(),jnc.get_ss_3p())

            self.exons = {}
            self.exons['chrom']  = jnc.get_gene().get_chromosome()
            self.exons['strand'] = jnc.get_gene().get_strand()
            if jnc.get_donor() is None:
                self.exons['coord1'] = [0,0]
            else:
                self.exons['coord1'] = jnc.get_donor().get_coordinates()

            if jnc.get_acceptor() is None:
                self.exons['coord2'] = [0,0]
            else:
                self.exons['coord2'] = jnc.get_acceptor().get_coordinates()


            self.coverage = jnc.coverage[exp_idx,:]
            #self.coverage = jnc.coverage[exp_idx,:].toarray()
#            self.gc_index = jnc.get_gc_factors()[0][exp_idx,:].toarray()[0]
            self.gc_factor = scipy.sparse.lil_matrix((1,(mglobals.readLen-16)+1),dtype=np.float)
#            self.gc_factor = np.zeros(shape=((mglobals.readLen-16)+1))
            for jj in range(mglobals.readLen-16+1):
                dummy = jnc.gc_content[jj]
                self.gc_factor[0,jj] = dummy


    def set_gc_factor(self, exp_idx):
        for jj in xrange(self.gc_factor.shape[0]):
            dummy = self.gc_factor [0,jj]
            if dummy == 0 :
                gc_f = 0
            else:
                gc_f = mglobals.gc_factor[exp_idx]( dummy )
            self.gc_factor[0,jj] = gc_f
#            self.gc_factor[jj] = mglobals.gc_content[exp_idx](dummy)

#                if dummy > 0 :
#                    self.gc_factor [jj] = mglobals.gc_bins_val[exp_idx][ dummy -1 ]

