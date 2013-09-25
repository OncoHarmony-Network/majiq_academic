import globals
import numpy as np
from utils.decorators import *


class Junction:

    __eq__ = lambda self, other: self.start == other.start  and self.end == other.end
    __ne__ = lambda self, other: self.start != other.start  or  self.end != other.end
    __lt__ = lambda self, other: self.start < other.start   or (self.start == other.start and self.end < other.end)
    __le__ = lambda self, other: self.start <= other.start  or (self.start == other.start and self.end <= other.end)
    __gt__ = lambda self, other: self.start > other.start   or (self.start == other.start and self.end > other.end)
    __ge__ = lambda self, other: self.start >= other.start  or (self.start == other.start and self.end >= other.end)

    def __init__ (self, start, end,donor, acceptor,gene,readN=0):
        ''' The start and end in junctions are the last exon in '''
        
        readLength = globals.readLen

        self.start = start
        self.end = end
        self.donor = donor
        self.acceptor = acceptor
        self.gene = gene
        self.txN = 1
        self.readN = np.zeros(shape=(globals.num_experiments),dtype=np.int)
#        self.readN = 0
        self.gccontent_x_pos = np.zeros(shape=(globals.num_experiments,(readLength-16)+1),dtype=np.float)
        self.coverage = np.zeros(shape=(globals.num_experiments,(readLength-16)+1),dtype=np.int)
        self.out_info = {}

        self.gc_index = np.zeros(shape=(globals.num_experiments,(readLength-16)+1),dtype=np.int)
        self.gc_factor = np.zeros(shape=(globals.num_experiments,(readLength-16)+1),dtype=np.float)
 

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
        return (self.gccontent_x_pos)
    def get_readN(self, idx):
        return (self.readN[idx])
    def get_gc_factors(self):
        return(self.gc_index, self.gc_factor)

    #MODIFIERs

    def add_gc_factor_positions(self, exp_idx, gc_index, gc_factor):
        self.gc_index[exp_idx] = gc_index
        self.gc_factor[exp_idx] = gc_factor

    def add_donor(self, donor):
        self.donor = donor

    def add_acceptor(self, acceptor):
        self.acceptor = acceptor

    def update_junction_read( self, exp_idx, readN, start,gc,unique ) :
        self.readN[exp_idx] += readN
        left_ind = globals.readLen - (self.start - start) - 8 +1
#        print left_ind
        if unique :
            self.coverage[exp_idx,left_ind]+= readN
        else:
            self.coverage[exp_idx,left_ind]= -1
        self.gccontent_x_pos[exp_idx,left_ind] = gc

    def add_read_number(self,exp_idx,readN):
        self.readN[exp_idx] += readN

    @deprecated
    def out_junc_information(self,readLength=76):
        ltemp = (readLength-8)
        set_OL = False
        set_OR = False

        self.out_info['pos']=np.zeros(shape=(1,2*ltemp),dtype=np.int)
        self.out_info['val']= self.coverage
        for ii in range(ltemp):
            self.out_info['pos'][0,ii] = self.start-ltemp+ii
            inv_index = (2*ltemp)-1 -ii
            self.out_info['pos'][0,inv_index] = self.end + ltemp - ii
            #set OL
            if not set_OL  and self.out_info['val'][0,ii] >0 :
                set_OL = True
                self.out_info['OL']= ii
            #set OR
            if not set_OR  and self.out_info['val'][0, inv_index] >0 :
                set_OR = True
                self.out_info['OR']= inv_index
#        print self.coverage
        self.out_info['GC']= 30
       
        return self.out_info


