import numpy as np
import mglobals


SSOURCE = 'source'
STARGET = 'target'

class LSV(object):

    def __init__ ( self, coords, id, junctions, type, exp_idx ):
        if type != SSOURCE and type != STARGET: raise RuntimeError('Incorrect LSV type %s'%type)
        self.coords = coords
        self.id = id
        #print "INIT",junctions
        junction_list = [ x for x in  junctions if x is not None] 
        if len(junction_list) == 0 : return
        self.junctions = np.asarray(junction_list)
        self.type = type
        self.ext_type = self.set_type(junction_list, type)
#        self.exp_idx = exp_idx
        print "LSV :::::::::::::::::", self.ext_type, self.id

    def get_coordinates(self):
        return self.coords

    def get_junctions_list(self):
        return self.junction_list

    def is_Ssource(self):
        return bool(self.type==SSOURCE)

    def is_Starget(self):
        return bool(self.type == STARGET)


    def set_type( self, jlist, type ):
        if type == SSOURCE:
            ex_id = jlist[0].donor.get_id()
            spsite = sorted(list(set(jlist[0].donor.ss_5p_list)))
        else:
            ex_id = jlist[0].acceptor.get_id()
            spsite = sorted(list(set(jlist[0].acceptor.ss_3p_list)))
        ex_set = set()
        intron_ret = False
        for junc in jlist:
            if type == SSOURCE:
                #print "MMM",ex_id, junc.donor.get_id()
                if junc.donor.get_id() != ex_id : 
                #    print "ERRORRR set_type SS:", ex_id, junc.donor.get_id(), junc
                    intron_ret = True
                    break
                if junc.acceptor is None: continue
                ex_set.add(junc.acceptor.get_id())
            else:
                #print "MMM",ex_id, junc.acceptor.get_id()
                if junc.acceptor.get_id() != ex_id : 
                #    print "ERRORRR set_type ST:", ex_id, junc.acceptor.get_id(), junc
                    intron_ret = True
                    break
                if junc.donor is None: continue
                ex_set.add(junc.donor.get_id())
        if intron_ret : return "intron ret"
        #print type, spsite
        ex_list = sorted(list(ex_set))
        ext_type = "%s"%(self.type[0])
        type_set = set()
        for junc in jlist:
            if type == SSOURCE:
                #print junc.start, junc.acceptor
                if junc.acceptor is None: continue
                s3 = sorted(list(set(junc.acceptor.ss_3p_list)))
                type_set.add("|%se%s.%s"%(spsite.index(junc.start)+1,ex_list.index(junc.acceptor.get_id())+1,s3.index(junc.end)+1))
            else:
                if junc.donor is None: continue
                s5 = sorted(list(set(junc.donor.ss_5p_list)))
                #print junc.end, spsite
                #print junc.donor.get_id(), ex_list
                #print junc.start, s5
                type_set.add("|%se%s.%s"%(spsite.index(junc.end)+1,ex_list.index(junc.donor.get_id())+1,s5.index(junc.start)+1))
        
        for tt in sorted(list(type_set)):
            ext_type += tt

        print "LSV::::::::::::::::::", ext_type 
        return ext_type

    def is_equivalent (self, variant):
        if self.type == variant.type : return False
        return np.array_equal(self.junctions, variant.junctions)

    def to_majiqLSV (self):
        return Majiq_LSV(self)

def print_lsv_extype ( list_lsv, filename ):
    fp = open(filename,'w+')
    print list_lsv.shape
    for idx in range(list_lsv.shape[0]):
        lsv = list_lsv[idx]
        fp.write("%s\n"%(lsv.type))
    fp.close()

class LSV_IR ( object):

    def __init__ (self, start, end, exon_list, gene):
        self.start = start
        self.end = end
        self.exonlist = exon_list
        self.gene = gene
        gene.add_intron_retention(self)




class Majiq_LSV(object):

    def __init__ (self, LSV):

        self.coords = LSV.coords
        self.id = LSV.id
        self.type = LSV.ext_type
        exp_idx = LSV.exp_idx
        self.junction_list = np.zeros( shape=(len(LSV.junctions)),dtype=np.dtype('object'))
        for idx,jj in enumerate(LSV.junctions):
            self.junction_list[idx] = jj.coverage[exp_idx,:]
            self.gc_index  = jj.get_gc_factors()[0][exp_idx,:].toarray()[0]
            self.gc_factor = np.zeros(shape=((mglobals.readLen-16)+1))
            for jidx in range(mglobals.readLen-16+1):
                dummy = self.gc_index [jidx]
                if dummy > 0 :
                    self.gc_factor [jidx] = mglobals.gc_bins_val[exp_idx][ dummy -1 ]

