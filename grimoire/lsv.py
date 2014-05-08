import numpy as np
import mglobals


SSOURCE = 'source'
STARGET = 'target'

class LSV(object):

    def __init__ ( self, coords, id, junctions, type ):
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

    def to_majiqLSV (self, exp_idx):
        return Majiq_LSV(self, exp_idx)

def extract_SE_events( list_lsv_per_gene ):

    sslist = list_lsv_per_gene[0]
    stlist = list_lsv_per_gene[1]
    
    for ss in sslist:
        slist = ss.junction_list
        if len(slist) != 2: continue
        if slist[0].acceptor == slist[1].acceptor: continue
        
        sindx = None
        tindx = None
        C2 = None

        for st in stlist:
            tlist = st.junction_list
            if len(tlist) != 2: continue
            if tlist[0].donor == tlist[1].donor: continue

            for ii in range(2):
                for jj in range(2):
                    if slist[ii] == tlist[jj]:
                        sindx = 1-ii
                        tindx = 1-jj
                        break
                else:
                    continue
                break
            else:
                continue
            if slist[sindx].acceptor == tlist[tindx].donor: 
                C2 = st.coords
                break
        else:
            continue
        
        C1 = ss.coords
        A = slist[sindx].acceptor.get_coordinates()
        #ret_list.append( (C1.)



def lsv_to_gff( list_lsv ):


    gtf = []
    for lsv in list_lsv[0]:
        trans = []
        jlist = sorted(lsv.junctions)
        lsv_coord = lsv.get_coordinates()

        gne = jlist[0].gene
        chrom = gne.get_chromosome()
        strand = gne.get_strand()
        gene = '%s\tscript\tgene\t'%chrom
        if lsv.type==SSOURCE:
            if jlist[-1].acceptor is None: continue
            gene += '%d\t%d\t'%(lsv.get_coordinates()[0],jlist[-1].acceptor.get_coordinates()[1])
        else:
            if jlist[0].donor is None: continue
            gene += '%d\t%d\t'%(jlist[0].donor.get_coordinates()[0],lsv.get_coordinates()[1])

        gene += '.\t%s\t.\tName=%s;Parent=%s;ID=%s'%(strand,lsv.id, lsv.id, lsv.id)
        trans.append(gene)  

        for jidx,junc in enumerate(jlist):
            mrna = '%s\tscript\tmRNA\t'%chrom
            mrna_id = '%s.%d'%(lsv.id,jidx)
            ex1 = '%s\tscript\texon\t%d\t%d\t.\t%s\t.\tName=%s.lsv;Parent=%s;ID=%s.lsv'%(chrom, lsv_coord[0], lsv_coord[1], strand, mrna_id, mrna_id,mrna_id)
            ex2 = '%s\tscript\texon\t'%chrom

            if lsv.type == SSOURCE:
                if junc.acceptor is None: break
                excoord = junc.acceptor.get_coordinates()
                mrna +='%d\t%d\t'%(lsv_coord[0],excoord[1])
                ex2 += '%d\t%d\t'%(excoord[0],excoord[1])
            else:
                if junc.donor is None: break
                excoord = junc.acceptor.get_coordinates()
                mrna +='%d\t%d\t'%(junc.donor.get_coordinates()[0],lsv.get_coordinates()[1])
                ex2 += '%d\t%d\t'%(excoord[0],excoord[1])
            mrna += '.\t%s\t.\tName=%s;Parent=%s;ID=%s'%(strand,mrna_id,lsv.id,mrna_id)
            ex2  += '.\t%s\t.\tName=%s.ex;Parent=%s;ID=%s.ex'%(strand, mrna_id, mrna_id, mrna_id)
        
            trans.append(mrna)  
            trans.append(ex1)
            trans.append(ex2)
        else:
            gtf.extend(trans)

    return gtf


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

    def __init__ (self, LSV, exp_idx):

        self.coords = LSV.coords
        self.id = LSV.id
        self.type = LSV.ext_type
#        exp_idx = LSV.exp_idx
        self.junction_list = np.zeros( shape=(len(LSV.junctions)),dtype=np.dtype('object'))
        for idx,jj in enumerate(LSV.junctions):
            self.junction_list[idx] = jj.coverage[exp_idx,:]
            self.gc_index  = jj.get_gc_factors()[0][exp_idx,:].toarray()[0]
            self.gc_factor = np.zeros(shape=((mglobals.readLen-16)+1))
            for jidx in range(mglobals.readLen-16+1):
                dummy = self.gc_index [jidx]
                if dummy > 0 :
                    self.gc_factor [jidx] = mglobals.gc_bins_val[exp_idx][ dummy -1 ]

