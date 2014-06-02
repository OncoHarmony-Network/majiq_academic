import numpy as np
import grimoire.mglobals as mglobals
import scipy

import pdb

SSOURCE = 'source'
STARGET = 'target'

class LSV(object):

    def __init__ ( self, exon, id, junctions, type ):
        if type != SSOURCE and type != STARGET: raise RuntimeError('Incorrect LSV type %s'%type)
        self.coords = exon.get_coordinates()
        self.id = id
        junction_list = [ x for x in  junctions if x is not None] 
        if len(junction_list) < 2 or exon.ir : raise ValueError
        self.type = type
        self.exon = exon

        self.tlb_junc = {}
        self.ext_type = self.set_type(junction_list, self.tlb_junc)
        if self.ext_type == 'intron':
            raise ValueError

        if len(junction_list) > len(self.ext_type.split('|')) -1 :
            print " ERROR_LSV :: with inconsistent junction-type %s, %s"%(len(junction_list), len(self.ext_type.split('|')))

        for kk,vv in self.tlb_junc.items():
            count = np.sum(junction_list[vv].coverage[0].toarray())
            if kk.find('e0') != -1  and count != 0:
                raise ValueError


        self.junctions = []
        order = self.ext_type.split('|')[1:]
        for idx,jj in enumerate(order):
            if jj[-2:] == 'e0': continue
            self.junctions.append(LSV.tlb_junc[jj])

        self.junctions=np.array(self.junctions)

#        if self.check_type(self.ext_type) == -1:
        #    pdb.set_trace()
#            pass


    def check_type(self, type):
        tab = type.split('|')[1:]
        
        exss = []
        targ = {}
        for tt in tab:
            dum = tt.split('e')
            exss.append(int(dum[0]))
            tr = dum[1].split('.')
            if len(tr) == 1 and dum[1]=='0' : continue
            if int(tr[0]) not in targ:
                targ[int(tr[0])] = []
            targ[int(tr[0])].append(int(tr[1]))

        exss.sort()
        for iidx, ii in enumerate(exss[1:]):
            if ii != exss[iidx]+1 and ii != exss[iidx]:
                    print "ERROR 1", type
                    return -1

        for kk, vlist in targ.items():
            vlist.sort()
            for vidx, vv in enumerate(vlist[1:]):
                if vv != vlist[vidx]+1 and vv != vlist[vidx]:
                    print "ERROR 2", type
                    return -1

        return 0

    def get_coordinates(self):
        return self.coords

    def get_junctions_list(self):
        return self.junctions

    def is_Ssource(self):
        return bool(self.type==SSOURCE)

    def is_Starget(self):
        return bool(self.type == STARGET)


    def set_type( self, jlist, tlb_junc ):
        ex_id = self.exon.get_id()
        if self.type == SSOURCE:
            spsite = sorted(set(self.exon.ss_5p_list))
        else:
            spsite = sorted(set(self.exon.ss_3p_list))
        ex_set = set()
        skip = False
        for junc in jlist:
            if self.type == SSOURCE:
                lsv_exon = junc.donor
                if lsv_exon.get_id() != ex_id:
                    skip = True
                    break
#                assert lsv_exon.get_id() == ex_id , "SOURCE, Gene: %s, junc_id %s is different than lsv id %s\n lsv_exon coords %s, ex_idx coords %s, %s"%(junc.get_gene().get_id(), lsv_exon.get_id(), ex_id, lsv_exon.get_coordinates(),jlist[0].donor.get_coordinates(), self.coords)
                
                if not junc.acceptor is None: 
                    ex_set.add( junc.acceptor.get_id() )
            else:
                lsv_exon = junc.acceptor
                if lsv_exon.get_id() != ex_id:
                    skip = True
                    break
#                assert lsv_exon.get_id() == ex_id , "TARGET, Gene: %s, junc_id %s is different than lsv id %s\n lsv_exon coords %s, ex_idx coords %s, %s"%(junc.get_gene().get_id(),junc.acceptor.get_id(), ex_id, lsv_exon.get_coordinates(), jlist[0].acceptor.get_coordinates(),self.coords)
                if not junc.donor is None: 
                    ex_set.add( junc.donor.get_id() )
        if skip : return 'intron'
        ex_list = sorted(list(ex_set))
        ext_type = "%s"%(self.type[0])
        type_set = set()
        for jidx, junc in enumerate(jlist):
            if self.type == SSOURCE:
                if junc.acceptor is None: 
                    exs3 = ''
                    ex = '0'
                else:
                    s3 = sorted(list(set(junc.acceptor.ss_3p_list)))
                    ex1 = ex_list.index(junc.acceptor.get_id())+1
                    ex = '%s.%s'%(ex1,s3.index(junc.end)+1)
                jtype="|%se%s"%(spsite.index(junc.start)+1,ex)
            else:
                if junc.donor is None:
                    exs5 = ''
                    ex = '0'
                else:
                    s5 = sorted(list(set(junc.donor.ss_5p_list)))
                    ex1 = ex_list.index(junc.donor.get_id())+1
                    ex = '%s.%s'%(ex1,s5.index(junc.start)+1)
                jtype = "|%se%s"%(spsite.index(junc.end)+1,ex)
            type_set.add( jtype )
            tlb_junc[jtype[1:]] = jidx
        for tt in sorted(list(type_set)):
            ext_type += tt

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
            gene += '%d\t%d\t'%(lsv_coord[0],jlist[-1].acceptor.get_coordinates()[1])
        else:
            if jlist[0].donor is None: continue
            gene += '%d\t%d\t'%(jlist[0].donor.get_coordinates()[0],lsv_coord[1])

        gene += '.\t%s\t.\tName=%s;Parent=%s;ID=%s'%(strand,lsv.id, lsv.id, lsv.id)
        trans.append(gene)  

        for jidx,junc in enumerate(jlist):

            mrna = '%s\tscript\tmRNA\t'%chrom
            mrna_id = '%s.%d'%(lsv.id,jidx)
            ex1 = '%s\tscript\texon\t'%chrom 
            ex2 = '%s\tscript\texon\t'%chrom

            if lsv.type == SSOURCE:
                if junc.acceptor is None: break
                excoord = junc.acceptor.get_coordinates()
                variant = junc.get_coordinates()
                mrna +='%d\t%d\t'%(lsv_coord[0], excoord[1])
                ex1 += '%d\t%d\t'%(lsv_coord[0], variant[0])
                ex2 += '%d\t%d\t'%(variant[1],excoord[1])
            else:
                if junc.donor is None: break
                excoord = junc.donor.get_coordinates()
                variant = junc.get_coordinates()
                mrna +='%d\t%d\t'%(excoord[0], lsv_coord[1])
                ex1 += '%d\t%d\t'%(variant[1], lsv_coord[1])
                ex2 += '%d\t%d\t'%(excoord[0], variant[0])
            mrna += '.\t%s\t.\tName=%s;Parent=%s;ID=%s'%(strand,mrna_id,lsv.id,mrna_id)
            ex1  += '.\t%s\t.\tName=%s.lsv;Parent=%s;ID=%s.lsv'%(strand, mrna_id, mrna_id,mrna_id)
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
        self.junction_list = scipy.sparse.lil_matrix((len(ind_list),(mglobals.readLen-16)+1),dtype=np.int)
        self.gc_factor = scipy.sparse.lil_matrix( (len(ind_list),(mglobals.readLen-16)+1), dtype=np.dtype('float') )

        for idx,junc in enumerate(LSV.junctions):
            self.junction_list[idx,:] = junc.coverage[exp_idx,:]
            for jidx in range(mglobals.readLen-16+1):
                dummy = junc.get_gc_content()[jidx]
                self.gc_factor[idx,jidx] = dummy


    def set_gc_factor( self , exp_idx):
        for idx in xrange(self.gc_factor.shape[0]):
            for jidx in xrange(self.gc_factor.shape[1]):
                dummy = self.gc_factor[idx,jidx]
                if dummy == 0 :
                    gc_f = 0
                else:
                    gc_f = mglobals.gc_factor[exp_idx]( dummy )
                self.gc_factor[idx,jidx] = gc_f
