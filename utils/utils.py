import scipy.io
import numpy as np
#from scipy.sparse import lil_matrix
import scipy.sparse
import pickle,sys


'''
    @create_junction_matlab_file
    @param list_events: this is a list of a triplet exon (ex1,ex2,ex3). 
        Each exon is a Exon object correctly identified.
    @param fname: The name of the matlab output file
    @return: None
    @io: Create a matlab file following the name on fname param
'''
def create_junction_matlab_file(gene_list, fname,chr_list,order):

    n_genes = 0
    for chr,gl in gene_list.items():
        n_genes += len(gl)

    notskip =0
    skipped = 0
    skipped_gene = 0
    gidx = 0
    genes = np.zeros( shape=n_genes, dtype=np.dtype('object'))
    for chr in chr_list:
        for gk in order[chr]:
            gn = gene_list[chr][gk]
#        for gn in gene_list[chr].values():
#            print gn
    
            dict_ss = gn.get_all_ss()
            header = []
            juncs = []
            tlb = {}
            head_idx = 0
            for ex_id, ss_lst in dict_ss.items():
                idx_3ss = 0
                idx_5ss = 0
                for ss,type in ss_lst:
                    if type == '3prime':
                        headline = "%s:exon%s"%(idx_3ss,ex_id)
                        idx_3ss+=1
                    else:
                        headline = "exon%s:%s"%(ex_id,idx_5ss)
                        idx_5ss+=1
                    juncs.append(ss)
                    header.append(headline)
                    tlb[int(ss)] = head_idx
                    head_idx +=1

            num_reads = len(gn.get_RNAread_list())
            num_ss = len(header)

#            if num_reads == 0 or num_ss==0 : continue
            if num_ss == 0:
                header = ['None']
                juncs = ['None']
                num_ss = 1
            genes[gidx] = {}
            genes[gidx]['id'] = gn.get_id()
            genes[gidx]['strand'] = gn.get_strand()
            genes[gidx]['chrom']  = gn.get_chromosome()
            genes[gidx]['header']= np.asarray(header,dtype=np.object)
            genes[gidx]['RPKM']= gn.get_RPKM()
            genes[gidx]['total_reads']= gn.get_read_count()
            genes[gidx]['reads']= np.zeros(shape=(num_reads,num_ss),dtype=np.dtype('int'))
            if num_reads == 0 :
                skipped_gene +=1
                sparse_mat = []
            else:
                sparse_mat = np.zeros(shape=(num_reads,num_ss),dtype=np.dtype('int'))

            genes[gidx]['GC']  = np.zeros(shape=(num_reads),dtype=np.dtype('float'))
            genes[gidx]['coordinates'] = np.zeros(shape=(num_reads,2),dtype=np.dtype('int'))
            genes[gidx]['junc'] = np.asarray(juncs,dtype=np.dtype('int'))
            r_idx = 0
            found = True
            for read in gn.get_RNAread_list():
                coord = read.get_coordinates()
                for jss5,jss3 in  read.get_junctions_info():
                    if not jss5 in tlb or not jss3 in tlb :
                        skipped +=1
                        found = False
                        continue
                    
                    notskip +=1
                    ss_idx = tlb[jss5]
                    if read.unique :
                        val = read.get_read_count()
                    else:
                        val = -1*read.get_read_count()
                    sparse_mat[r_idx,ss_idx] = val
#                    genes[gidx]['reads'][r_idx,ss_idx] = val
                    ss_idx = tlb[jss3]
                    sparse_mat[r_idx,ss_idx] = val
#                    genes[gidx]['reads'][r_idx,ss_idx] = val
                if not found : break

                genes[gidx]['GC'][r_idx] = read.get_GCcontent()
                genes[gidx]['coordinates'][r_idx,0] = coord[0]
                genes[gidx]['coordinates'][r_idx,1] = coord[1]
                r_idx += 1
            if r_idx > 0:
                genes[gidx]['reads']= scipy.sparse.coo_matrix(sparse_mat[:r_idx])
            else:
                genes[gidx]['reads']= []

            gidx +=1
    print "SKIPPED READS :",skipped
    print "NOT SKIPPED READS :",notskip
    print "SKIPPED GENES :",skipped_gene

#    with open(fname+'.dat', 'wb') as outfile:
#        pickle.dump(genes, outfile, pickle.HIGHEST_PROTOCOL)

    scipy.io.savemat(fname,{'RNASEQ_INFO':genes},oned_as='row')
    return

def print_junc_matrices(mat, tlb=None,fp=None):
    if fp is None:
        out = open ('./junc_matrix.tab', 'a+')
    else:
        out = sys.stdout
    out.write("\n=== BEGIN %s === \n\n"%id)

    (N,M)= mat.shape
    header = [0]*N
    if not tlb is None:
        out.write("Nan\t")
        for ex,(p1,p2) in tlb.items():
            for nid,n in enumerate(p1):
                out.write("%d\t"%(ex+1))
            for nid,n in enumerate(p2):
#                header[nid] = "%d:%d\t"%(ex,nid)
                header[n] = "%d"%(ex+1)
    out.write("\n") 
    for ii in np.arange(N):
        out.write("%s\t"%header[ii])
        for jj in np.arange(M):
            val = mat[ii,jj]
            out.write("%s\t"%val)
        out.write("\n")
    out.write("\n=== END %s === \n\n"%id)
    if fp is None: out.close()


'''
 @deprecated 
'''
def prepare_const_set(wei,T,M,user_lab=''):
    '''
    (W) Weijun set
    (T) Transcript const set
    (M) MGP ReadSeq const set
    '''

    TM = set(T).intersection(set(M))
    fp1 = open("./new","w+")
    fp2 = open("./new.coord","w+")
    for ii in TM:
        if ii.gene.strand == '+':
            invj = ii.gene.exonNum - (ii.id-1)
            j = ii.id
        else:
            j = ii.gene.exonNum - (ii.id -1)
            invj = ii.id
        fp1.write("%s.%s:%s:%s:%s\n"%(ii.gene.id,j,str(j-1),invj,ii.gene.exonNum))
        ii.print_triplet_coord(fp2)
        
    fp1.close()
    fp2.close()
    fp1 = open("./wei","w+")
    for ii in wei:
        if ii[1].gene.strand == '+':
            invj = ii[1].gene.exonNum - (ii[1].id-1)
            j = ii[1].id
        else:
            j = ii[1].gene.exonNum  - (ii[1].id-1)
            invj = ii[1].id
        fp1.write("%s.%s:%s:%s:%s\n"%(ii[1].gene.id,j,str(j-1),invj,ii[1].gene.exonNum))
    fp1.close()

    return


'''
@deprecated
'''
def store_binaries(foldername, all_gene_by_chrom):

    if not os.path.exists(foldername):
            os.makedirs(foldername)
    for chr, list in all_gene_by_chrom.items():
        file_pi = open("%s/%s.data"%(foldername,chr), 'wb')
        pickle.dump(list, file_pi)
        file_pi.close()

'''
@deprecated
'''
def load_binaries(foldername, all_gene_by_chrom):

    if not os.path.exists(foldername):
        print "incorrect bin folder"
        return None

    data = {}
    for ff in  os.listdir(foldername) :
        inp = open("%s/%s"%(foldername,ff),'rb')
        tab = ff.split('.')
        if tab[1] == 'data' and tab[0].starswith('chr'):
            chr = tab[0]
        else:
            continue
        pk = pickle.Unpickler(inp)
        data[chr] = pk.load()
        inp.close()

    return data

