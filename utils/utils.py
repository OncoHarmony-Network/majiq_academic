import scipy.io
import numpy as np
#from scipy.sparse import lil_matrix
import scipy.sparse
import pickle,sys
import globals
import random

def __gc_factor_ind(val, exp_idx):
    res = 0
    for ii,jj in enumerate(globals.gc_bins[exp_idx]):
        if val < jj:
            res = ii

    return res


def prepare_MAJIQ_matlab_table(junc_list_per_chr, rand_junc_list_per_chr):

    for chr  in junc_list_per_chr.keys():
        junc_set = junc_list_per_chr[chr]
        rand10k = rand_junc_list_per_chr[chr]
        for name, ind_list in globals.tissue_repl.items() :
            aggr_set = set()
            aggr_rand = set()

            for exp_idx in ind_list :
                jun = set(junc_set[exp_idx])
                print rand10k[exp_idx]
                rand10k[exp_idx].difference(jun)
                list = junc_set[exp_idx]

                aggr_set.union(set(list))
                aggr_rand.union(rand10k[exp_idx])

                mat_file = {}
                mat_file ['experiment'] = globals.exp_list[exp_idx]
                mat_file ['GC_bins'] = globals.gc_bins[exp_idx]
                mat_file ['GC_bins_val'] = globals.gc_bins_val[exp_idx]
                mat_file ['weigh_factor'] = globals.weigh_factor
                mat_file ['Inc'] = {}
                mat_file ['Exc'] = {}
                mat_file ['Inc']['cov'] = np.zeros(shape=(len(list),globals.readLen - 16+1))
                mat_file ['Exc']['cov'] = np.zeros(shape=(len(list),globals.readLen - 16+1))
                mat_file ['Inc']['gc_idx'] = np.zeros(shape=(len(list),globals.readLen - 16+1))
                mat_file ['Exc']['gc_idx'] = np.zeros(shape=(len(list),globals.readLen - 16+1))
                mat_file ['Inc']['gc_val'] = np.zeros(shape=(len(list),globals.readLen - 16+1))
                mat_file ['Exc']['gc_val'] = np.zeros(shape=(len(list),globals.readLen - 16+1))
    
    #            print "SET:", junc_set[exp_idx]
                for ii in range(len(list)) :
                    jinc = list[ii][0]
                    jexc = list[ii][1]
                    gci_inc = np.zeros(shape=(globals.num_experiments,(globals.readLen-16)+1),dtype=np.int)
                    gci_exc = np.zeros(shape=(globals.num_experiments,(globals.readLen-16)+1),dtype=np.int)
                    if not jinc is None:
                        mat_file['Inc']['cov'][ii] = jinc.coverage[exp_idx]
                        gci_inc, gc_inc = jinc.get_gc_factors()
                    if not jexc is None:
                        mat_file['Exc']['cov'][ii] = jexc.coverage[exp_idx]
                        gci_exc, gc_exc = jexc.get_gc_factors()
                    mat_file['Inc']['gc_idx'][ii] = gci_inc[exp_idx,ii]
                    mat_file['Exc']['gc_idx'][ii] = gci_exc[exp_idx,ii]
        #            mat_file['Inc']['gc_val'][ii] = gc_inc[exp_idx,ii]
        #            mat_file['Exc']['gc_val'][ii] = gc_exc[exp_idx,ii]
                    for jj in range(globals.readLen-16+1):
                        print exp_idx, ii, jj
                        dummy = gci_inc[exp_idx,jj]
                        if dummy > 0 :
                            mat_file['Inc']['gc_val'][ii,jj] = globals.gc_bins_val[exp_idx][ dummy -1 ]
                        dummy = gci_exc[exp_idx,jj]
                        if dummy > 0 :
                            mat_file['Exc']['gc_val'][ii,jj] = globals.gc_bins_val[exp_idx][ dummy -1 ]
                if len(rand10k[exp_idx]) < 10000:
                    rand_size = len(rand10k[exp_idx])
                else:
                    rand_size = 10000
                print "RAND", rand_size, "LEN"
                mat_file['rand10k'] = {}
                mat_file['rand10k']['cov'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
                mat_file['rand10k']['gc_idx'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
                mat_file['rand10k']['gc_val'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
                for ii, j in enumerate(random.sample(rand10k[exp_idx], rand_size)) :
                    mat_file['rand10k']['cov'][ii] = j.coverage[exp_idx]
                    for jj in range(globals.readLen-16+1):
                        dummy = __gc_factor_ind(j.get_gc_content()[exp_idx,jj],exp_idx)
                        if dummy>0:
                            mat_file['rand10k']['gc_idx'][ii,jj] = dummy
                            mat_file['rand10k']['gc_val'][ii,jj] = globals.gc_bins_val[exp_idx][dummy-1]

                scipy.io.savemat("./test%s"%globals.exp_list[exp_idx],mat_file,oned_as='row')
            #END for exp_idx

            agreg_file = {}
            agreg_file ['Inc'] = {}
            agreg_file ['Exc'] = {}
            agreg_file ['Inc']['cov'] = np.zeros(shape=(len(aggr_set),globals.readLen - 16+1))
            agreg_file ['Exc']['cov'] = np.zeros(shape=(len(aggr_set),globals.readLen - 16+1))
            for ii,nn in enumerate(aggr_set):
                for exp_idx in ind_list:
                    for jj in range(globals.readLen-16+1):
                        dummy = __gc_factor_ind(j.get_gc_content()[exp_idx,jj],exp_idx)
                        if not nn[0] is None and dummy > 0:
                            agreg_file['Inc']['cov'][ii] += nn[0].coverage[exp_idx] * globals.gc_bins_val[exp_idx][ dummy -1 ]
                        if not nn[1] is None and dummy > 0:
                            agreg_file['Exc']['cov'][ii] += nn[1].coverage[exp_idx]
            if len(aggr_rand) < 10000:
                rand_size = len(aggr_rand)
            else:
                rand_size = 10000
            print "RAND", rand_size, "LEN"
            agreg_file['rand10k'] = {}
            agreg_file['rand10k']['cov'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
            agreg_file['rand10k']['gc_idx'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
            agreg_file['rand10k']['gc_val'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
            for ii, j in enumerate(random.sample(aggr_rand, rand_size)) :
                agreg_file['rand10k']['cov'][ii] = j.coverage[exp_idx]
                for jj in range(globals.readLen-16+1):
                    dummy = __gc_factor_ind(j.get_gc_content()[exp_idx,jj],exp_idx)
                    if dummy>0:
                        agreg_file['rand10k']['gc_idx'][ii,jj] = dummy
                        agreg_file['rand10k']['gc_val'][ii,jj] = globals.gc_bins_val[exp_idx][dummy-1]

            scipy.io.savemat("./agreg_%s"%globals.exp_list[exp_idx],agreg_file,oned_as='row')
    return

def prepare_junctions_gc( junc , exp_idx):
    
    gc = np.zeros(shape=(globals.readLen - 16+1))
    gci = np.zeros(shape=(globals.readLen - 16+1))
    for jj in range(globals.readLen - 16+1) :
        if not junc is None and junc.get_gc_content()[exp_idx,jj] != 0:
            gci[jj] = __gc_factor_ind(junc.get_gc_content()[exp_idx,jj],exp_idx)
            gc[jj] = globals.gc_factor[exp_idx](junc.get_gc_content()[exp_idx,jj])

    if not junc is None:
        junc.add_gc_factor_positions(exp_idx, gci,gc)
    return


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

