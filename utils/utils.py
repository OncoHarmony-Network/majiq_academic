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



def _create_mat_dict(exp_idx, junc_list,rand_list):

    mat_file = {}
    mat_file ['experiment'] = globals.exp_list[exp_idx]
    mat_file ['GC_bins'] = globals.gc_bins[exp_idx]
    mat_file ['GC_bins_val'] = globals.gc_bins_val[exp_idx]
    mat_file ['weigh_factor'] = globals.weigh_factor

    for lab_idx, label in enumerate(('Inc','Exc')):
        mat_file [label] = {}
        mat_file [label]['cov']    = np.zeros(shape=(len(junc_list),globals.readLen - 16+1))
        mat_file [label]['gc_idx'] = np.zeros(shape=(len(junc_list),globals.readLen - 16+1))
        mat_file [label]['gc_val'] = np.zeros(shape=(len(junc_list),globals.readLen - 16+1))

        for ii in range(len(junc_list)) :
            jnc = junc_list[ii][lab_idx]
            gci = np.zeros(shape=((globals.readLen-16)+1),dtype=np.int)
            if not jnc is None:
                mat_file[label]['cov'][ii] = jnc.coverage[exp_idx,:].toarray()
                #gci_inc, gc_inc = jinc.get_gc_factors()
                gci = jnc.get_gc_factors()[0][exp_idx,:].toarray()[0]
            mat_file[label]['gc_idx'][ii] = gci
            #mat_file['Inc']['gc_val'][ii] = gc_inc[exp_idx,ii]
            for jj in range(globals.readLen-16+1):
                dummy = gci[jj]
                if dummy > 0 :
                    mat_file[label]['gc_val'][ii,jj] = globals.gc_bins_val[exp_idx][ dummy -1 ]

    rand_size = len(rand_list)
    mat_file['rand10k'] = {}
    mat_file['rand10k']['cov']    = np.zeros(shape=(rand_size,globals.readLen - 16+1))
    mat_file['rand10k']['gc_idx'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
    mat_file['rand10k']['gc_val'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
    for ii, j in enumerate(rand_list):
        mat_file['rand10k']['cov'][ii] = j.coverage[exp_idx,:].toarray()
        for jj in range(globals.readLen-16+1):
            dummy = __gc_factor_ind(j.get_gc_content()[exp_idx,jj],exp_idx)
            if dummy>0:
                mat_file['rand10k']['gc_idx'][ii,jj] = dummy
                mat_file['rand10k']['gc_val'][ii,jj] = globals.gc_bins_val[exp_idx][dummy-1]

    return mat_file


def prepare_MAJIQ_matlab_table(junc_set, rand10k):

    result_per_tiss = {}
    result_per_exp = {}
    result_p2p  = {}

    for name, ind_list in globals.tissue_repl.items() :
        aggr_set = set()
        aggr_rand = set()

        for idx,exp_idx in enumerate(ind_list) :
            jun = set(junc_set[exp_idx])
            rand10k[exp_idx].difference(jun)
            jlist = junc_set[exp_idx]

            aggr_set.union(jun)
            aggr_rand.union(rand10k[exp_idx])

            if len(rand10k[exp_idx]) < 10000: rand_size = len(rand10k[exp_idx])
            else: rand_size = 10000
            rand_list = random.sample(rand10k[exp_idx], rand_size)

            result_per_exp[exp_idx]= _create_mat_dict(exp_idx,jlist,rand_list)
#            scipy.io.savemat("./test%s"%globals.exp_list[exp_idx],mat_file,oned_as='row')
            
            #paired
            rand_pair = set(rand10k[exp_idx])
            jidx = idx +1
            while jidx < len(ind_list) :
                exp_idx2 = ind_list[jidx]
                rand_pair.intersection(set(rand10k[exp_idx2]))
                jun.intersection(set(junc_set[exp_idx2]))

                if len(rand10k[exp_idx]) < 10000: rand_size = len(rand10k[exp_idx])
                else: rand_size = 10000
                rand_pair_list = random.sample(rand_pair, rand_size)

                mat_file = {}
                mat_file[globals.exp_list[exp_idx]] = _create_mat_dict(exp_idx,list(jun), rand_pair_list)
                mat_file[globals.exp_list[exp_idx2]] = _create_mat_dict(exp_idx2,list(jun), rand_pair_list)
                p2p_id = "%s_%s"%(globals.exp_list[exp_idx],globals.exp_list[exp_idx2])
                result_p2p[p2p_id] = mat_file
                jidx += 1
        #END for exp_idx
        
            

        agreg_file = {}
        for lab_idx, label in enumerate(('Inc','Exc')):
            agreg_file [label] = {}
            agreg_file [label]['cov'] = np.zeros(shape=(len(aggr_set),globals.readLen - 16+1), dtype = np.dtype("float"))
            for ii,nn in enumerate(aggr_set):
                for exp_idx in ind_list:
                    for jj in range(globals.readLen-16+1):
                        dummy = __gc_factor_ind(j.get_gc_content()[exp_idx,jj],exp_idx)
                        if not nn[lab_idx] is None and dummy > 0:
                            agreg_file[label]['cov'][ii] += nn[lab_idx].coverage[exp_idx,jj]*  dummy 

        if len(aggr_rand) < 10000:
            rand_size = len(aggr_rand)
        else:
            rand_size = 10000
        print "RAND", rand_size, "LEN"
        agreg_file['rand10k'] = {}
        agreg_file['rand10k']['cov'] = np.zeros(shape=(rand_size,globals.readLen - 16+1))
        for ii, j in enumerate(random.sample(aggr_rand, rand_size)) :
            for jj in range(globals.readLen-16+1):
                dummy = __gc_factor_ind(j.get_gc_content()[exp_idx,jj],exp_idx)
                agreg_file['rand10k']['cov'][ii] += j.coverage[exp_idx,jj]*dummy

        result_per_tiss[name]=agreg_file
        #scipy.io.savemat("./agreg_%s"%globals.exp_list[exp_idx],agreg_file,oned_as='row')
    return (result_per_tiss, result_per_exp,result_p2p)


def merge_and_create_MAJIQ_matlab(res_x_tiss,res_x_exp, res_p2p):

    header = ['Inc','Exc','rand10k']

    label = ['cov','gc_idx','gc_val']

    for name, ind_list in globals.tissue_repl.items():
        print name
        agreg_file = {}
        for hh in header:
            agreg_file[hh] = {}
            agreg_file[hh]['cov'] = []
            for chrom in res_x_tiss.keys():
                for row in res_x_tiss[chrom][name][hh]['cov']: 
                    agreg_file[hh]['cov'].append(row)
        scipy.io.savemat("./agreg_%s"%name,agreg_file,oned_as='row')

        for idx, exp_idx in enumerate(ind_list) :
            mat_file = {}
            for hh in header:
                mat_file[hh] = {}
                for ll in label:
                    mat_file[hh][ll] = []
                    for chrom in res_x_exp.keys():
                        for row in res_x_exp[chrom][exp_idx][hh][ll]: mat_file[hh][ll].append(row)
            scipy.io.savemat("./%s"%globals.exp_list[exp_idx],mat_file,oned_as='row')

            jidx = idx +1

            while jidx < len(ind_list) :
                exp_idx2 = ind_list[jidx]
                p2p_id = "%s_%s"%(globals.exp_list[exp_idx],globals.exp_list[exp_idx2])
                mat_file = {}
                for idx_e in (exp_idx,exp_idx2):
                    exp = globals.exp_list[idx_e]
                    mat_file[exp] = {}
                    for hh in header:
                        mat_file[exp][hh] = {}
                        for ll in label:
                            mat_file[exp][hh][ll] = []
                            for chrom in res_p2p.keys():
#                                p2p_temp = res_p2p[chrom][p2p_id][exp][hh][ll]
                                for row in res_p2p[chrom][p2p_id][exp][hh][ll]: mat_file[exp][hh][ll].append(row)

                scipy.io.savemat("./%s_%s"%(globals.exp_list[exp_idx],globals.exp_list[exp_idx2]),mat_file,oned_as='row')
                jidx +=1

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


