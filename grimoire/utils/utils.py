import scipy.io
import numpy as np
import scipy.sparse
import pickle,sys
import mglobals
import random

def __gc_factor_ind(val, exp_idx):
    res = 0
    for ii,jj in enumerate(mglobals.gc_bins[exp_idx]):
        if val < jj:
            res = ii

    return res

def prepare_MAJIQ_table(junc_set, non_as, temp_file):

    print mglobals.tissue_repl
    for name, ind_list in mglobals.tissue_repl.items() :

        for idx,exp_idx in enumerate(ind_list) :
            info = {}
            info ['weigh_factor'] = mglobals.weigh_factor
            info ['experiment']   = mglobals.exp_list[exp_idx]
            info ['GC_bins']      = mglobals.gc_bins[exp_idx]
            info ['GC_bins_val']  = mglobals.gc_bins_val[exp_idx]

            jun = set(junc_set[exp_idx])
            non_as[exp_idx].difference(jun)
            jlist = junc_set[exp_idx]
            majiq_table_as    = np.zeros( shape=(len(j_list),2), dtype=np.dtype('object'))
            majiq_table_nonas = np.zeros( shape=(len(non_as),1), dtype=np.dtype('object'))

            # We iterate over the inc and exc in order to fill the majiq_junc_matrix
            for ii in range(len(junc_list)) :
                for lab_idx in range(2):
                    majiq_table_as = majiq_junc( junc_set[ii][lab_idx], exp_idx)
            for jn in non_as :
                    majiq_table_as = majiq_junc( jn , exp_idx)

            file_pi = open(temp_file, 'w+') 
            pickle.dump((info,majiq_table_as, majiq_table_nonas), file_pi)
            file_pi.close()


def prepare_junctions_gc( junc , exp_idx):
    
    gc = np.zeros(shape=(mglobals.readLen - 16+1))
    gci = np.zeros(shape=(mglobals.readLen - 16+1))
    for jj in range(mglobals.readLen - 16+1) :
        if not junc is None and junc.get_gc_content()[exp_idx,jj] != 0:
            gci[jj] = __gc_factor_ind(junc.get_gc_content()[exp_idx,jj],exp_idx)
            gc[jj] = mglobals.gc_factor[exp_idx](junc.get_gc_content()[exp_idx,jj])

    if not junc is None:
        junc.add_gc_factor_positions(exp_idx, gci,gc)
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


