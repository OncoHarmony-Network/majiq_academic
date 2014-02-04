import matplotlib
matplotlib.use('Agg')
import os
import sys
import random
import logging
import scipy.io
import numpy as np
import pickle,sys
import grimoire.mglobals as mglobals

from grimoire.junction import majiq_junc
from itertools import izip
import scipy.sparse
import scipy.io
from scipy.stats.mstats import mquantiles
from scipy import interpolate
from matplotlib import pyplot

def create_if_not_exists(my_dir, logger=False):
    "Create a directory path if it does not exist"
    if not os.path.exists(my_dir):
        if logger: logger.info("\nCreating directory %s..."%my_dir)
        os.makedirs(my_dir)   


def get_logger(logger_name, silent=False, debug=False): 
    """
    Returns a logger instance. verbose = False will silence the logger, debug will give 
    more information intended for debugging purposes.
    """
    logging_format= "%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s"
    logging.basicConfig(filename=logger_name, format=logging_format)
    logger = logging.getLogger(logger_name)
    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    if debug: 
        ch.setLevel(logging.DEBUG)
    elif not silent:
        ch.setLevel(logging.INFO)
    else:
        ch.setLevel(logging.WARNING)

    formatter = logging.Formatter("%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger



def __gc_factor_ind(val, exp_idx):
    res = 0
    for ii,jj in enumerate(mglobals.gc_bins[exp_idx]):
        if val < jj:
            res = ii

    return res

def prepare_MAJIQ_table(junc_set, non_as, temp_file):

    print "PREPAAAER",  mglobals.tissue_repl
    for name, ind_list in mglobals.tissue_repl.items() :

        for idx,exp_idx in enumerate(ind_list) :
#            info = {}
#            info ['weigh_factor'] = mglobals.weigh_factor
#            info ['experiment']   = mglobals.exp_list[exp_idx]
#            info ['GC_bins']      = mglobals.gc_bins[exp_idx]
#            info ['GC_bins_val']  = mglobals.gc_bins_val[exp_idx]

            jun = set(junc_set[exp_idx])
            non_as[exp_idx].difference(jun)
            majiq_table_as    = np.zeros( shape=(len(junc_set[exp_idx]),2), dtype=np.dtype('object'))
            majiq_table_nonas = np.zeros( shape=(len(non_as[exp_idx]),1), dtype=np.dtype('object'))

            # We iterate over the inc and exc in order to fill the majiq_junc_matrix
            
            for iix, jn_lst in enumerate(junc_set[exp_idx]) :
                for lab_idx in range(2):
                    #print jn_lst
                    majiq_table_as[iix, lab_idx] = majiq_junc( jn_lst[lab_idx], exp_idx)
            for jix, jn in enumerate(non_as[exp_idx]) :
                    majiq_table_nonas[jix] = majiq_junc( jn , exp_idx)

            file_pi = open("%s/temp_%s.%s"%(mglobals.temp_oDir[exp_idx],mglobals.exp_list[exp_idx], temp_file), 'w+')
            pickle.dump((majiq_table_as, majiq_table_nonas), file_pi)
            file_pi.close()


def merge_and_create_MAJIQ ( chr_list, ofile ):
    for name, ind_list in mglobals.tissue_repl.items() :
        for idx,exp_idx in enumerate(ind_list) :
            as_table = []
            nonas_table = []
            info = {}
#            info ['weigh_factor'] = mglobals.weigh_factor
            info ['experiment']   = mglobals.exp_list[exp_idx]
            info ['GC_bins']      = mglobals.gc_bins[exp_idx]
            info ['GC_bins_val']  = mglobals.gc_bins_val[exp_idx]
            for chrom in chr_list:
                filename = '%s/temp_%s.%s.obj'%(mglobals.temp_oDir[exp_idx], mglobals.exp_list[exp_idx], chrom)
                if not os.path.exists(filename): continue
                file_pi2 = open(filename, 'rb')
                as_t,non_as = pickle.load(file_pi2)
                as_table.append(as_t)
                nonas_table.append(non_as)
            if len(as_table)==0: continue
            AT  = np.concatenate((as_table))
            NAT = np.concatenate((nonas_table))
            file_pi = open('%s/toJuan.%s.majiq'%(mglobals.outDir,mglobals.exp_list[exp_idx]), 'w+')
            pickle.dump((info,AT, NAT), file_pi)
            file_pi.close()



def set_exons_gc_content(chrom, exon_list ):

    out_list = {}
    current_chrom = None
    loaded_chrom = ''
    fastadir_path = "%s/Genomes/goldenPath/%s/"%(os.environ["ASP_DATA_ROOT"],mglobals.genome)

    print "Loading chromosome... %s"%chrom
    chrom_path = fastadir_path + chrom + ".fa"
    if not os.path.exists(chrom_path): return 
    chrom_file = open(chrom_path)
    loaded_chrom = []
    for chrom_line in chrom_file:
        if not chrom_line.startswith(">"):
            loaded_chrom.append(chrom_line.strip("\n"))
    loaded_chrom  = ''.join(loaded_chrom)

    for exon in exon_list:
        strt,end = exon.get_coordinates()
        sequence = loaded_chrom[strt:end]
        #reverse the sequence if the strand is reverse
        sequence = sequence.lower()
        if exon.get_strand() == "-":
            new_seq = []
            for char in sequence[::-1]:
                if char == 'g':
                    new_seq.append('c')
                elif char == 'c':
                    new_seq.append('g')
                elif char == 'a':
                    new_seq.append('t')
                elif char == 't':
                    new_seq.append('a')
                else:
                    new_seq.append(char)
            sequence = ''.join(new_seq)
        if len(sequence) == 0 : 
            print "KKKKseq",exon.get_coordinates()
        exon.set_gc_content(sequence)


def prepare_junctions_gc( junc , exp_idx):

    gc = np.zeros(shape=(mglobals.readLen - 16+1))
    gci = np.zeros(shape=(mglobals.readLen - 16+1))
    for jj in range(mglobals.readLen - 16+1) :
        if not junc is None and junc.get_gc_content()[exp_idx,jj] != 0:
            #gci[jj] = __gc_factor_ind(junc.get_gc_content()[exp_idx,jj],exp_idx)
            pass
            #gc[jj] = mglobals.gc_factor[exp_idx](junc.get_gc_content()[exp_idx,jj])

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


def gc_factor_calculation(exon_list, nb):

    local_bins     = np.zeros( shape=(mglobals.num_experiments,nb+1),dtype=np.dtype('float'))
    local_meanbins = np.zeros( shape=(mglobals.num_experiments,nb),dtype=np.dtype('float'))
    local_factor   = np.zeros( shape=(mglobals.num_experiments,nb),dtype=np.dtype('float'))
    
    dummy_counter = 0

    print mglobals.tissue_repl
    for tissue, list_idx in mglobals.tissue_repl.items():
        for exp_n in list_idx :
            print "EXP", exp_n
            count = []
            gc = []
            for idx, ex in enumerate(exon_list):
                gc_val = ex.get_gc_content()
                st,end = ex.get_coordinates()
                cov = ex.get_coverage(exp_n)

                # TEST AND CHECK
#                if gc_val is None or cov == 0:
#                    print ex.strand, st, end


                if  gc_val is None or end-st < 30  or cov < 5: continue
                count.append( cov )
                gc.append( gc_val )
            if len(gc) == 0 : continue
            print "cont", len(count)
            print count
            print gc
            count,gc = izip(*sorted(izip(count, gc), key=lambda x: x[1]))
            print count
            print gc

            num_regions = len(count)
            nperbin =  num_regions / nb

            quant_median =[0.0]*8
            mean_bins = [0]*nb
            bins = [0]*nb

            for ii in range(nb):
                lb = ii* nperbin
                if ii == nb-1:
                    ub = num_regions
                else:
                    ub = (ii+1)*nperbin
#                print "LB",lb , ub

                a = np.asarray(count[lb:ub])
                t = np.asarray(gc[lb:ub])
#                print "a",a
#                print "t",t
                local_bins[exp_n,ii] = t.min()
                if ii == nb -1 :
                    local_bins[exp_n,ii+1] = np.max(t)

                mean_bins[ii] = np.median(t)
#                mean_bins[ii] = np.mean(t)
                bins[ii] = mquantiles(a,prob=np.arange(0.1,0.9,0.1))
                print "quantiles",bins[ii]
            print bins
            for qnt in range(8):
                qnt_bns = np.ndarray(len(bins))
                for idx,bb in enumerate(bins):
                    qnt_bns[idx] = bb[qnt]
    #            print "BINS",qnt_bns
                quant_median[qnt]=np.median(qnt_bns)

            print quant_median
            gc_factor = np.zeros(nb,dtype=np.dtype('float'))
            for ii in range(nb):
                offst = np.zeros(len(quant_median),dtype=np.dtype('float'))
                for idx,xx in enumerate(quant_median):
                    offst[idx] = float(bins[ii][idx]) / float(xx+1)
                gc_factor[ii] = np.median(offst)

            print 'MMMMM', gc_factor
            local_meanbins[exp_n] = mean_bins
            local_factor[exp_n] = gc_factor

    mglobals.set_gc_factors(local_bins, local_factor, local_meanbins)


def plot_gc_content():

    idx = 0
    for tissue, list_idx in mglobals.tissue_repl.items():
        pyplot.figure(idx)
        for exp_n in list_idx :
#            f = interpolate.interp1d(mglobals.gc_means[exp_n], mglobals.gc_bins_vaL[exp_n])
            print mglobals.gc_means[exp_n]
            mn = mglobals.gc_means[exp_n].min()
            mx = mglobals.gc_means[exp_n].max()
            xx = np.arange(mn, mx ,0.001)
            yy = mglobals.gc_factor[exp_n](xx)
            print "XX",xx
            print "Yy",yy
            pyplot.plot(xx,yy,label=mglobals.exp_list[exp_n])
            pyplot.title("Gc factor")
            pyplot.grid()
            pyplot.legend(loc='upper left')
#        pyplot.show()
        pyplot.savefig('./gcontent_%s.png'%tissue)
        idx += 1
