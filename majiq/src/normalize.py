import majiq.src.io_utils

__author__ = 'jordi@biociphers.org'

from itertools import izip
import os
import sys
import numpy as np
from numpy.ma import masked_less
import scipy.sparse
from scipy.stats.mstats_basic import mquantiles

from majiq.src import polyfitnb as majiqfit
import majiq.src.config as majiq_config
import majiq.src.io_utils as majiq_io_utils


def mark_stacks(lsv_list, fitfunc_r, pvalue_limit, logger=None):

    if pvalue_limit < 0:
        return
    logger.debug("Marking and masking stacks")
    minstack = sys.maxint
    # the minimum value marked as stack
    numstacks = 0
    for lidx, junctions in enumerate(lsv_list[0]):

        for i, junction in enumerate(junctions):
            if np.count_nonzero(junction) == 0:
                continue
            for j, value in enumerate(junction):
                if value > 0:
                    # TODO Use masker, and marking stacks will probably be faster.
                    copy_junc = list(junction)
                    copy_junc.pop(j)
                    copy_junc = np.array(copy_junc)
                    copy_junc = copy_junc[copy_junc > 0]
                    nzpos = np.count_nonzero(copy_junc)

                    #FINISH TODO
                    mean_rest = np.mean(copy_junc) * nzpos
                    pval = majiqfit.get_negbinom_pval(fitfunc_r, mean_rest, value)
                    if pval < pvalue_limit:
                        lsv_list[0][lidx][i, j] = -2
                        minstack = min(minstack, value)
                        numstacks += 1
        masked_less(lsv_list[0][lidx], 0)

    if logger:
        logger.debug("Out of %s values, %s marked as stacks with a p-value threshold of %s (%.3f%%)"
                    % (junctions.size, numstacks, pvalue_limit, (float(numstacks) / junctions.size) * 100))
    return lsv_list


def __gc_factor_ind(val, exp_idx):
    res = 0
    for ii, jj in enumerate(majiq_config.gc_bins[exp_idx]):
        if val < jj:
            res = ii
    return res


def gc_content_norm(self, lsv_list, const_list):
    """Normalize the matrix using the gc content"""
    self.logger.debug("GC content normalization...")
    if self.gcnorm:
        for lidx, lsv in enumerate(lsv_list[0]):
            lsv_list[0][lidx] = np.multiply(lsv, lsv_list[2][lidx])
        const_list[0] = np.multiply(const_list[0], const_list[2])
    return lsv_list, const_list


def gc_factor_calculation(nb):

    local_bins = np.zeros(shape=(majiq_config.num_experiments, nb+1), dtype=np.dtype('float'))
    local_meanbins = np.zeros(shape=(majiq_config.num_experiments, nb),   dtype=np.dtype('float'))
    local_factor = np.zeros(shape=(majiq_config.num_experiments, nb),   dtype=np.dtype('float'))

    gc_pairs = {'GC': [[] for xx in xrange(majiq_config.num_experiments)],
                'COV': [[] for xx in xrange(majiq_config.num_experiments)]}

    # read local files
    for chnk in range(majiq_config.num_final_chunks):
        temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
        yfile = '%s/gccontent.temppkl' % temp_dir
        if not os.path.exists(yfile):
            continue
        gc_c = majiq.src.io_utils.load_bin_file(yfile)
        for exp_n in xrange(majiq_config.num_experiments):
            gc_pairs['GC'][exp_n].extend(gc_c['GC'][exp_n])
            gc_pairs['COV'][exp_n].extend(gc_c['COV'][exp_n])

    #print mglobals.tissue_repl
    for tissue, list_idx in majiq_config.tissue_repl.items():
        for exp_n in list_idx:
            count = gc_pairs['COV'][exp_n]
            gc = gc_pairs['GC'][exp_n]

            if len(gc) == 0:
                continue

            count, gc = izip(*sorted(izip(count, gc), key=lambda x: x[1]))

            num_regions = len(count)
            nperbin = num_regions / nb

            quant_median = [0.0]*8
            mean_bins = [0]*nb
            bins = [0]*nb

            for ii in range(nb):
                lb = ii * nperbin
                if ii == nb-1:
                    ub = num_regions
                else:
                    ub = (ii+1) * nperbin

                a = np.asarray(count[lb:ub])
                t = np.asarray(gc[lb:ub])

                try:
                    local_bins[exp_n, ii] = t.min()
                except ValueError:
                    local_bins[exp_n, ii] = 0
                if ii == nb - 1:
                    local_bins[exp_n, ii+1] = np.max(t)

                #mean_bins[ii] = np.median(t)
                mean_bins[ii] = np.mean(t)
                bins[ii] = mquantiles(a, prob=np.arange(0.1, 0.9, 0.1))
                #print "quantiles", bins[ii]

            for qnt in range(8):
                qnt_bns = np.ndarray(len(bins))
                for idx, bb in enumerate(bins):
                    qnt_bns[idx] = bb[qnt]
                #print "BINS", qnt_bns
                #quant_median[qnt]=np.median(qnt_bns)
                quant_median[qnt] = np.mean(qnt_bns)

            #print quant_median
            gc_factor = np.zeros(nb, dtype=np.dtype('float'))
            for ii in range(nb):
                offst = np.zeros(len(quant_median), dtype=np.dtype('float'))
                for idx, xx in enumerate(quant_median):
                    offst[idx] = float(bins[ii][idx]) / float(xx)
                gc_factor[ii] = 1/np.mean(offst)

            #print 'MMMMM', gc_factor
            local_meanbins[exp_n] = mean_bins
            local_factor[exp_n] = gc_factor

    majiq_config.set_gc_factors(local_bins, local_factor, local_meanbins)


def prepare_gc_content(gene_list, temp_dir):
    gc_pairs = {'GC': [[] for xx in xrange(majiq_config.num_experiments)],
                'COV': [[] for xx in xrange(majiq_config.num_experiments)]}

    for gn in gene_list:
        for ex in gn.get_exon_list():
            gc_val = ex.get_gc_content()
            st, end = ex.get_coordinates()
            if gc_val == 0 or end - st < 30:
                continue
            for exp_n in xrange(majiq_config.num_experiments):
                cov = ex.get_coverage(exp_n)
                if cov < 1:
                    continue
                gc_pairs['GC'][exp_n].append(gc_val)
                gc_pairs['COV'][exp_n].append(cov)

    fname = '%s/gccontent.temppkl' % temp_dir
    majiq_io_utils.dump_bin_file(gc_pairs, fname)