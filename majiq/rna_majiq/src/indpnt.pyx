import os
import shutil
import sys
import rna_majiq.src.io as majiq_io
cimport rna_majiq.src.io as majiq_io
import psutil
import rna_majiq.src.logger as majiq_logger
from rna_majiq.src.basic_pipeline import BasicPipeline, pipeline_run
import rna_majiq.src.constants as constants
from rna_majiq.src.internals.mtypes cimport (
    pair_int_t,
    psi_distr_t,
)
from rna_majiq.src.internals.HetStats cimport HetStats
from rna_majiq.src.internals.qLSV cimport hetLSV, qLSV
from rna_majiq.src.internals.psi cimport get_samples_from_psi, get_psi_border, test_calc, mt19937

from rna_voila.api import Matrix
from rna_voila.constants import ANALYSIS_HETEROGEN, VOILA_FILE_VERSION

from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from cython.parallel import prange, threadid

cimport numpy as np
import numpy as np

def calc_independent(args):
    pipeline_run(independent(args))

cdef int _statistical_test_computation(object out_h5p, dict comparison, list list_of_lsv, vector[string] stats_list,
                                       int psi_samples, map[string, qLSV*] lsv_vec, str outDir, int nthreads,
                                       float test_percentile, object logger )  except -1 :
    cdef int nlsv = len(list_of_lsv)
    cdef vector[np.float32_t*] cond1_smpl
    cdef vector[np.float32_t*] cond2_smpl

    cdef np.ndarray[np.float32_t, ndim=2, mode="c"]  k
    cdef object cc, cc_memmap
    cdef int cond, xx, i, lsv_chunk_idx
    cdef int lsv_min_idx, lsv_max_idx, junction_min, junction_max
    cdef str cond_name, statsnames
    cdef int index, nways, lsv_index
    cdef list file_list = []
    cdef list statlist
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] np_psisamples_pval_percentile
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] np_psimean_pval
    cdef map[string, vector[psi_distr_t]] psisamples_pval_percentile
    cdef map[string, vector[psi_distr_t]] psimean_pval
    cdef map[string,psi_distr_t] oScore
    cdef np.ndarray[np.float32_t, ndim=1, mode="c"] score
    cdef string lsv_id, stname
    cdef HetStats* StatsObj = new HetStats()
    cdef int nstats

    cdef hetLSV* hetObj_ptr
    cdef hetLSV* lsv_min
    cdef hetLSV* lsv_max

    if not StatsObj.initialize_statistics(stats_list):
        logger.error('Unable to initialize statistics object')
        return -1

    statlist = []
    for stname in StatsObj.names:
        statlist.append(stname.decode('utf-8'))

    out_h5p.stat_names = statlist

    logger.info("Using statistics: %s" % " ".join(statlist))
    nstats = StatsObj.get_number_stats()

    logger.info(f"Comparison: {comparison}")
    # XXX NOTE: chunking of LSVs assumes that `list_of_lsv` LSVs are
    # aligned/share order with junctions in psi samples files

    # get number of LSVs to read in at a time
    # 250000 magic number -- roughly 0.5-1.5G total depending on average number
    # of junctions per LSV
    cdef int lsv_chunksize = 1 + 250000 // sum(comparison.values())
    # get indices of LSVs to read in at a time
    cdef list lsv_chunks = list(range(0, nlsv, lsv_chunksize))
    if lsv_chunks[-1] < nlsv:
        lsv_chunks.append(nlsv)

    for lsv_chunk_idx in range(len(lsv_chunks) - 1):
        # get indices of junctions in LSVs to process in this chunk
        lsv_min_idx = lsv_chunks[lsv_chunk_idx]
        lsv_max_idx = lsv_chunks[lsv_chunk_idx + 1]
        lsv_min = <hetLSV*> lsv_vec[list_of_lsv[lsv_min_idx]]
        lsv_max = <hetLSV*> lsv_vec[list_of_lsv[lsv_max_idx - 1]]
        junction_min = lsv_min.get_junction_index()
        junction_max = lsv_max.get_junction_index() + lsv_max.get_num_ways()
        if len(lsv_chunks) > 2:
            logger.info(
                f"Calculating statistics for LSVs {lsv_min_idx} to {lsv_max_idx}"
                f" out of {nlsv} (progress:"
                f" {lsv_min_idx / nlsv:.1%}-{lsv_max_idx / nlsv:.1%})"
            )
        # load LSVs to process in this chunk
        index = 0
        for cond_name, cond in comparison.items():
            file_list.append([])
            for xx in range(cond):
                cc_memmap = np.load(
                    constants.get_tmp_psisample_file(outDir, "%s_%s" %(cond_name, xx)), "r"
                )
                # load contiguous chunk of the array all at once, add to filelist
                cc = np.array(cc_memmap[junction_min:junction_max])
                file_list[index].append(cc)
                # close the memory map (note cc_memmap unusable until new np.load)
                cc_memmap._mmap.close()
            index +=1
        # for each of the LSVs in this chunk... (in parallel)
        for i in prange(lsv_min_idx, lsv_max_idx, nogil=True, num_threads=nthreads):
            with gil:
                lsv_id = list_of_lsv[i]
                hetObj_ptr = <hetLSV*> lsv_vec[lsv_id]
                hetObj_ptr.create_condition_samples(len(file_list[0]), len(file_list[1]), psi_samples)
                # get index into chunked arrays (offset by junction_min)
                lsv_index = hetObj_ptr.get_junction_index() - junction_min
                nways = hetObj_ptr.get_num_ways()

                for fidx, cc in enumerate(file_list[0]):
                    k = cc[lsv_index:lsv_index+nways]
                    hetObj_ptr.add_condition1(<np.float32_t *> k.data, fidx, nways, psi_samples)

                for fidx, cc in enumerate(file_list[1]):
                    k = cc[lsv_index:lsv_index+nways]
                    hetObj_ptr.add_condition2(<np.float32_t *> k.data, fidx, nways, psi_samples)

                # insertions into maps are not thread-safe
                psisamples_pval_percentile[lsv_id] = vector[psi_distr_t](nways, psi_distr_t(nstats))
                psimean_pval[lsv_id] = vector[psi_distr_t](nways, psi_distr_t(nstats))
                oScore[lsv_id] = psi_distr_t(nways)
            test_calc(psimean_pval[lsv_id], psisamples_pval_percentile[lsv_id], oScore[lsv_id], StatsObj, hetObj_ptr, psi_samples, test_percentile)
            hetObj_ptr.clear()
        # clear out loaded arrays
        for cc_list in file_list:
            cc_list.clear()
        file_list.clear()


    logger.info('Storing Voila file statistics')
    for lsv in list_of_lsv:
        nways =lsv_vec[lsv].get_num_ways()
        np_psisamples_pval_percentile = np.zeros(shape=(nways, nstats), dtype=np.float32)
        np_psimean_pval = np.zeros(shape=(nways, nstats), dtype=np.float32)
        score = np.zeros(shape=nways, dtype=np.float32)
        for ii in range(nways):
            for jj in range(nstats):
                np_psisamples_pval_percentile[ii, jj] = psisamples_pval_percentile[lsv][ii][jj]
                np_psimean_pval[ii, jj] = psimean_pval[lsv][ii][jj]
            score[ii] = oScore[lsv][ii]
        out_h5p.heterogen(lsv.decode('utf-8')).add(
            junction_stats=np_psimean_pval,
            junction_psisamples_stats=np_psisamples_pval_percentile,
            tnom_score=score
        )


cdef int _het_computation(
    object out_h5p, dict file_cond, list list_of_lsv,
    map[string, qLSV*] lsv_vec, dict lsv_type_dict, dict junc_info,
    int psi_samples, int nthreads, int nbins, str outdir, int minreads,
    int minnonzero, float visualization_var_max, object logger
) except -1:
    cdef string lsv
    cdef list cond_list, conditions
    cdef str f, cond_name, fname ;
    cdef int cidx, fidx, i, msamples, nways

    # cdef map[string, qLSV] cov_dict
    cdef int nlsv = len(list_of_lsv)

    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] osamps
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] mupsi
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] postpsi

    cdef psi_distr_t psi_border = psi_distr_t(nbins+1)
    cdef int max_nfiles = 0
    cdef int total_njuncs = 0
    cdef bint is_ir
    cdef hetLSV* hetObj_ptr
    cdef int visualization_samples  # done per condition given var_max and experiments

    # cdef list narray_mu, na_postpsi
    # cdef ArrayWrapper mu_w0, mu_w1, ppsi_w0, ppsi_w1

    get_psi_border(psi_border, nbins)

    for lsv in list_of_lsv:
        nways = lsv_vec[lsv].get_num_ways()
        total_njuncs += nways

    # initialize source of randomness
    cdef vector[mt19937] generators = vector[mt19937](nthreads)
    cdef int thread_idx
    for thread_idx in range(nthreads):
        generators[thread_idx].seed(constants.HET_SAMPLING_SEED + thread_idx)

    conditions = list(file_cond.keys())
    for cidx, cond_name in enumerate(conditions):
        cond_list = file_cond[cond_name]
        # posterior estimate has variance < 1 / (4 * num_experiments * num_samples_per_exp)
        visualization_samples = int(max(np.ceil(1 / (4 * len(cond_list) * visualization_var_max)), 1))
        logger.info(
            f"Group {cond_name} performing at least {visualization_samples} psi"
            " samples per experiment to bound visualization error"
        )
        max_nfiles = max(max_nfiles, len(cond_list))
        for fidx, f in enumerate(cond_list):
            logger.info(
                f"(group {cidx + 1}/2, experiment {fidx + 1}/{len(cond_list)})"
                f" Group {cond_name} sampling PSI from '{f}'"
            )
            # osamps[..., 0] corresponds to posterior mean, osamps[..., 1:] to psisamples
            osamps = np.zeros(shape=(total_njuncs, psi_samples + 1), dtype=np.float32)
            # lsvs should be made unpassed, require get_coverage_mat_lsv to pass LSVs using minreads/minnonzero
            for i in prange(nlsv, nogil=True, num_threads=nthreads):
                with gil:
                    lsv = list_of_lsv[i]
                lsv_vec[lsv].set_bool(False)
            majiq_io.get_coverage_mat_lsv(lsv_vec, [f], nthreads, False, minreads, minnonzero)
            for i in prange(nlsv, nogil=True, num_threads=nthreads):
                with gil:
                    lsv = list_of_lsv[i]

                thread_idx = threadid()

                get_samples_from_psi(
                    <np.float32_t *> osamps.data, <hetLSV*> lsv_vec[lsv],
                    psi_samples, visualization_samples,
                    psi_border, nbins, cidx, fidx, generators[thread_idx]
                )
                lsv_vec[lsv].reset_samps()
            fname = constants.get_tmp_psisample_file(outdir, "%s_%s" %(cond_name, fidx) )
            majiq_io.dump_hettmp_file(fname, osamps)


    logger.info("Store Voila LSV information")

    for lsv in list_of_lsv:
        nways = lsv_vec[lsv].get_num_ways()
        mupsi = np.ndarray(shape=(nways, len(file_cond), max_nfiles), dtype=np.float32, order="c")
        postpsi = np.ndarray(shape=(nways, len(file_cond), nbins), dtype=np.float32, order="c")
        mupsi.fill(-1)
        hetObj_ptr = <hetLSV*> lsv_vec[lsv]
        for z in range(nways):
            for x in range(len(file_cond)):
                for y in range(len(file_cond[conditions[x]])):
                    mupsi[z, x, y] = hetObj_ptr.mu_psi[x][y][z]

        for y in range(nways):
            for x in range(len(file_cond)):
                for z in range(nbins):
                    postpsi[y, x, z] = hetObj_ptr.post_psi[x][y][z]
                postpsi[y, x] /= postpsi[y, x].sum()

        out_h5p.heterogen(lsv.decode('utf-8')).add(lsv_type=lsv_type_dict[lsv][0].decode('utf-8'), mu_psi=mupsi,
                                                   mean_psi=postpsi, junctions=junc_info[lsv])

cdef void _core_independent(object self):

    cdef dict junc_info = {}
    cdef dict lsv_type_dict = {}
    cdef object logger
    cdef int nbins = 40
    cdef bint is_ir
    cdef string lsv
    cdef int nways, msamples, i, j_offset
    cdef list list_of_lsv
    # cdef map[string, pair_int_t] lsv_vec
    cdef map[string, qLSV*] lsv_map
    cdef dict file_cond = {self.names[0]: self.files1, self.names[1]: self.files2}
    cdef pair_int_t tpair
    cdef vector[string] stats_vec ;
    cdef dict comparison = {self.names[0]: len(self.files1), self.names[1]: len(self.files2)}
    cdef hetLSV* m
    cdef int max_nfiles = 0
    cdef float visualization_var_max = self.visualization_std * self.visualization_std

    majiq_logger.create_if_not_exists(self.outDir)
    logFile = self.logger if self.logger else f"{self.outDir}/het_majiq.log"
    logger = majiq_logger.get_logger(logFile, silent=self.silent, debug=self.debug)

    logger.info(f"Majiq deltapsi heterogeneous v{constants.VERSION}")
    logger.info("Command: %s" % " ".join(sys.argv))
    logger.info("GROUP1: %s" % self.files1)
    logger.info("GROUP2: %s" % self.files2)

    for stats_name in self.stats:
        stats_vec.push_back(stats_name.upper().encode('utf-8'))

    list_of_lsv1, exps1 = majiq_io.extract_lsv_summary(self.files1, types_dict=lsv_type_dict,
                                                       minnonzero=self.minpos, min_reads=self.minreads,
                                                       junc_info=junc_info, percent=self.min_exp, logger=logger)

    logger.info("Group %s: %s LSVs" % (self.names[0], len(list_of_lsv1)))

    list_of_lsv2, exps2 = majiq_io.extract_lsv_summary(self.files2, types_dict=lsv_type_dict,
                                                       minnonzero=self.minpos, min_reads=self.minreads,
                                                       junc_info=junc_info, percent=self.min_exp, logger=logger)
    logger.info("Group %s: %s LSVs" % (self.names[1], len(list_of_lsv1)))

    list_of_lsv = list(set(list_of_lsv1).intersection(set(list_of_lsv2)))
    nlsv = len(list_of_lsv)
    if nlsv == 0:
        logger.info("There is no LSVs that passes the filters")
        return

    logger.info("Number quantifiable LSVs: %s" % nlsv)
    nthreads = min(self.nthreads, nlsv)

    for cond_name, cond_list in file_cond.items():
        max_nfiles = max(max_nfiles, len(cond_list))

    with Matrix(constants.get_quantifier_voila_filename(self.outDir, self.names, het=True), 'w') as out_h5p:
        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_HETEROGEN
        out_h5p.group_names = self.names
        out_h5p.experiment_names = [exps1, exps2]
        out_h5p.psi_samples = self.psi_samples
        out_h5p.test_percentile = self.test_percentile

        j_offset = 0
        for lsv in list_of_lsv:
            nways = lsv_type_dict[lsv][1]
            is_ir = b'i' in lsv_type_dict[lsv][0]
            m = new hetLSV(nways, j_offset, max_nfiles, nbins, is_ir, len(file_cond))
            lsv_map[lsv] = <qLSV*> m
            j_offset += nways

        logger.info('Sampling from PSI')
        _het_computation(
            out_h5p, file_cond, list_of_lsv, lsv_map, lsv_type_dict, junc_info,
            self.psi_samples, nthreads, nbins, self.outDir, self.minreads,
            self.minpos, visualization_var_max, logger
        )

        logger.info('Calculating statistics pvalues')
        _statistical_test_computation(out_h5p, comparison, list_of_lsv, stats_vec, self.psi_samples, lsv_map,
                                      self.outDir, nthreads, self.test_percentile, logger)


    if not self.keep_tmpfiles:
        shutil.rmtree(constants.get_tmp_dir(self.outDir))

    if self.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        logger.info("Max Memory used %.2f MB" % mem_allocated)
    logger.info(f"Majiq Heterogeneous calculation for {self.names[0]}{constants.GROUP_NAME_SEP}{self.names[1]} ended successfully! "
                f"Result can be found at {self.outDir}")



class independent(BasicPipeline):

    def run(self):
        _core_independent(self)

