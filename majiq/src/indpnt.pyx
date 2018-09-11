import sys
import majiq.src.io as majiq_io
cimport majiq.src.io as majiq_io
import psutil
import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.internals.HetStats cimport HetStats
from majiq.src.internals.qLSV cimport hetLSV, qLSV

from voila.api import Matrix
from voila.constants import ANALYSIS_HETEROGEN, VOILA_FILE_VERSION
from majiq.src.internals.psi cimport psi_distr_t, get_samples_from_psi, get_psi_border, pair_int_t, test_calc

from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from cython.parallel import prange

cimport numpy as np
import numpy as np

def calc_independent(args):
    pipeline_run(independent(args))

cdef int _statistical_test_computation(object out_h5p, dict comparison, list list_of_lsv, vector[string] stats_list,
                                        int psi_samples, map[string, qLSV*] lsv_vec, str outDir, int nthreads,
                                        object logger )  except -1 :
    cdef int nlsv = len(list_of_lsv)
    cdef vector[np.float32_t*] cond1_smpl
    cdef vector[np.float32_t*] cond2_smpl

    cdef np.ndarray[np.float32_t, ndim=2, mode="c"]  k
    cdef object cc
    cdef int cond, xx, i
    cdef str cond_name, lsv, statsnames
    cdef int index, nways, lsv_index
    cdef list file_list = []
    cdef list statlist
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"]  oPvals
    # cdef dict output = {}
    cdef map[string, vector[psi_distr_t]] output
    cdef string lsv_id, stname
    cdef HetStats* StatsObj = new HetStats()
    cdef int nstats

    cdef hetLSV* hetObj_ptr

    if not StatsObj.initialize_statistics(stats_list):
        print('ERROR stats')
        return -1

    statlist = []
    for stname in StatsObj.names:
        statlist.append(stname.decode('utf-8'))

    out_h5p.stat_names = statlist

    logger.info("Using statistics: %s" % " ".join(statlist))
    nstats = StatsObj.get_number_stats()

    index = 0

    print(comparison)

    for cond_name, cond in comparison.items():
        file_list.append([])
        for xx in range(cond):
            cc = np.load(open(get_tmp_psisample_file(outDir, "%s_%s" %(cond_name, xx)), 'rb'))
            file_list[index].append(cc)
        index +=1

    # for lsv in list_of_lsv:
    #     lsv_id = lsv.encode('utf-8')
    #     nways = lsv_vec[lsv_id].get_num_ways()
    #     output[lsv_id] = np.zeros(shape=(nways, nstats), dtype=np.float32)

    for i in prange(nlsv, nogil=True, num_threads=nthreads):
        with gil:
            lsv = list_of_lsv[i]
            lsv_id = lsv.encode('utf-8')
            hetObj_ptr = <hetLSV*> lsv_vec[lsv_id]
            hetObj_ptr.create_condition_samples(len(file_list[0]), len(file_list[1]), psi_samples)
            lsv_index = hetObj_ptr.get_junction_index()
            nways = hetObj_ptr.get_num_ways()

            for fidx, cc in enumerate(file_list[0]):
                k = cc[lsv_index:lsv_index+nways]
                hetObj_ptr.add_condition1(<np.float32_t *> k.data, fidx, nways, psi_samples)

            for fidx, cc in enumerate(file_list[1]):
                k = cc[lsv_index:lsv_index+nways]
                hetObj_ptr.add_condition2(<np.float32_t *> k.data, fidx, nways, psi_samples)

        output[lsv_id] = vector[psi_distr_t](nways, psi_distr_t(nstats))
        test_calc(output[lsv_id], StatsObj, hetObj_ptr, psi_samples, 0.95)
        # test_calc(<np.float32_t *> oPvals.data, StatsObj, hetObj_ptr, psi_samples, 0.95)
        hetObj_ptr.clear()


    logger.info('Storing Voila file statistics')
    for lsv in list_of_lsv:
        lsv_id = lsv.encode('utf-8')
        nways =lsv_vec[lsv_id].get_num_ways()
        oPvals = np.zeros(shape=(nways, nstats), dtype=np.float32)
        for ii in range(nways):
            for jj in range(nstats):
                # print(output[lsv_id][ii][jj])
                oPvals[ii, jj] = output[lsv_id][ii][jj]
        out_h5p.heterogen(lsv).add(junction_stats=oPvals)


cdef int _het_computation(object out_h5p, dict file_cond, list list_of_lsv, map[string, qLSV*] lsv_vec,
                          dict lsv_type_dict, dict junc_info, int psi_samples, int nthreads, int nbins, str outdir,
                          object logger ) except -1:
    cdef string lsv_id
    cdef list cond_list
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

    # cdef list narray_mu, na_postpsi
    # cdef ArrayWrapper mu_w0, mu_w1, ppsi_w0, ppsi_w1

    get_psi_border(psi_border, nbins)

    for lsv in list_of_lsv:
        lsv_id = lsv.encode('utf-8')
        nways = lsv_vec[lsv_id].get_num_ways()
        total_njuncs += nways


    for cidx, (cond_name, cond_list) in enumerate(file_cond.items()):
        max_nfiles = max(max_nfiles, len(cond_list))
        for fidx, f in enumerate(cond_list):
            osamps = np.zeros(shape=(total_njuncs, psi_samples), dtype=np.float32)
            majiq_io.get_coverage_mat_lsv(lsv_vec, [f], "", nthreads)
            for i in prange(nlsv, nogil=True, num_threads=nthreads):
                with gil:
                    lsv = list_of_lsv[i]
                    lsv_id = lsv.encode('utf-8')
                get_samples_from_psi(<np.float32_t *> osamps.data, <hetLSV*> lsv_vec[lsv_id], psi_samples, psi_border,
                                     nbins, cidx, fidx)
            fname = get_tmp_psisample_file(outdir, "%s_%s" %(cond_name, fidx) )
            majiq_io.dump_hettmp_file(fname, osamps)


    logger.info("Store Voila LSV information")

    for lsv in list_of_lsv:
        lsv_id = lsv.encode('utf-8')
        nways = lsv_vec[lsv_id].get_num_ways()
        mupsi = np.ndarray(shape=(len(file_cond), max_nfiles, nways), dtype=np.float32, order="c")
        postpsi = np.ndarray(shape=(len(file_cond), nways, nbins), dtype=np.float32, order="c")
        mupsi.fill(-1)
        hetObj_ptr = <hetLSV*> lsv_vec[lsv_id]
        for x in range(len(file_cond)):
            for y in range(len(file_cond[x])):
                for z in range(nways):
                    mupsi[x,y,z] = hetObj_ptr.mu_psi[x][y][z]

        for x in range(len(file_cond)):
            for y in range(nways):
               for z in range(nbins):
                   postpsi[x,y,z] = hetObj_ptr.post_psi[x][y][z]

        out_h5p.heterogen(lsv).add(lsv_type=lsv_type_dict[lsv], mu_psi=mupsi,mean_psi=postpsi, junctions=junc_info[lsv])

cdef void _core_independent(object self):


    cdef dict junc_info = {}
    cdef dict lsv_type_dict = {}
    cdef object logger
    cdef int nbins = 40
    cdef bint is_ir
    cdef string lsv_id
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

    majiq_logger.create_if_not_exists(self.outDir)
    logger = majiq_logger.get_logger("%s/het_majiq.log" % self.outDir, silent=self.silent,
                                     debug=self.debug)

    logger.info("Majiq deltapsi heterogeneous v%s" % VERSION)
    logger.info("Command: %s" % " ".join(sys.argv))
    logger.info("GROUP1: %s" % self.files1)
    logger.info("GROUP2: %s" % self.files2)

    try:
        for stats_name in self.stats:
            stats_vec.push_back(stats_name.upper().encode('utf-8'))
            # module_ = __import__('majiq.src.stats.' + stats_name.lower(), fromlist=stats_name.title())
            # class_ = getattr(module_, stats_name.title())
            # operator[stats_name] = class_()
    except ImportError:
        logger.error("The %s statistic is not one of the available statistics, ")
                     #"in  [ %s ]" % (stats_name, ' | '.join(all_stats)))
        return

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

    logger.info('Computation done, saving results....')
    with Matrix(get_quantifier_voila_filename(self.outDir, self.names, het=True), 'w') as out_h5p:
        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_HETEROGEN
        out_h5p.group_names = self.names
        out_h5p.experiment_names = [exps1, exps2]

        j_offset = 0
        for lsv in list_of_lsv:
            nways = len(lsv_type_dict[lsv].split('|')) -1
            # tpair =  pair_int_t(nways, j_offset)
            is_ir = 'i' in lsv_type_dict[lsv]
            m = new hetLSV(nways, j_offset, max_nfiles, nbins, is_ir, len(file_cond))
            lsv_map[lsv.encode('utf-8')] = <qLSV*> m
            j_offset += nways

        logger.info('Sampling from PSI')
        _het_computation(out_h5p, file_cond, list_of_lsv, lsv_map, lsv_type_dict, junc_info, self.psi_samples,
                         nthreads, nbins, self.outDir, logger)

        logger.info('Calculating statistics pvalues')
        _statistical_test_computation(out_h5p, comparison, list_of_lsv, stats_vec, self.psi_samples, lsv_map,
                                      self.outDir, nthreads, logger)

    if self.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        logger.info("Max Memory used %.2f MB" % mem_allocated)
    logger.info("Majiq Heterogeneous calculation for %s_%s ended succesfully! "
                "Result can be found at %s" % (self.names[0], self.names[1], self.outDir))



class independent(BasicPipeline):

    def run(self):
        _core_independent(self)

