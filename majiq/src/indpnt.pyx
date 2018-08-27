import sys
import majiq.src.io as majiq_io
cimport majiq.src.io as majiq_io
import psutil
from majiq.src.psi import heterogen_posterior
import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
# from majiq.src.stats import operator, all_stats
from majiq.src.internals.HetStats cimport HetStats

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

cdef void _statistical_test_computation(object out_h5p, dict comparison, list list_of_lsv, vector[string] stats_list,
                                        int psi_samples, map[string, pair_int_t] lsv_vec, str outDir, int nthreads ) :
    cdef int nlsv = len(list_of_lsv)
    cdef vector[float*] cond1_smpl
    cdef vector[float*] cond2_smpl

    cdef np.ndarray[np.float32_t, ndim=2, mode="c"]  k
    cdef object cc
    cdef int cond, xx, i
    cdef str cond_name, lsv
    cdef int index, nways, lsv_index
    cdef list file_list = []
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"]  oPvals
    cdef dict output = {}
    cdef string lsv_id
    cdef HetStats* StatsObj = new HetStats()
    cdef int nstats = stats_list.size()


    if not StatsObj.initialize_statistics(stats_list):
        print('ERROR stats')
        return

    index = 0

    print(comparison)

    for cond_name, cond in comparison.items():
        file_list.append([])
        for xx in range(cond):
            cc = np.load(open(get_tmp_psisample_file(outDir, "%s_%s" %(cond_name, xx)), 'rb'))
            file_list[index].append(cc)
        index +=1

    for lsv in list_of_lsv:
        lsv_id = lsv.encode('utf-8')
        nways = lsv_vec[lsv_id].first
        output[lsv_id] = np.zeros(shape=(nways, nstats), dtype=np.float32)



    for i in prange(nlsv, nogil=True, num_threads=nthreads):
        with gil:
            lsv = list_of_lsv[i]
            # print (i, lsv)
            lsv_id = lsv.encode('utf-8')
            lsv_index = lsv_vec[lsv_id].second
            nways = lsv_vec[lsv_id].first
            oPvals = output[lsv_id]
            # oPvals = np.zeros(shape=(nways, nstats), dtype=np.float32)
            # output[lsv_id] = oPvals
            for cc in file_list[0]:
                k = cc[lsv_index:lsv_index+nways]
                cond1_smpl.push_back(<np.float32_t *> k.data)

            for cc in file_list[1]:
                k = cc[lsv_index:lsv_index+nways]
                cond2_smpl.push_back(<np.float32_t *> k.data)

        test_calc(<np.float32_t *> oPvals.data, cond1_smpl, cond2_smpl, StatsObj, nways, psi_samples, 0.95)
        cond1_smpl.clear()
        cond2_smpl.clear()
        # with gil:
        #     print('END KOLA', i)

    print('DUMP VOILA FILE')
    for lsv in list_of_lsv:
        lsv_id = lsv.encode('utf-8')
        out_h5p.heterogen(lsv).add(junction_stats=output[lsv_id])



cdef int _het_computation(object out_h5p, dict file_cond, list list_of_lsv, map[string, pair_int_t] lsv_vec,
                           dict lsv_type_dict, dict junc_info, int psi_samples, int nthreads, int nbins, str outdir) except -1:
    cdef string lsv_id
    cdef str f, cond_name, fname ;
    cdef int cidx, fidx
    cdef map[string, vector[psi_distr_t]] cov_dict
    cdef int nways
    cdef int nlsv = len(list_of_lsv)
    cdef map[string, int] lsv_map

    cdef dict out_mupsi_d = {}
    cdef dict out_postpsi_d = {}
    # cdef np.float32_t * o_postpsi
    # cdef np.float32_t * o_mupsi

    cdef int j_offset = 0
    cdef psi_distr_t psi_border
    cdef max_nfiles = 0
    cdef list cond_list
    cdef int total_njuncs = 0
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] osamps
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] o_mupsi
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] o_postpsi
    cdef int i, msamples
    cdef bint is_ir

    for cond_name, cond_list in file_cond.items():
        max_nfiles = max(max_nfiles, len(cond_list))

    psi_border = get_psi_border(nbins)
    for lsv in list_of_lsv:
        lsv_id = lsv.encode('utf-8')
        nways = len(lsv_type_dict[lsv].split('|')) -1
        out_mupsi_d[lsv_id]= np.ndarray(shape=(len(file_cond), max_nfiles, nways), dtype=np.float32)
        out_mupsi_d[lsv_id].fill(-1)
        out_postpsi_d[lsv_id] = np.zeros(shape=(len(file_cond), nways, nbins), dtype=np.float32)
        lsv_map[lsv.encode('utf-8')] = nways

        j_offset += nways
    total_njuncs = j_offset

    for cidx, (cond_name, cond_list) in enumerate(file_cond.items()):


        for fidx, f in enumerate(cond_list):
            osamps = np.zeros(shape=(total_njuncs, psi_samples), dtype=np.float32)
            majiq_io.get_coverage_mat(cov_dict, lsv_map, [f], "", nthreads)
            for i in prange(nlsv, nogil=True, num_threads=nthreads):
                with gil:
                    lsv = list_of_lsv[i]
                    lsv_id = lsv.encode('utf-8')

                    nways = cov_dict[lsv_id].size()
                    msamples = cov_dict[lsv_id][0].size()
                    # print(i, lsv, lsv_type_dict[lsv], nways, cidx)
                    o_mupsi = out_mupsi_d[lsv_id][cidx]
                    o_postpsi = out_postpsi_d[lsv_id][cidx]
                    is_ir = 'i' in lsv_type_dict[lsv]

                get_samples_from_psi(cov_dict[lsv_id], <np.float32_t *> osamps.data, <np.float32_t *> o_mupsi.data,
                                     <np.float32_t *> o_postpsi.data, psi_samples, lsv_vec[lsv_id].second, psi_border,
                                     nways, msamples, nbins, is_ir)

            fname = get_tmp_psisample_file(outdir, "%s_%s" %(cond_name, fidx) )
            majiq_io.dump_hettmp_file(fname, osamps)

    # print("Dump psi_samples")
    for lsv in list_of_lsv:
        lsv_id = lsv.encode('utf-8')
        out_h5p.heterogen(lsv).add(lsv_type=lsv_type_dict[lsv], mu_psi=out_mupsi_d[lsv_id], mean_psi=out_postpsi_d[lsv_id],
                                      junctions=junc_info[lsv])


cdef void _core_independent(object self):


    cdef dict junc_info = {}
    cdef dict lsv_type_dict = {}
    cdef object logger
    cdef int nbins = 40
    cdef bint is_ir
    cdef string lsv_id
    cdef int nways, msamples, i, j_offset
    cdef list list_of_lsv
    cdef map[string, pair_int_t] lsv_vec
    cdef dict file_cond = {self.names[0]: self.files1, self.names[1]: self.files2}
    cdef pair_int_t tpair
    cdef vector[string] stats_vec ;
    cdef dict comparison = {self.names[0]: len(self.files1), self.names[1]: len(self.files2)}

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


    # list_of_lsv = ['ENSMUSG00000032735:s:61823867-61823938']
    nlsv = len(list_of_lsv)
    if nlsv == 0:
        logger.info("There is no LSVs that passes the filters")
        return

    logger.info("Number quantifiable LSVs: %s" % nlsv)
    nthreads = min(self.nthreads, nlsv)


    logger.info('Computation done, saving results....')
    with Matrix(get_quantifier_voila_filename(self.outDir, self.names, het=True), 'w') as out_h5p:
        out_h5p.file_version = VOILA_FILE_VERSION
        out_h5p.analysis_type = ANALYSIS_HETEROGEN
        out_h5p.group_names = self.names
        out_h5p.experiment_names = [exps1, exps2]
        out_h5p.stat_names = self.stats
        j_offset = 0
        for lsv in list_of_lsv:
            nways = len(lsv_type_dict[lsv].split('|')) -1
            tpair =  pair_int_t(nways, j_offset)
            lsv_vec[lsv.encode('utf-8')] = tpair
            j_offset += nways

        logger.info('Sampling from PSI')
        _het_computation(out_h5p, file_cond, list_of_lsv, lsv_vec, lsv_type_dict, junc_info, self.psi_samples,
                         nthreads, nbins, self.outDir)

        logger.info('Calculating statistics pvalues')
        _statistical_test_computation(out_h5p, comparison, list_of_lsv, stats_vec, self.psi_samples, lsv_vec,
                                      self.outDir, nthreads)

        #
        #
        # for lsv in list_of_lsv:
        #     lsv_id = lsv.encode('utf-8')
        #     lsv_id = results[2]
        #     lsv_type = lsv_type_dict[lsv_id]
        #     mu_psi = results[1][0]
        #     mean_psi = results[1][1]
        #     out_h5p.heterogen(lsv_id).add(lsv_type=lsv_type, mu_psi=mu_psi, mean_psi=mean_psi, junction_stats=results[0],
        #                                  junctions=junc_info[lsv_id])



    if self.mem_profile:
        mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
        logger.info("Max Memory used %.2f MB" % mem_allocated)
    logger.info("Majiq Heterogeneous calculation for %s_%s ended succesfully! "
                "Result can be found at %s" % (self.names[0], self.names[1], self.outDir))



class independent(BasicPipeline):

    def run(self):
        _core_independent(self)

