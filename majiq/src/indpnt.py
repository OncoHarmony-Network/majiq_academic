import multiprocessing as mp
import sys

import majiq.src.io as majiq_io
import psutil
from majiq.src.psi import heterogen_posterior

import majiq.src.logger as majiq_logger
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper, chunks
from majiq.src.stats import operator, all_stats
from voila.api import Matrix
from voila.constants import ANALYSIS_HETEROGEN, VOILA_FILE_VERSION


def het_quantification(list_of_lsv, chnk, conf, logger):
    logger.info("Quantifying LSVs PSI.. %s" % chnk)
    num_exp = [len(conf.files1), len(conf.files2)]

    f_list = [None, None]
    f_list[0] = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files1)
    f_list[1] = majiq_io.get_extract_lsv_list(list_of_lsv, conf.files2)

    for lidx, lsv_id in enumerate(list_of_lsv):
        if lidx % 50 == 0:
            print("Event %d ..." % lidx)
            sys.stdout.flush()

        boots = [np.array(f_list[0][lsv_id].coverage), np.array(f_list[1][lsv_id].coverage)]
        msamples = boots[0].shape[2]
        # stat_het = []
        grp_het = []
        assert boots[0].shape[2] == boots[1].shape[2], "LSV %s, has different types in %s and %s (%s vs %s). " \
                                                       "Please check that the conditions has been build together." \
                                                       % (lsv_id, conf.names[0], conf.names[1],
                                                          boots[0].shape[2], boots[1].shape[2])
        samps = heterogen_posterior(boots, grp_het, msamples, conf.nsamples, conf.vwindow, num_exp,
                                    conf.nbins, conf.lsv_type_dict[lsv_id])

        out_stats = do_test_stats(samps, conf.stats, conf.minsamps)
        # for stat_idx in range(out_stats.shape[1]):
        #     stat_het.append(out_stats[:, stat_idx])

        # if num_ways == 2:
        #     break

        qm = QueueMessage(QUEUE_MESSAGE_HETER_DELTAPSI, (out_stats, grp_het, lsv_id), chnk)
        conf.queue.put(qm, block=True)


def do_test_stats(insamps, stats, minsamps):
    # def do_test_stats(lsv_attrs, samples, stats, minsamps):
    """
    Log P values for the given set of statistical tests
    :param lsv_attrs: 2-tuple (lsv_id, lsv_type)
    :param samples: 3d numpy array with dimensions (nfiles, njunc, nsamps)
    :param stats: iterable of test statistics to perform
    :param minsamps: minimum number of samples to perform the test
    :return: 3d numpy array with dimensions
             (njunc, nsamps, len(stats)) containing log P values
    """

    cts = [insamps[0].shape[0], insamps[1].shape[0]]
    # lsv_id, lsv_type = lsv_attrs
    # samps = []
    # cts = []
    # for grp in samples.values():
    #     cts.append(len(grp))
    #     samps.append(grp)
    samps = np.concatenate(insamps, axis=0)
    labels = np.concatenate((np.zeros(cts[0]), np.ones(cts[1])), axis=0).astype(int)
    nfiles, njuncs, nsamples = samps.shape
    outstats = np.ones(shape=(njuncs, len(stats)))
    if all([nsamps >= minsamps for nsamps in cts]):
        for jn_idx in range(njuncs):
            cur_output = np.zeros(shape=(len(stats), nsamples))
            for samp_idx in range(nsamples):
                csamps = samps[:, jn_idx, samp_idx]
                clabels = labels[csamps != -1]
                npos = clabels.sum()
                nneg = clabels.size - npos
                if npos < minsamps or nneg < minsamps:
                    continue
                csamps = csamps[csamps != -1]
                asort = csamps.argsort()
                csamps = csamps[asort]
                clabels = clabels[asort]
                for stat_idx, stat_name in enumerate(stats):
                    cur_output[stat_idx, samp_idx] = operator[stat_name].operator(csamps, clabels)

            outstats[jn_idx, :] = np.percentile(cur_output, 95, axis=1)

    return outstats


def calc_independent(args):
    pipeline_run(independent(args))


class independent(BasicPipeline):

    def store_results(self, output, results, msg_type, extra):
        lsv_id = results[2]
        lsv_type = self.lsv_type_dict[lsv_id]
        mu_psi = results[1][0]
        mean_psi = results[1][1]
        output.heterogen(lsv_id).add(lsv_type=lsv_type, mu_psi=mu_psi, mean_psi=mean_psi, junction_stats=results[0],
                                     junctions=extra['junc_info'][lsv_id])

    def run(self):
        self.independent()

    def independent(self):

        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """

        majiq_logger.create_if_not_exists(self.outDir)
        logger = majiq_logger.get_logger("%s/het_majiq.log" % self.outDir, silent=self.silent,
                                         debug=self.debug)

        logger.info("Majiq deltapsi heterogeneous v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("GROUP1: %s" % self.files1)
        logger.info("GROUP2: %s" % self.files2)

        self.nbins = 40

        manager = mp.Manager()
        self.lsv_type_dict = manager.dict()
        self.queue = mp.Queue()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]
        junc_info = {}

        try:
            for stats_name in self.stats:
                module_ = __import__('majiq.src.stats.' + stats_name.lower(), fromlist=stats_name.title())
                class_ = getattr(module_, stats_name.title())
                operator[stats_name] = class_()
        except ImportError:
            logger.error("The %s statistic is not one of the available statistics, "
                         "in  [ %s ]" % (stats_name, ' | '.join(all_stats)))
            return
        pool = mp.Pool(processes=self.nthreads, initializer=process_conf,
                       initargs=[het_quantification, self],
                       maxtasksperchild=1)
        print(operator)
        list_of_lsv1, exps1 = majiq_io.extract_lsv_summary(self.files1, types_dict=self.lsv_type_dict,
                                                           minnonzero=self.minpos, min_reads=self.minreads,
                                                           junc_info=junc_info, percent=self.min_exp, logger=logger)
        logger.info("Group %s: %s LSVs" % (self.names[0], len(list_of_lsv1)))

        list_of_lsv2, exps2 = majiq_io.extract_lsv_summary(self.files2, types_dict=self.lsv_type_dict,
                                                           minnonzero=self.minpos, min_reads=self.minreads,
                                                           junc_info=junc_info, percent=self.min_exp, logger=logger)
        logger.info("Group %s: %s LSVs" % (self.names[1], len(list_of_lsv1)))

        list_of_lsv = list(set(list_of_lsv1).intersection(set(list_of_lsv2)))
        logger.info("Number quantifiable LSVs: %s" % len(list_of_lsv))

        if len(list_of_lsv) > 0:
            nthreads = min(self.nthreads, len(list_of_lsv))
            [xx.acquire() for xx in self.lock]
            pool.map_async(process_wrapper, chunks(list_of_lsv, nthreads))
            pool.close()
            with Matrix(get_quantifier_voila_filename(self.outDir, self.names, het=True), 'w') as out_h5p:
                out_h5p.file_version = VOILA_FILE_VERSION
                out_h5p.analysis_type = ANALYSIS_HETEROGEN
                out_h5p.group_names = self.names
                out_h5p.experiment_names = [exps1, exps2]
                out_h5p.stat_names = self.stats
                queue_manager(out_h5p, self.lock, self.queue, num_chunks=nthreads, func=self.store_results,
                              junc_info=junc_info, logger=logger, lsv_type=self.lsv_type_dict)
            pool.join()

        if self.mem_profile:
            mem_allocated = int(psutil.Process().memory_info().rss) / (1024 ** 2)
            logger.info("Max Memory used %.2f MB" % mem_allocated)
        logger.info("DeltaPSI Het calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                         self.names[1],
                                                                                                         self.outDir))
