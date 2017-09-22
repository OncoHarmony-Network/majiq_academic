import multiprocessing as mp
import sys
import traceback

import h5py

import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.utils as majiq_utils
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.constants import *
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper
from majiq.src.psi import heterogen_posterior
from majiq.src.stats import operator, all_stats
from voila.api import Voila
from voila.constants import ANALYSIS_HETEROGEN
from voila.vlsv import Het


def het_quantification(list_of_lsv, chnk, process_conf, logger):
    logger.info("Quantifying LSVs PSI.. %s" % chnk)
    num_exp = [len(process_conf.files1), len(process_conf.files2)]

    f_list = [None, None]
    f_list[0] = majiq_io.get_extract_lsv_list(list_of_lsv, process_conf.files1)
    f_list[1] = majiq_io.get_extract_lsv_list(list_of_lsv, process_conf.files2)

    for lidx, lsv_id in enumerate(list_of_lsv):
        if lidx % 50 == 0:
            print("Event %d ..." % lidx)
            sys.stdout.flush()

        boots = [np.array(f_list[0][lidx].coverage), np.array(f_list[1][lidx].coverage)]
        lsv_type = [f_list[0][lidx].type, f_list[1][lidx].type]
        lsv_het = Het()
        msamples = boots[0].shape[2]
        assert boots[0].shape[2] == boots[1].shape[2], "LSV %s, has different types in %s and %s (%s vs %s). " \
                                                       "Please check that the conditions has been build together." \
                                                       % (lsv_id, process_conf.names[0], process_conf.names[1],
                                                          boots[0].shape[2], boots[1].shape[2])
        samps = heterogen_posterior(boots, lsv_het, msamples, process_conf.nsamples, process_conf.vwindow, num_exp,
                                    process_conf.nbins, lsv_type[0])

        out_stats = do_test_stats(samps, process_conf.stats, process_conf.minsamps)
        for stat_idx in range(out_stats.shape[1]):
            lsv_het.add_junction_stats(out_stats[:, stat_idx])

            # if num_ways == 2:
            #     break

        qm = QueueMessage(QUEUE_MESSAGE_HETER_DELTAPSI, (lsv_het, lsv_id), chnk)
        process_conf.queue.put(qm, block=True)


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
    outstats = np.zeros(shape=(njuncs, len(stats)))
    if all([nsamps > minsamps for nsamps in cts]):
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
    def run(self):
        self.independent()

    def independent(self):
        majiq_utils.create_if_not_exists(self.logger_path)
        logger = majiq_utils.get_logger("%s/deltapsi_majiq.log" % self.logger_path, silent=self.silent,
                                        debug=self.debug)

        logger.info("Majiq deltapsi v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))
        logger.info("GROUP1: %s" % self.files1)
        logger.info("GROUP2: %s" % self.files2)
        self.nbins = 40

        try:
            for stats_name in self.stats:
                module_ = __import__('majiq.src.stats.' + stats_name.lower(), fromlist=stats_name.title())
                class_ = getattr(module_, stats_name.title())
                operator[stats_name] = class_()
        except ImportError as i_err:
            logger.error("The %s statistic is not one of the available statistics, "
                         "in  [ %s ]" % (stats_name, ' | '.join(all_stats)))
            return

        lsv_dict1, lsv_types1, lsv_summarized1, meta1, lsv_dict_graph1 = majiq_io.extract_lsv_summary(self.files1)
        list_of_lsv_group1 = majiq_filter.merge_files_hdf5(lsv_dict1, lsv_summarized1, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=logger)

        lsv_dict2, lsv_types2, lsv_summarized2, meta2, lsv_dict_graph2 = majiq_io.extract_lsv_summary(self.files2)
        list_of_lsv_group2 = majiq_filter.merge_files_hdf5(lsv_dict2, lsv_summarized2, self.minpos,
                                                           self.minreads, percent=self.min_exp, logger=logger)

        list_of_lsv = list(set(list_of_lsv_group1).intersection(set(list_of_lsv_group2)))
        lchnksize = max(len(list_of_lsv) / self.nthreads, 1) + 1
        if len(list_of_lsv) > 0:
            self.queue = mp.Queue()
            self.lock = [mp.Lock() for xx in range(self.nthreads)]
            [xx.acquire() for xx in self.lock]
            pool = mp.Pool(processes=self.nthreads, initializer=process_conf,
                           initargs=[het_quantification, self],
                           maxtasksperchild=1)

            pool.map_async(process_wrapper, majiq_utils.chunks2(list_of_lsv, lchnksize, extra=range(self.nthreads)))
            pool.close()
            with Voila(get_quantifier_voila_filename(self.outDir, self.names, deltapsi=True), 'w') as out_h5p:
                # out_h5p.add_genome(meta1['genome'])
                out_h5p.set_analysis_type(ANALYSIS_HETEROGEN)
                out_h5p.add_experiments(self.names[0], experiment_names=meta1['experiments'])
                out_h5p.add_experiments(self.names[1], experiment_names=meta2['experiments'])
                out_h5p.add_stat_names(self.stats)

                in_h5p = h5py.File(self.files1[0], 'r')
                queue_manager(in_h5p, out_h5p, self.lock, self.queue, num_chunks=self.nthreads, logger=logger,
                              list_of_lsv_graphics=lsv_dict_graph1)
                in_h5p.close()
            pool.join()
