import majiq.src.io_utils
from majiq.src.basic_pipeline import BasicPipeline, _pipeline_run, get_clean_raw_reads
#import pickle
from multiprocessing import Pool, Process
import os
import sys
import gc
import numpy as np
import scipy.misc
from majiq.src.psi import combine_for_priormatrix
from majiq.src.utils.utils import create_if_not_exists, get_logger
import majiq.src.filter as majiq_filter
import majiq.src.io as majiq_io
import majiq.src.psi as majiq_psi
import majiq.src.pipe as pipe
import majiq.src.sample as majiq_sample
import majiq.src.normalize as majiq_norm


def deltapair(args):
    _pipeline_run(DeltaPair(args))


def multi_dpsi(args):
    _pipeline_run(Multi_Deltapair(args))


class DeltaPair(BasicPipeline):
    def run(self):
        self.delta_groups()

    def prepare_lsvs(self, conf, nchunks, logger=None):

        if logger is None:
            logger = get_logger("%s/majiq.log" % self.output, silent=False)
        self.logger = logger
        exec_id = '%s_%s' % (self.names[0], self.names[1])
        # tempfile = '%s/%s_temp_mid_exec.pickle' % (dpsi_obj.output, exec_id)
        num_exp = [len(self.files1), len(self.files2)]

        filtered_lsv1 = [None] * num_exp[0]
        fitfunc = [[None] * num_exp[0], [None] * num_exp[1]]
        meta_info = [[0] * num_exp[0], [0] * num_exp[1]]

        for ii, fname in enumerate(self.files1):
            meta_info[0][ii], lsv_junc, const = majiq_io.load_data_lsv(fname, self.names[0], logger)

            #fitting the function
            #lsv_junc, const = self.gc_content_norm(lsv_junc, const)
            fitfunc[0][ii] = self.fitfunc(const[0])
            filtered_lsv1[ii] = majiq_norm.mark_stacks(lsv_junc, fitfunc[0][ii], self.markstacks, self.logger)
        filtered_lsv1 = majiq_filter.quantifiable_in_group(filtered_lsv1, self.minpos, self.minreads,
                                                           logger=logger)
        logger.info("Group1: %s quantifiable in group" % str(len(filtered_lsv1[0])))

        filtered_lsv2 = [None] * num_exp[1]
        for ii, fname in enumerate(self.files2):
            meta_info[1][ii], lsv_junc, const = majiq_io.load_data_lsv(fname, self.names[1], logger)

            #fitting the function
            #lsv_junc, const = self.gc_content_norm(lsv_junc, const)
            fitfunc[1][ii] = self.fitfunc(const[0])
            filtered_lsv2[ii] = majiq_norm.mark_stacks(lsv_junc, fitfunc[1][ii], self.markstacks, self.logger)
        filtered_lsv2 = majiq_filter.quantifiable_in_group(filtered_lsv2, self.minpos, self.minreads,
                                                           logger=logger)
        logger.info("Group2: %s quantifiable in group" % str(len(filtered_lsv2[0])))

        matched_lsv, matched_info = majiq_filter.lsv_intersection(filtered_lsv1, filtered_lsv2)
        logger.info("After intersection:  %d/(%d, %d)" % (len(matched_info), len(filtered_lsv1[0]),
                                                          len(filtered_lsv2[0])))

        get_clean_raw_reads(matched_info, matched_lsv[0], self.output, self.names[0], num_exp[0])
        get_clean_raw_reads(matched_info, matched_lsv[1], self.output, self.names[1], num_exp[1])

        group1, group2 = combine_for_priormatrix(matched_lsv[0], matched_lsv[1], matched_info, num_exp)
        psi_space, prior_matrix = majiq_psi.gen_prior_matrix(self, group1, group2, self.output, numbins=20,
                                                             defaultprior=self.default_prior)

        outfdir = '%s/tmp/chunks/' % self.output
        if not os.path.exists(outfdir):
            os.makedirs(outfdir)

        logger.info("Saving prior matrix for %s..." % self.names)
        majiq_io.dump_bin_file(prior_matrix, "%s/%s_priormatrix.pickle" % (self.output, exec_id))

        # tout = open("%s/%s_priormatrix.pickle" % (self.output, exec_id), 'w+')
        # pickle.dump(prior_matrix, tout)
        # tout.close()

        logger.info("Saving meta info for %s..." % self.names)
        majiq_io.dump_bin_file(meta_info, "%s/tmp/%s_metainfo.pickle" % (self.output, exec_id))
        # tout = open("%s/tmp/%s_metainfo.pickle" % (self.output, exec_id), 'w+')
        # pickle.dump(meta_info, tout)
        # tout.close()

        csize = len(matched_lsv[0]) / nchunks

        logger.info("Creating %s chunks with <= %s lsv" % (nchunks, csize))
        for nthrd in xrange(nchunks):
            lb = nthrd * csize
            ub = min((nthrd + 1) * csize, len(matched_lsv[0]))
            if nthrd == nchunks - 1:
                ub = len(matched_lsv[0])
            lsv_list = [matched_lsv[0][lb:ub], matched_lsv[1][lb:ub]]
            lsv_info = matched_info[lb:ub]

            out_file = '%s/chunk_%d.pickle' % (outfdir, nthrd)
            majiq_io.dump_bin_file([lsv_list, lsv_info, num_exp, conf, fitfunc, psi_space], out_file)
            # tout = open(out_file, 'w+')
            # pickle.dump([lsv_list, lsv_info, num_exp, conf, fitfunc, psi_space], tout)
            # tout.close()

        gc.collect()

    def delta_groups(self):
        logger = get_logger("%s/majiq.log" % self.logger_path, silent=self.silent, debug=self.debug)
        logger.info("")
        logger.info("Running deltagroups new model ...")
        logger.info("GROUP 1: %s" % self.files1)
        logger.info("GROUP 2: %s" % self.files2)

        exec_id = '%s_%s' % (self.names[0], self.names[1])

        num_exp = [len(self.files1), len(self.files2)]
        nchunks = int(self.nthreads)

        conf = {'minnonzero': self.minpos,
                'minreads': self.minreads,
                'm': self.m,
                'k': self.k,
                'discardzeros': self.discardzeros,
                'trimborder': self.trimborder,
                'debug': self.debug,
                'plotpath': self.plotpath,
                'names': self.names}
        if self.nthreads > 1:
            p = Process(target=self.prepare_lsvs, args=(conf, nchunks))
            p.start()
            p.join()
        else:
            self.prepare_lsvs(conf, nchunks)

        pool = Pool(processes=self.nthreads)
        for nthrd in xrange(nchunks):
            chunk_fname = '%s/tmp/chunks/chunk_%d.pickle' % (self.output, nthrd)
            delta_prior_path = "%s/%s_priormatrix.pickle" % (self.output, exec_id)
            if self.nthreads == 1:
                pipe.parallel_lsv_child_calculation(deltapsi_quantify,
                                                    [chunk_fname, delta_prior_path, True],
                                                    '%s/tmp' % self.output,
                                                    '%s_%s' % (self.names[0], self.names[1]),
                                                    nthrd)

            else:
                pool.apply_async(pipe.parallel_lsv_child_calculation, [deltapsi_quantify,
                                                                       [chunk_fname, delta_prior_path, True],
                                                                       '%s/tmp' % self.output,
                                                                       '%s_%s' % (self.names[0], self.names[1]),
                                                                       nthrd])

        if self.nthreads > 1:
            pool.close()
            pool.join()

        posterior_matrix = []
        names = []
        psi_list1 = []
        psi_list2 = []
        logger.info("GATHER pickles")
        for nthrd in xrange(self.nthreads):

            ptempt = majiq.src.io_utils.load_bin_file("%s/tmp/%s_%s_th%s.%s.pickle" % (self.output, self.names[0], self.names[1],
                                                                                       nthrd, deltapsi_quantify.__name__))
            # tempfile = open("%s/tmp/%s_%s_th%s.%s.pickle" % (self.output, self.names[0],
            #                                                  self.names[1], nthrd, deltapsi_quantify.__name__))
            # ptempt = pickle.load(tempfile)
            posterior_matrix.extend(ptempt[0])
            names.extend(ptempt[1])
            psi_list1.extend(ptempt[2])
            psi_list2.extend(ptempt[3])

        logger.debug("Getting meta info for %s..." % self.names)
        meta_info = majiq.src.io_utils.load_bin_file("%s/tmp/%s_metainfo.pickle" % (self.output, exec_id))
        # tin = open("%s/tmp/%s_metainfo.pickle" % (self.output, exec_id))
        # meta_info = pickle.load(tin)
        # tin.close()

        pickle_path = "%s/%s_%s.%s.pickle" % (self.output, self.names[0], self.names[1], deltapsi_quantify.__name__)
        # pickle.dump([posterior_matrix, names, meta_info, psi_list1, psi_list2], open(pickle_path, 'w'))
        majiq_io.dump_lsvs_voila(pickle_path, posterior_matrix, names, meta_info, psi_list1, psi_list2)
        logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                     self.names[1],
                                                                                                     self.output))

        logger.info("Alakazam! Done.")


class Multi_Deltapair(BasicPipeline):
    def run(self):
        self.multi_dpsi()

    def _auxiliar_multi(self, fname, name, conf, logger=None):
        if logger is None:
            logger = get_logger("%s/majiq.log" % self.output, silent=False)
            self.logger = logger

        meta_info, lsv_junc, const = majiq_io.load_data_lsv(fname, name, logger)

        # fitting the function
        # lsv_junc, const = self.gc_content_norm(lsv_junc, const)
        fitfunc = self.fitfunc(const[0])
        filtered_lsv = majiq_norm.mark_stacks(lsv_junc, fitfunc, self.markstacks, self.logger)

        lsv_info = np.zeros(shape=(len(filtered_lsv[0])), dtype=np.dtype('object'))
        lsv_samples = np.zeros(shape=(len(filtered_lsv[0])), dtype=np.dtype('object'))
        filt_vals = np.zeros(shape=(len(filtered_lsv[0])), dtype=np.dtype('bool'))

        for lidx, lsv in enumerate(filtered_lsv[0]):
            for jj in lsv:
                filt_vals[lidx] |= (np.count_nonzero(jj) >= conf['minnonzero'] and jj.sum() >= conf['minreads'])

            m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(junction_list=filtered_lsv[0][lidx],
                                                                       m=conf['m'],
                                                                       k=conf['k'],
                                                                       discardzeros=conf['discardzeros'],
                                                                       trimborder=conf['trimborder'],
                                                                       fitted_one_over_r=fitfunc,
                                                                       debug=conf['debug'])

            lsv_samples[lidx] = s_lsv
            lsv_info[lidx] = filtered_lsv[1][lidx]

        outfdir = '%s/tmp/samples/' % self.output
        if not os.path.exists(outfdir):
            os.makedirs(outfdir)
        out_file = '%s/%s.pickle' % (outfdir, name)
        majiq_io.dump_bin_file([meta_info, (lsv_samples, lsv_info), filt_vals], out_file)
        # tout = open(out_file, 'w+')
        # pickle.dump([meta_info, (lsv_samples, lsv_info), filt_vals], tout)
        # tout.close()

    def multi_dpsi(self):

        logger = get_logger("%s/majiq.log" % self.logger_path, silent=self.silent, debug=self.debug)
        logger.info("")
        logger.info("Running multi_deltas...")
        logger.info("Config file: %s" % self.deltapairs)
        groups, dict_files, list_deltas = majiq_io.read_multi_dpsi_conf(self.deltapairs)
        print list_deltas
        # exec_id = '%s_%s' % (self.names[0], self.names[1])

        #num_exp = [len(self.files1), len(self.files2)]
        nchunks = int(self.nthreads)

        conf = {'minnonzero': self.minpos,
                'minreads': self.minreads,
                'm': self.m,
                'k': self.k,
                'discardzeros': self.discardzeros,
                'trimborder': self.trimborder,
                'debug': self.debug,
                'plotpath': self.plotpath}

        outfdir = '%s/tmp/samples/' % self.output
        if not os.path.exists(outfdir):
            for name, fname in dict_files.items():

                if self.nthreads > 1:
                    p = Process(target=self._auxiliar_multi, args=(fname, name, conf))
                    p.start()
                    p.join()
                else:
                    self._auxiliar_multi(fname, name, conf)

        pool = Pool(processes=self.nthreads)

        onlygather = True

        for group1, group2 in list_deltas:
            dpsi_name = '%s_%s' % (group1, group2)

            num_exp = [len(groups[group1]), len(groups[group2])]
            meta_info = [[0] * num_exp[0], [0] * num_exp[1]]
            self.names = [group1, group2]
            self.logger = logger
            if not onlygather:

                matched_files = [None] * num_exp[0]
                for idx, ii in enumerate(groups[group1]):
                    outfdir = '%s/tmp/samples/' % self.output
                    infile = '%s/%s.pickle' % (outfdir, ii)
                    meta_info[0][idx], matched_files[idx], filt_vals = majiq.src.io_utils.load_bin_file(infile)
                    #pickle.load(open(infile))
                filtered_lsv1 = majiq_filter.quantifiable_in_group(matched_files, self.minpos, self.minreads,
                                                                   filt_vals, logger)

                matched_files = [None] * num_exp[1]
                for idx, ii in enumerate(groups[group2]):
                    outfdir = '%s/tmp/samples/' % self.output
                    infile = '%s/%s.pickle' % (outfdir, ii)
                    meta_info[1][idx], matched_files[idx], filt_vals = majiq.src.io_utils.load_bin_file(infile)
                    #pickle.load(open(infile))
                filtered_lsv2 = majiq_filter.quantifiable_in_group(matched_files, self.minpos, self.minreads,
                                                                   filt_vals, logger)

                matched_lsv, matched_info = majiq_filter.lsv_intersection(filtered_lsv1, filtered_lsv2, bycol=True)
                logger.info("After intersection:  %d/(%d, %d)" % (len(matched_info), len(filtered_lsv1[0]),
                                                                  len(filtered_lsv2[0])))


                group1, group2 = combine_for_priormatrix(matched_lsv[0], matched_lsv[1], matched_info, num_exp)
                psi_space, prior_matrix = majiq_psi.gen_prior_matrix(self, group1, group2, self.output, numbins=20,
                                                                     defaultprior=self.default_prior)


                outfdir = '%s/tmp/%s/chunks/' % (self.output, dpsi_name)
                if not os.path.exists(outfdir):
                    os.makedirs(outfdir)

                logger.debug("Saving prior matrix for %s..." % dpsi_name)
                dpsi_prior_name = "%s/%s_priormatrix.pickle" % (self.output, dpsi_name)
                majiq_io.dump_bin_file(prior_matrix, dpsi_prior_name)
                # tout = open(dpsi_prior_name, 'w+')
                # pickle.dump(prior_matrix, tout)
                # tout.close()

                logger.debug("Saving meta info for %s..." % dpsi_name)
                majiq_io.dump_bin_file(meta_info, "%s/tmp/%s_metainfo.pickle" % (self.output, dpsi_name))
                # tout = open("%s/tmp/%s_metainfo.pickle" % (self.output, dpsi_name), 'w+')
                # pickle.dump(meta_info, tout)
                # tout.close()

                csize = len(matched_lsv[0]) / nchunks

                logger.info("Creating %s chunks with <= %s lsv" % (nchunks, csize))
                for nthrd in xrange(nchunks):
                    lb = nthrd * csize
                    ub = min((nthrd + 1) * csize, len(matched_lsv[0]))
                    if nthrd == nchunks - 1:
                        ub = len(matched_lsv[0])
                    lsv_list = [matched_lsv[0][lb:ub], matched_lsv[1][lb:ub]]
                    lsv_info = matched_info[lb:ub]

                    out_file = '%s/chunk_%d.pickle' % (outfdir, nthrd)
                    majiq_io.dump_bin_file([lsv_list, lsv_info, num_exp, conf, None, psi_space], out_file)
                    # tout = open(out_file, 'w+')
                    # pickle.dump([lsv_list, lsv_info, num_exp, conf, None, psi_space], tout)
                    # tout.close()

                    if self.nthreads == 1:
                        pipe.parallel_lsv_child_calculation(deltapsi_quantify,
                                                            [out_file, dpsi_prior_name, False],
                                                            outfdir,
                                                            dpsi_name,
                                                            nthrd)

                    else:
                        pool.apply_async(pipe.parallel_lsv_child_calculation, [deltapsi_quantify,
                                                                               [out_file, dpsi_prior_name, False],
                                                                               outfdir,
                                                                               dpsi_name,
                                                                               nthrd])

                if self.nthreads > 1:
                    pool.close()
                    pool.join()

                    gc.collect()

            posterior_matrix = []
            names = []
            psi_list1 = []
            psi_list2 = []
            logger.info("GATHER pickles")
            for nthrd in xrange(self.nthreads):
                ptempt = majiq.src.io_utils.load_bin_file("%s/tmp/%s/chunks/%s_%s_th%s.%s.pickle" % (self.output, dpsi_name,
                                                                                                     self.names[0], self.names[1],
                                                                                                     nthrd,
                                                                                                     deltapsi_quantify.__name__))
                # tempfile = open("%s/tmp/%s/chunks/%s_%s_th%s.%s.pickle" % (self.output, dpsi_name, self.names[0],
                #                                                  self.names[1], nthrd, deltapsi_quantify.__name__))
                # ptempt = pickle.load(tempfile)
                posterior_matrix.extend(ptempt[0])
                names.extend(ptempt[1])
                psi_list1.extend(ptempt[2])
                psi_list2.extend(ptempt[3])

            logger.debug("Getting meta info for %s..." % self.names)
            meta_info = majiq.src.io_utils.load_bin_file("%s/tmp/%s_metainfo.pickle" % (self.output, dpsi_name))

            # tin = open("%s/tmp/%s_metainfo.pickle" % (self.output, dpsi_name))
            # pickle.load(tin)
            # tin.close()

            pickle_path = "%s/%s_%s.%s.pickle" % (self.output, self.names[0], self.names[1], deltapsi_quantify.__name__)
            # pickle.dump([posterior_matrix, names, meta_info, psi_list1, psi_list2], open(pickle_path, 'w'))
            majiq_io.dump_lsvs_voila(pickle_path, posterior_matrix, names, meta_info, psi_list1, psi_list2)
            logger.info("DeltaPSI calculation for %s_%s ended succesfully! Result can be found at %s" % (self.names[0],
                                                                                                         self.names[1],
                                                                                                         self.output))

        logger.info("Alakazam! Done.")


def deltapsi_quantify(fname, delta_prior_path, boots_sample=True, logger=None):

    matched_lsv, info, num_exp, conf, fitfunc, psi_space, prior_matrix = pipe.__load_execution_chunk(fname,
                                                                                                     delta_prior_path)

    if boots_sample:

        lsv_samples1 = np.zeros(shape=(len(info), num_exp[0]), dtype=np.dtype('object'))
        lsv_samples2 = np.zeros(shape=(len(info), num_exp[1]), dtype=np.dtype('object'))

        logger.info("Bootstrapping for all samples...")
        for grp_idx, group in enumerate(matched_lsv):
            for lidx, lsv_all in enumerate(group):
                for eidx, lsv in enumerate(lsv_all):
                    m_lsv, var_lsv, s_lsv = majiq_sample.sample_from_junctions(junction_list=lsv,
                                                                               m=conf['m'],
                                                                               k=conf['k'],
                                                                               discardzeros=conf['discardzeros'],
                                                                               trimborder=conf['trimborder'],
                                                                               fitted_one_over_r=fitfunc[grp_idx][eidx],
                                                                               debug=conf['debug'])
                    if grp_idx == 0:
                        lsv_samples1[lidx, eidx] = s_lsv
                    else:
                        lsv_samples2[lidx, eidx] = s_lsv

    else:
        lsv_samples1 = np.array(matched_lsv[0])
        lsv_samples2 = np.array(matched_lsv[1])

    nbins = len(psi_space)
    logger.info("Calculating deltas...")

    post_matrix = []
    new_info = []
    ones_n = np.ones(shape=(1, nbins), dtype=np.float)

    # pickle.dump([lsv_samples1, info], open('./lsv_binomproblem.pkl', 'w+b'))
    posterior_psi1 = []
    posterior_psi2 = []
    #print lsv_samples1
    for lidx, lsv_info in enumerate(info):
        num_ways = len(lsv_samples1[lidx][0])
        if lidx % 50 == 0:
            print "Event %d ..." % lidx,
            sys.stdout.flush()

        alpha_prior, beta_prior = majiq_psi.__get_prior_params(lsv_info, num_ways)
        if 'i' in lsv_info[2]:
            prior_idx = 1
        else:
            prior_idx = 0
        post_matrix.append([])
        new_info.append(lsv_info)

        posterior_psi1.append([])
        posterior_psi2.append([])
        psi1 = lsv_samples1[lidx, :]
        psi2 = lsv_samples2[lidx, :]

        for p_idx in xrange(num_ways):

            posterior = np.zeros(shape=(nbins, nbins), dtype=np.float)
            post_psi1 = np.zeros(shape=nbins, dtype=np.float)
            post_psi2 = np.zeros(shape=nbins, dtype=np.float)
            for m in xrange(conf['m']):
                # log(p(D_T1(m) | psi_T1)) = SUM_t1 T ( log ( P( D_t1 (m) | psi _T1)))
                junc = [psi1[xx][p_idx][m] for xx in xrange(num_exp[0])]
                junc = np.array(junc)
                all_sample = [psi1[xx][yy][m].sum() for xx in xrange(num_exp[0]) for yy in xrange(num_ways)]
                all_sample = np.array(all_sample)
                data_given_psi1 = np.log(majiq_psi.prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                    alpha_prior[p_idx], beta_prior[p_idx]))

                psi_v1 = data_given_psi1.reshape(nbins, -1)
                post_psi1 += (data_given_psi1 - scipy.misc.logsumexp(data_given_psi1))

                junc = [psi2[xx][p_idx][m] for xx in xrange(num_exp[1])]
                junc = np.array(junc)
                all_sample = [psi2[xx][yy][m].sum() for xx in xrange(num_exp[1]) for yy in xrange(num_ways)]
                all_sample = np.array(all_sample)
                data_given_psi2 = np.log(majiq_psi.prob_data_sample_given_psi(junc.sum(), all_sample.sum(), nbins,
                                                                    alpha_prior[p_idx], beta_prior[p_idx]))
                post_psi2 += (data_given_psi2 - scipy.misc.logsumexp(data_given_psi2))
                psi_v2 = data_given_psi2.reshape(-1, nbins)

                A = (psi_v1 * ones_n + psi_v2 * ones_n.T) + np.log(prior_matrix[prior_idx])

                posterior += np.exp(A - scipy.misc.logsumexp(A))

            post_matrix[-1].append(posterior / conf['m'])
            posterior_psi1[-1].append(post_psi1 / conf['m'])
            posterior_psi2[-1].append(post_psi2 / conf['m'])
            if num_ways == 2:
                break

    return post_matrix, new_info, posterior_psi1, posterior_psi2
