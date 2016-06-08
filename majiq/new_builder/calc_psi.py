
import os
import majiq.src.utils.utils as majiq_utils
from majiq.src.basic_pipeline import BasicPipeline, _pipeline_run
import majiq_io as majiq_io
import majiq.src.normalize as majiq_norm
import majiq.src.filter as majiq_filter

################################
# PSI calculation pipeline     #
################################

def calcpsi(args):
    return _pipeline_run(CalcPsi(args))


class CalcPsi(BasicPipeline):

    def run(self):
        self.calcpsi()

    def pre_psi(self, nchunks, logger=None):

        if logger is None:
            logger = majiq_utils.get_logger("%s/majiq.log" % self.output, silent=False)

        self.logger = logger

        num_exp = len(self.files)
        meta_info = [0] * num_exp
        filtered_lsv = [None] * num_exp
        fitfunc = [None] * num_exp
        for ii, fname in enumerate(self.files):
            meta_info[ii], lsv_junc, const = majiq_io.load_data_lsv(fname, self.name, logger)
            fitfunc[ii] = self.fitfunc(const[0])
            filtered_lsv[ii] = majiq_norm.mark_stacks(lsv_junc, fitfunc[ii], self.markstacks, self.logger)

        matched_lsv, matched_info = majiq_filter.quantifiable_in_group(filtered_lsv, self.minpos, self.minreads,
                                                                       logger=logger)

        #get_clean_raw_reads(matched_info, matched_lsv, self.output, self.name, num_exp)

        csize = len(matched_lsv) / nchunks
        outfdir = '%s/tmp/chunks/' % self.output
        if not os.path.exists(outfdir):
            os.makedirs(outfdir)

        logger.debug("Saving meta info for %s..." % self.name)
        majiq_io.dump_bin_file(meta_info, "%s/tmp/%s_metainfo.pickle" % (self.output, self.name))


        logger.info("Creating %s chunks with <= %s lsv" % (nchunks, csize))
        for nthrd in xrange(nchunks):
            lb = nthrd * csize
            ub = min((nthrd + 1) * csize, len(matched_lsv))
            if nthrd == nchunks - 1:
                ub = len(matched_lsv)
            lsv_list = matched_lsv[lb:ub]
            lsv_info = matched_info[lb:ub]

            out_file = '%s/chunk_%d.pickle' % (outfdir, nthrd)
            majiq_io.dump_bin_file([lsv_list, lsv_info, num_exp, fitfunc], out_file)

        gc.collect()

    def calcpsi(self):
        """
        Given a file path with the junctions, return psi distributions.
        write_pickle indicates if a .pickle should be saved in disk
        """

        logger = majiq_utils.get_logger("%s/majiq.log" % self.logger_path, silent=self.silent, debug=self.debug)
        logger.info("")
        logger.info("Running Psi ...")
        logger.info("GROUP: %s" % self.files)



        conf = {'minnonzero': self.minpos,
                'minreads': self.minreads,
                'm': self.m,
                'k': self.k,
                'discardzeros': self.discardzeros,
                'trimborder': self.trimborder,
                'debug': self.debug,
                'nbins': 40,
                'name': self.name}

        nchunks = self.nthreads
        if self.nthreads > 1:
            p = Process(target=self.pre_psi, args=[nchunks])
            p.start()
            p.join()

            pool = Pool(processes=self.nthreads)
        else:
            self.pre_psi(nchunks)

        for nthrd in xrange(nchunks):
            chunk_fname = '%s/tmp/chunks/chunk_%d.pickle' % (self.output, nthrd)
            if self.nthreads == 1:
                pipe.parallel_lsv_child_calculation(pipe.calcpsi,
                                                    [chunk_fname, conf],
                                                    '%s/tmp' % self.output,
                                                    self.name,
                                                    nthrd)

            else:
                pool.apply_async(pipe.parallel_lsv_child_calculation, [pipe.calcpsi,
                                                                       [chunk_fname, conf],
                                                                       '%s/tmp' % self.output,
                                                                       self.name,
                                                                       nthrd])

        if self.nthreads > 1:
            pool.close()
            pool.join()

        posterior_matrix = []
        names = []
        logger.info("GATHER pickles")
        for nt in xrange(self.nthreads):
            tempfile = open("%s/tmp/%s_th%s.calcpsi.pickle" % (self.output, self.name, nt))
            ptempt = majiq.src.io_utils.load_bin_file(tempfile)
            posterior_matrix.extend(ptempt[0])
            names.extend(ptempt[1])

        logger.debug("Getting meta info for %s..." % self.name)
        tin = open("%s/tmp/%s_metainfo.pickle" % (self.output, self.name))
        meta_info = majiq.src.io_utils.load_bin_file(tin)
        tin.close()


        pickle_path = "%s/%s_psigroup.pickle" % (self.output, self.name)
        majiq_io.dump_lsvs_voila(pickle_path, posterior_matrix, names, meta_info)
        # pickle.dump([posterior_matrix, names, meta_info], open(pickle_path, 'w'))
        logger.info("PSI calculation for %s ended succesfully! Result can be found at %s" % (self.name, self.output))
        logger.info("Alakazam! Done.")
