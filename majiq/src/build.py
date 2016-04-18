import os
import multiprocessing as mp
from majiq.src.basic_pipeline import BasicPipeline, _pipeline_run
import majiq.src.config as majiq_config
import majiq.utils as majiq_utils
import h5py
def build(args):
    _pipeline_run(Builder(args))


class Builder(BasicPipeline):
    def run(self):
        if not os.path.exists(self.conf):
            raise RuntimeError("Config file %s does not exist" % self.conf)
        majiq_config.global_conf_ini(self.conf, self)
        self.builder()

    def builder(self):

        logger = majiq_utils.get_logger("%s/majiq.log" % majiq_config.outDir, silent=majiq_utils.silent,
                                        debug=majiq_utils.debug)
        logger.info("")
        logger.info("Command: %s" % self)

        import yappi
        yappi.start()

        p = majiq_multi.parallel_lsv_child_calculation(majiq_io.read_gff, [self.transcripts, self.pcr_filename,
                                                                     self.nthreads, self.tlogger])

        logger.debug("... waiting gff3 parsing")
        p.start()
        p.join()

        logger.debug("Get samfiles")
        sam_list = []
        for exp_idx, exp in enumerate(majiq_config.exp_list):
            samfile = "%s/%s.bam" % (majiq_config.sam_dir, exp)
            if not os.path.exists(samfile):
                raise RuntimeError("Skipping %s.... not found" % samfile)
            baifile = "%s/%s.bam.bai" % (majiq_config.sam_dir, exp)
            if not os.path.exists(baifile):
                raise RuntimeError("Skipping %s.... not found ( index file for bam file is required)" % baifile)
            sam_list.append(samfile)
            majiq_config.exp_list[exp_idx] = os.path.split(exp)[1]

        if self.nthreads > 1:
            lock_array = [mp.Lock() for xx in range(majiq_config.num_final_chunks)]
            q = mp.Queue()
            pool = mp.Pool(processes=self.nthreads, initializer=__builder_init,
                           initargs=[q, lock_array], maxtasksperchild=1)

        for chnk in range(majiq_config.num_final_chunks):
            temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
            majiq_utils.create_if_not_exists(temp_dir)
            if self.nthreads == 1:
                majiq_builder(sam_list, chnk, pcr_validation=self.pcr_filename, gff_output=self.gff_output,
                              only_rna=self.only_rna, nondenovo=self.non_denovo, logging=logger)
            else:
                lock_array[chnk].acquire()
                pool.apply_async(majiq_multi.parallel_lsv_child_calculation, [majiq_builder, [sam_list, chnk,
                                                                                              self.pcr_filename,
                                                                                              self.gff_output,
                                                                                              self.only_rna,
                                                                                              self.non_denovo,
                                                                                              self.silent,
                                                                                              self.debug]
                                                                                             ])

        count = 0

        lsv_list = []
        junc_list = []
        splicegraph = []
        file_list = []
        junc_idx = [0] * majiq_config.num_experiments
        for exp_idx, exp in enumerate(majiq_config.exp_list):
            fname = "%s/%s.majiq.hdf5" % (majiq_config.outDir, majiq_config.exp_list[exp_idx])
            f = h5py.File(fname, 'w', compression='gzip', compression_opts=9)
            file_list.append(f)

            fname_sg = "%s/%s.splicegraph.hdf5" % (majiq_config.outDir, majiq_config.exp_list[exp_idx])
            f_splicegraph = h5py.File(fname_sg, 'w', compression='gzip', compression_opts=9)

            #as_table = f.create_group('LSVs')
            lsv_list.append(f)
            effective_readlen = (majiq_config.readLen - 16) + 1
            non_as_table = f.create_dataset("/const_junctions",
                                            (majiq_config.nrandom_junctions, effective_readlen),
                                            maxshape=(None, effective_readlen))
            junc_list.append(non_as_table)

            splicegraph.append(f_splicegraph)

        #logger.debug("... waiting childs")
        pool.close()
        while True:
            try:
                val = q.get(block=True, timeout=10)
                if val[0] == 0:
                    for exp_idx in majiq_config.tissue_repl[val[2]]:
                        val[1].to_hdf5(lsv_list[exp_idx])
                elif val[0] == 1:
                    for jdx, exp_idx in enumerate(majiq_config.tissue_repl[val[2]]):
                        if junc_idx[exp_idx] >= majiq_config.nrandom_junctions:
                            continue
                        junc_list[exp_idx][junc_idx[exp_idx], :] = val[1][jdx, :].toarray()
                        junc_idx[exp_idx] += 1

                elif val[0] == 2:
                    val[1].to_hdf5(splicegraph[val[2]])
                    pass
                elif val[0] == 3:
                    lock_array[val[1]].release()
                    count += 1

            except Queue.Empty:
                if count < majiq_config.num_final_chunks:
                    continue
                break
        pool.join()

        for ff in file_list:
            ff.close()

        for ff in splicegraph:
            ff.close()


        #GATHER
    #    gather_files(majiq_config.outDir, '', params.gff_output, params.pcr_filename,
    #                       nthreads=params.nthreads, logger=logger)

        yappi.get_func_stats().print_all()

        majiq_config.print_numbers()
        logger.info("End of execution")

def majiq_builder(samfiles_list, chnk, pcr_validation=None, gff_output=None, create_tlb=True, only_rna=False,
                  nondenovo=False, logging=None):

    if not logging is None:
        logging.debug("Building chunk %s" % chnk)

    temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
    annot_file = '%s/annot_genes.pkl' % temp_dir
    if not os.path.exists(annot_file):
        majiq_builder.queue.put([3, chnk], block=True)
        majiq_builder.lock_arr[chnk].acquire()
        majiq_builder.lock_arr[chnk].release()
        return

    # temp_file = open(annot_file, 'rb')
    gene_list = majiq.src.io_utils.load_bin_file(annot_file)

    if create_tlb:
        if not logging is None:
            logging.debug("[%s] Recreatin Gene TLB" % chnk)
        majiq_gene.recreate_gene_tlb(gene_list)

    if not logging is None:
        logging.info("[%s] Reading BAM files" % chnk)
    majiq_io.read_sam_or_bam(samfiles_list, gene_list, chnk,
                             nondenovo=nondenovo, logging=logging)
    if not logging is None:
        logging.info("[%s] Detecting intron retention events" % chnk)
    majiq_io.rnaseq_intron_retention(samfiles_list, gene_list, chnk, permissive=majiq_config.permissive_ir,
                                     nondenovo=nondenovo, logging=logging)
    if not logging is None:
        logging.info("[%s] Detecting LSV" % chnk)
    const = analize.lsv_detection(gene_list, chnk, only_real_data=only_rna, out_queue=majiq_builder.queue,
                                  logging=logging)


    prepare_gc_content(gene_list, temp_dir)

    # if pcr_validation:
    #     utils.get_validated_pcr_lsv(lsv, temp_dir)
    # if gff_output:
    #     majiq_lsv.extract_gff(lsv, temp_dir)
    majiq.src.voila_wrapper.generate_visualization_output(gene_list, temp_dir, majiq_builder.queue)
    if not logging is None:
        logging.info("[%s] Preparing output" % chnk)
    lsv = None
    utils.send_output(lsv, const, temp_dir, majiq_builder.queue, chnk, majiq_builder.lock_arr[chnk])


def __builder_init(out_queue, lock_arr):
    majiq_builder.queue = out_queue
    majiq_builder.lock_arr = lock_arr


