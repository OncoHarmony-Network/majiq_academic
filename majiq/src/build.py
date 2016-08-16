import Queue
import multiprocessing as mp
import os
import sys
import traceback
import types
import gc
import h5py

import io as majiq_io
import majiq.grimoire.gene
import majiq.src.config as majiq_config
import majiq.src.utils as majiq_utils
import multiproc as majiq_multi
import majiq.src.normalize as majiq_norm
from analize import lsv_detection
from constants import *
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run


def build(args):
    pipeline_run(Builder(args))


def builder_init(out_queue, lock_arr, sam_list, pcr_filename, gff_output, only_rna,
                 non_denovo, dbfile, silent, debug):
    majiq_builder.queue = out_queue
    majiq_builder.lock_arr = lock_arr
    majiq_builder.sam_list = sam_list
    majiq_builder.pcr_filename = pcr_filename
    majiq_builder.gff_output = gff_output
    majiq_builder.only_rna = only_rna
    majiq_builder.non_denovo = non_denovo
    majiq_builder.dbfile = dbfile
    majiq_builder.silent = silent
    majiq_builder.debug = debug


def majiq_builder(args_vals):

    list_of_genes, chnk = args_vals
    tlogger = majiq_utils.get_logger("%s/%s.majiq.log" % (majiq_config.outDir, chnk),
                                     silent=majiq_builder.silent, debug=majiq_builder.debug)

    tlogger.debug("[%s] Starting new chunk" % chnk)
    majiq_utils.monitor('CHILD %s:: CREATION' % chnk)

    db_f = h5py.File(majiq_builder.dbfile)
    if isinstance(list_of_genes, types.StringTypes):
        list_of_genes = [list_of_genes]

    counter = [0] * 6
    # samfile = [pysam.Samfile(xx, "rb") for xx in majiq_builder.sam_list]
    for gne_id in list_of_genes:
        majiq_utils.monitor('CHILD %s:: STARTLOOP' % chnk)

        try:
            loop_id = '%s - %s' % (chnk, gne_id)
            tlogger.info("[%s] Retrieving gene" % loop_id)
            gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f)

            tlogger.info("[%s] Reading BAM files" % loop_id)
            majiq_io.read_sam_or_bam(gene_obj, majiq_builder.sam_list, counter, h5py_file=db_f,
                                     nondenovo=majiq_builder.non_denovo, info_msg=loop_id, logging=tlogger)

            if gene_obj.get_read_count() == 0:
                continue

            tlogger.info("[%s] Detecting intron retention events" % loop_id)
            majiq_io.rnaseq_intron_retention(gene_obj, majiq_builder.sam_list, chnk,
                                             permissive=majiq_config.permissive_ir,
                                             nondenovo=majiq_builder.non_denovo, logging=tlogger)

            tlogger.info("[%s] Detecting LSV" % loop_id)
            lsv_detection(gene_obj, only_real_data=majiq_builder.only_rna,
                          out_queue=majiq_builder.queue, logging=tlogger)

            majiq_utils.monitor('CHILD %s:: ENDLOOP' % chnk)

        except Exception as e:
            majiq_utils.monitor('CHILD %s:: EXCEPT' % chnk)
            traceback.print_exc()
            sys.stdout.flush()
            raise

        finally:
            del gene_obj
            del majiq_config.gene_tlb[gne_id]
            gc.collect()

    db_f.close()
    majiq_utils.monitor('CHILD %s:: WAITING' % chnk)
    tlogger.info("[%s] Waiting to be freed" % chnk)
    majiq_builder.queue.put([-1, chnk], block=True)
    majiq_builder.lock_arr[chnk].acquire()
    majiq_builder.lock_arr[chnk].release()


class Builder(BasicPipeline):

    def run(self):
        if not os.path.exists(self.conf):
            raise RuntimeError("Config file %s does not exist" % self.conf)
        sam_list = majiq_config.global_conf_ini(self.conf, self)
        self.builder(sam_list)

    def queue_manager(self, lock_array, result_queue, first_id=0, vfunc_gc=None, logger=None):
        lsv_list = []
        splicegraph = []
        file_list = []

        junc_idx = [0] * majiq_config.num_experiments
        lsv_idx = [0] * majiq_config.num_experiments

        for exp_idx, exp in enumerate(majiq_config.exp_list):
            # Majiq file
            f = h5py.File(get_builder_majiq_filename(majiq_config.outDir, majiq_config.exp_list[exp_idx]),
                          'w', compression='gzip', compression_opts=9)
            file_list.append(f)

            effective_readlen = (majiq_config.readLen - 16) + 1
            f.create_dataset(LSV_JUNCTIONS_DATASET_NAME,
                             (majiq_config.nrandom_junctions, effective_readlen),
                             maxshape=(None, effective_readlen))
            f.create_dataset(CONST_JUNCTIONS_DATASET_NAME,
                             (majiq_config.nrandom_junctions, effective_readlen),
                             maxshape=(None, effective_readlen))
            lsv_list.append(f)

            # Splicegraph
            f_splicegraph = h5py.File(get_builder_splicegraph_filename(majiq_config.outDir,
                                                                       majiq_config.exp_list[exp_idx]),
                                      'w', compression='gzip', compression_opts=9)
            splicegraph.append(f_splicegraph)

        majiq_utils.monitor('AFTER CHILD CREATION AND FILES PREP')
        nthr_count = 0

        while True:
            try:
                val = result_queue.get(block=True, timeout=10)
                print "QUEUE SIZE", result_queue.qsize()
                if val[0] == 0:
                    for jdx, exp_idx in enumerate(majiq_config.tissue_repl[val[2]]):
                        lsv_idx[exp_idx] = val[1].to_hdf5(hdf5grp=lsv_list[exp_idx],
                                                          lsv_idx=lsv_idx[exp_idx],
                                                          exp_idx=jdx,
                                                          gc_func=vfunc_gc[exp_idx])

                elif val[0] == 1:
                    for jdx, exp_idx in enumerate(majiq_config.tissue_repl[val[3]]):
                        if junc_idx[exp_idx] >= majiq_config.nrandom_junctions:
                            continue
                        junc_group = lsv_list[exp_idx][CONST_JUNCTIONS_DATASET_NAME]
                        junc_group[junc_idx[exp_idx], :] = val[1][jdx, :].toarray()
                        junc_idx[exp_idx] += 1

                elif val[0] == 2:
                    val[1].to_hdf5(splicegraph[val[2]])

                elif val[0] == QUEUE_MESSAGE_END_WORKER:
                    lock_array[val[1]].release()
                    nthr_count += 1

            except Queue.Empty:
                if nthr_count < majiq_config.num_final_chunks:
                    continue
                break

        for exp_idx, exp in enumerate(majiq_config.exp_list):
            lsv_list[exp_idx].close()

        majiq_utils.monitor('MASTER END')

    def builder(self, sam_list):

        logger = majiq_utils.get_logger("%s/majiq.log" % majiq_config.outDir, silent=False,
                                        debug=self.debug)
        logger.info("")
        logger.info("Command: %s" % self)

        manager = mp.Manager()
        list_of_genes = manager.list()
        gc_pairs = manager.dict()

        p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
                       args=(majiq_io.read_gff, [self.transcripts, list_of_genes, gc_pairs, sam_list],
                             '%s/tmp' % majiq_config.outDir, 'db', 0, False))

        logger.info("... waiting gff3 parsing")
        p.start()
        p.join()

        if majiq_config.gcnorm:


            pool = mp.Pool(processes=self.nthreads, maxtasksperchild=1)
            lchnksize = max(len(sam_list)/self.nchunks, 1)
            lchnksize = lchnksize if len(sam_list) % self.nchunks == 0 else lchnksize + 1
            values = list(zip(range(len(sam_list)), sam_list))
            for vals in majiq_utils.chunks(values, lchnksize):
                pool.apply_async(majiq_io.gc_content_per_file, [vals, gc_pairs, majiq_config.outDir])
            pool.close()
            pool.join()
            vfunc_gc = majiq_norm.gc_normalization(gc_pairs, logger)

        else:
            vfunc_gc = [None] * majiq_config.num_experiments

        majiq_utils.monitor('AFTER READ GFF')

        lock_array = [mp.Lock() for xx in range(self.nthreads)]
        q = mp.Queue()

        pool = mp.Pool(processes=self.nthreads, initializer=builder_init,
                       initargs=[q, lock_array, sam_list, self.pcr_filename, self.gff_output,
                                 self.only_rna, self.non_denovo, get_build_temp_db_filename(majiq_config.outDir),
                                 self.silent, self.debug],
                       maxtasksperchild=1)

        lchnksize = max(len(list_of_genes)/self.nchunks, 1) + 1
        [xx.acquire() for xx in lock_array]

        pool.map_async(majiq_builder, majiq_utils.chunks(list_of_genes, lchnksize, extra=range(self.nchunks)))
        pool.close()
        self.queue_manager(lock_array, q, vfunc_gc=vfunc_gc, logger=logger)
        pool.join()
        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
