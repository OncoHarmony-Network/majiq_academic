import Queue
import multiprocessing as mp
import os
import sys
import traceback
import h5py
import numpy as np

import majiq.src.io as majiq_io
import majiq.grimoire.gene
import majiq.src.config as majiq_config
import majiq.src.utils as majiq_utils
import majiq.src.multiproc as majiq_multi
from majiq.grimoire.exon import detect_exons
from majiq.src.normalize import gc_factor_calculation, gc_normalization
from majiq.src.analize import lsv_detection
from majiq.src.constants import *
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run


def build(args):
    pipeline_run(Builder(args))


def builder_init(sam_list, pcr_filename, gff_output, only_rna,
                 non_denovo, dbfile, list_of_genes, silent, debug):

    builder_init.sam_list = sam_list
    builder_init.pcr_filename = pcr_filename
    builder_init.gff_output = gff_output
    builder_init.only_rna = only_rna
    builder_init.non_denovo = non_denovo
    builder_init.dbfile = dbfile
    builder_init.silent = silent
    builder_init.debug = debug
    builder_init.list_of_genes = list_of_genes


def parsing_files(args_vals):

    filesnames, chnk = args_vals
    tlogger = majiq_utils.get_logger("%s/%s.majiq.log" % (majiq_config.outDir, chnk),
                                     silent=builder_init.silent, debug=builder_init.debug)

    tlogger.debug("[%s] Starting new thread" % chnk)
    majiq_utils.monitor('CHILD %s:: CREATION' % chnk)

    db_f = h5py.File(builder_init.dbfile)

    # samfile = [pysam.Samfile(xx, "rb") for xx in majiq_builder.sam_list]
    counter = [0] * 6
    for vals in filesnames:
        chnk, sam_file = vals
        loop_id = sam_file
        tlogger.info("Analizing %s" % loop_id)

        out_f = h5py.File(get_builder_temp_majiq_filename(majiq_config.outDir, sam_file),
                          'w', compression='gzip', compression_opts=9)
        effective_readlen = (majiq_config.readLen - MIN_BP_OVERLAP*2) + 1
        out_f.create_dataset(CONST_JUNCTIONS_DATASET_NAME,
                             (majiq_config.nrandom_junctions, effective_readlen),
                             maxshape=(None, effective_readlen))

        majiq_utils.monitor('CHILD %s:: STARTLOOP' % chnk)
        try:
            samfl = majiq_io.open_rnaseq("%s/%s.bam" % (majiq_config.sam_dir, sam_file))
            gc_pairs = {'GC': [], 'COV': []}
            jnc_idx = 0
            for gne_id in builder_init.list_of_genes:
                out_f.create_group('%s/junctions' % gne_id)
                tlogger.info("[%s] Retrieving gene" % loop_id)
                gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f)

                tlogger.info("[%s] Reading BAM files" % loop_id)
                majiq_io.read_sam_or_bam(gene_obj, samfl, counter, h5py_file=db_f,
                                         nondenovo=builder_init.non_denovo, info_msg=loop_id, logging=tlogger)

                if gene_obj.get_read_count() == 0:
                    continue

                if majiq_config.gcnorm:
                    for ex in gene_obj.get_exon_list():
                        gc_pairs['GC'].append(ex.get_gc_content())
                        gc_pairs['COV'].append(ex.get_coverage())

                tlogger.info("[%s] Detecting intron retention events" % loop_id)
                majiq_io.rnaseq_intron_retention(gene_obj, samfl, chnk,
                                                 permissive=majiq_config.permissive_ir,
                                                 nondenovo=builder_init.non_denovo, logging=tlogger)

                for jnc in gene_obj.get_all_junctions():
                    jnc.to_rna_hdf5(out_f['%s/junctions' % gne_id], out_f[CONST_JUNCTIONS_DATASET_NAME],
                                    data_index=jnc_idx)
                    jnc_idx += 1

                del gene_obj
                del majiq_config.gene_tlb[gne_id]

            factor, meanbins = gc_factor_calculation(gc_pairs, nbins=10)
            out_f.attrs['gc_values'] = (factor, meanbins)
            out_f.close()
            majiq_utils.monitor('CHILD %s:: ENDLOOP' % chnk)

        except Exception:
            majiq_utils.monitor('CHILD %s:: EXCEPT' % chnk)
            traceback.print_exc()
            sys.stdout.flush()
            raise

        finally:
            majiq_io.close_rnaseq(samfl)

    db_f.close()
    # majiq_utils.monitor('CHILD %s:: WAITING' % chnk)
    # tlogger.info("[%s] Waiting to be freed" % chnk)
    # qm = majiq_multi.QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
    # majiq_builder.queue.put(qm, block=True)
    # majiq_builder.lock_arr[chnk].acquire()
    # majiq_builder.lock_arr[chnk].release()
    tlogger.info("[%s] Child work done." % chnk)


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
                # print "QUEUE SIZE", result_queue.qsize()
                if val.get_type() == QUEUE_MESSAGE_BUILD_LSV:
                    for jdx, exp_idx in enumerate(majiq_config.tissue_repl[val.get_value()[1]]):
                        lsvobj = val.get_value()[0]
                        lsv_idx[exp_idx] = lsvobj.to_hdf5(hdf5grp=lsv_list[exp_idx],
                                                          lsv_idx=lsv_idx[exp_idx],
                                                          exp_idx=jdx,
                                                          gc_func=vfunc_gc[exp_idx])

                elif val.get_type() == QUEUE_MESSAGE_BUILD_CONST_JUNCTION:
                    juncs = val.get_value()[0]
                    for jdx, exp_idx in enumerate(majiq_config.tissue_repl[val.get_value()[1]]):
                        if junc_idx[exp_idx] >= majiq_config.nrandom_junctions:
                            continue
                        elif junc_idx[exp_idx] + juncs.shape[0] >= majiq_config.nrandom_junctions:
                            shp = lsv_list[exp_idx][CONST_JUNCTIONS_DATASET_NAME].shape
                            shp_new = shp[0] + juncs.shape[0]
                            lsv_list[exp_idx][CONST_JUNCTIONS_DATASET_NAME].resize((shp_new, shp[1]))

                        lsv_list[exp_idx][CONST_JUNCTIONS_DATASET_NAME][junc_idx[exp_idx]:junc_idx[exp_idx]+juncs.shape[0], :] = juncs[:, jdx, :]
                        junc_idx[exp_idx] += juncs.shape[0]

                elif val.get_type() == QUEUE_MESSAGE_SPLICEGRAPH:
                    val[1].to_hdf5(splicegraph[val[2]])

                elif val.get_type() == QUEUE_MESSAGE_END_WORKER:
                    lock_array[val.get_chunk()].release()
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

        p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
                       args=(majiq_io.read_gff, [self.transcripts, list_of_genes, sam_list],
                             '%s/tmp' % majiq_config.outDir, 'db', 0, False))

        logger.info("... waiting gff3 parsing")
        p.start()
        p.join()
        majiq_utils.monitor('AFTER READ GFF')

        # if majiq_config.gcnorm:
        #     pool = mp.Pool(processes=self.nthreads, maxtasksperchild=1)
        #     lchnksize = max(len(sam_list)/self.nchunks, 1)
        #     lchnksize = lchnksize if len(sam_list) % self.nchunks == 0 else lchnksize + 1
        #     values = list(zip(range(len(sam_list)), sam_list))
        #     output_gc_vals = manager.dict()
        #     for vals in majiq_utils.chunks(values, lchnksize):
        #         pool.apply_async(majiq_io.gc_content_per_file, [vals, output_gc_vals, majiq_config.outDir])
        #     pool.close()
        #     pool.join()
        #     vfunc_gc = majiq_norm.gc_normalization(output_gc_vals, logger)
        #
        # else:
        #     vfunc_gc = [None] * majiq_config.num_experiments

        if self.prebam:
            pool = mp.Pool(processes=self.nthreads, initializer=builder_init,
                           initargs=[sam_list, self.pcr_filename, self.gff_output, self.only_rna, self.non_denovo,
                                     get_build_temp_db_filename(majiq_config.outDir), list_of_genes,
                                     self.silent, self.debug],
                           maxtasksperchild=1)
            lchnksize = max(len(sam_list)/self.nchunks, 1)
            lchnksize = lchnksize if len(sam_list) % self.nchunks == 0 else lchnksize + 1
            values = list(zip(range(len(sam_list)), sam_list))
            pool.map_async(parsing_files, majiq_utils.chunks(values, lchnksize, extra=range(self.nthreads)))
            pool.close()
            pool.join()


        # VALUES
        db_f = h5py.File(get_build_temp_db_filename(majiq_config.outDir))

        out_files = []
        rna_files = []
        vfunc_gc = []
        for exp_idx, sam_file in enumerate(sam_list):
            rnaf = h5py.File(get_builder_temp_majiq_filename(majiq_config.outDir, sam_file))
            rna_files.append(rnaf)
            if majiq_config.gcnorm:
                vfunc_gc.append(gc_normalization(rnaf.attrs['gc_values']))
            else:
                vfunc_gc.append(None)

            f = h5py.File(get_builder_majiq_filename(majiq_config.outDir, sam_file),
                          'w', compression='gzip', compression_opts=9)
            out_files.append(f)

            effective_readlen = (majiq_config.readLen - 16) + 1
            f.create_dataset(LSV_JUNCTIONS_DATASET_NAME,
                             (majiq_config.nrandom_junctions, effective_readlen),
                             maxshape=(None, effective_readlen))
            jj_list = rnaf[CONST_JUNCTIONS_DATASET_NAME][()]
            indx = np.arange(jj_list.shape[0])[jj_list.sum(axis=1) >= 2]

            f.create_dataset(CONST_JUNCTIONS_DATASET_NAME, data=jj_list[indx, :])

        out_files_idx = [0] * majiq_config.num_experiments
        for gne_id in list_of_genes:
            loop_id = 'MASTER - %s' % (gne_id)
            logger.info("[%s] Retrieving gene" % loop_id)
            gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f, all_exp=True)

            junction_list = {}

            splice_list = set()

            for exp_idx, filename in enumerate(sam_list):
                for jj_grp_id in rna_files[exp_idx]["%s/junctions" % gne_id]:
                    jj_grp = rna_files[exp_idx]["%s/junctions/%s" % (gne_id, jj_grp_id)]
                    annot = jj_grp.attrs['annotated']
                    junc = majiq.grimoire.gene.extract_junctions_hdf5(gene_obj, jj_grp, junction_list, annotated=annot,
                                                                      all_exp=True)
                    junc.set_coverage(exp_idx, rna_files[exp_idx][CONST_JUNCTIONS_DATASET_NAME][jj_grp.attrs['coverage_index'], :])
                    splice_list.add((junc.start, '5prime', junc))
                    splice_list.add((junc.end, '3prime', junc))

            detect_exons(gene_obj, list(splice_list), None)

            logger.info("[%s] Detecting LSV" % loop_id)
            lsv_detection(gene_obj, vfunc_gc, out_files, out_files_idx, only_real_data=self.only_rna, logging=logger)

        ''' Closing HDF5 files'''
        db_f.close()
        [xx.close() for xx in rna_files]
        [xx.close() for xx in out_files]

        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
