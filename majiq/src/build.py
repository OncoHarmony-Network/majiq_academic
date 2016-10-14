import Queue
import multiprocessing as mp
import os
import sys
import traceback
import h5py
import numpy as np
from progressbar import ProgressBar

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
from majiq.src.polyfitnb import fit_nb
from majiq.src.voila_wrapper import gene_to_splicegraph, init_splicegraph
import datetime


def build(args):
    pipeline_run(Builder(args))


def builder_init(idx_count, lock_array, sam_list, pcr_filename, gff_output, only_rna,
                 non_denovo, dbfile, list_of_genes, silent, debug):

    builder_init.idx_count = idx_count
    builder_init.files_locks = lock_array
    builder_init.sam_list = sam_list
    builder_init.pcr_filename = pcr_filename
    builder_init.gff_output = gff_output
    builder_init.only_rna = only_rna
    builder_init.non_denovo = non_denovo
    builder_init.dbfile = dbfile
    builder_init.silent = silent
    builder_init.debug = debug
    builder_init.list_of_genes = list_of_genes


def merging_files(args_vals):

    list_of_genes, chnk = args_vals
    logger = majiq_utils.get_logger("%s/%s.majiq.log" % (majiq_config.outDir, chnk),
                                    silent=builder_init.silent, debug=builder_init.debug)

    try:
        rna_files = []
        vfunc_gc = []
        for exp_idx, sam_file in enumerate(builder_init.sam_list):
            rnaf = h5py.File(get_builder_temp_majiq_filename(majiq_config.outDir, sam_file))
            rna_files.append(rnaf)
            if majiq_config.gcnorm:
                vfunc_gc.append(gc_normalization(rnaf.attrs['gc_values']))
            else:
                vfunc_gc.append(None)

        db_f = h5py.File(builder_init.dbfile)
        for gne_idx, gne_id in enumerate(list_of_genes):
            if gne_idx % 50 == 0:
                logger.info("[%s] Progress %s/%s" % (chnk, gne_idx, len(list_of_genes)))
            loop_id = '%s - %s' % (chnk, gne_id)
            logger.debug("[%s] Retrieving gene" % loop_id)
            gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f, all_exp=True)

            junction_list = {}
            splice_list = set()

            for exp_idx, filename in enumerate(builder_init.sam_list):
                for jj_grp_id in rna_files[exp_idx]["%s/junctions" % gne_id]:
                    jj_grp = rna_files[exp_idx]["%s/junctions/%s" % (gne_id, jj_grp_id)]
                    annot = jj_grp.attrs['annotated']
                    junc = majiq.grimoire.gene.extract_junctions_hdf5(gene_obj, jj_grp, junction_list, annotated=annot,
                                                                      all_exp=True)
                    junc.set_coverage(exp_idx,
                                      rna_files[exp_idx][CONST_JUNCTIONS_DATASET_NAME][jj_grp.attrs['coverage_index'], :])
                    splice_list.add((junc.start, '5prime', junc))
                    splice_list.add((junc.end, '3prime', junc))

            detect_exons(gene_obj, list(splice_list), retrieve=True)
            del splice_list

            logger.debug("[%s] Detecting LSV" % loop_id)
            lsv_detection(gene_obj, gc_vfunc=vfunc_gc, lsv_list=builder_init.sam_list, lsv_idx=builder_init.idx_count,
                          locks=builder_init.files_locks, logging=None)

            del majiq_config.gene_tlb[gne_id]
            del gene_obj

    except Exception:
        majiq_utils.monitor('CHILD %s:: EXCEPT' % chnk)
        traceback.print_exc()
        sys.stdout.flush()
        raise

    finally:
        [xx.close() for xx in rna_files]
        db_f.close()

    logger.info("[%s] End" % chnk)


def parsing_files(args_vals):

    filesnames, chnk = args_vals
    tlogger = majiq_utils.get_logger("%s/%s.majiq.log" % (majiq_config.outDir, chnk),
                                     silent=builder_init.silent, debug=builder_init.debug)

    tlogger.debug("[%s] Starting new thread" % chnk)
    majiq_utils.monitor('CHILD %s:: CREATION' % chnk)
    db_f = h5py.File(builder_init.dbfile)
    counter = [0] * 6

    for vals in filesnames:
        chnk, sam_file = vals
        loop_id = sam_file
        out_f = h5py.File(get_builder_temp_majiq_filename(majiq_config.outDir, sam_file),
                          'w', compression='gzip', compression_opts=9)
        effective_readlen = (majiq_config.readLen - MIN_BP_OVERLAP*2) + 1
        out_f.create_dataset(CONST_JUNCTIONS_DATASET_NAME,
                             (majiq_config.nrandom_junctions, effective_readlen),
                             maxshape=(None, effective_readlen))

        sgraph = init_splicegraph(get_builder_splicegraph_filename(majiq_config.outDir, sam_file))


        try:
            samfl = majiq_io.open_rnaseq("%s/%s.bam" % (majiq_config.sam_dir, sam_file))
            gc_pairs = {'GC': [], 'COV': []}
            jnc_idx = 0

            for gne_idx, gne_id in enumerate(builder_init.list_of_genes):
                if gne_idx % 50 == 0:
                    tlogger.info("[%s] Progress %s/%s" % (loop_id, gne_idx, len(builder_init.list_of_genes)))
                out_f.create_group('%s/junctions' % gne_id)
                tlogger.debug("[%s] Retrieving gene" % gne_id)
                gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f)

                tlogger.debug("[%s] Reading BAM files" % gne_id)
                majiq_io.read_sam_or_bam(gene_obj, samfl, counter, h5py_file=db_f,
                                         nondenovo=builder_init.non_denovo, info_msg=loop_id, logging=tlogger)

                if gene_obj.get_read_count() == 0:
                    continue

                if majiq_config.gcnorm:
                    for ex in gene_obj.get_exon_list():
                        gc_pairs['GC'].append(ex.get_gc_content())
                        gc_pairs['COV'].append(ex.get_coverage())

                tlogger.debug("[%s] Detecting intron retention events" % gne_id)
                majiq_io.rnaseq_intron_retention(gene_obj, samfl, chnk,
                                                 permissive=majiq_config.permissive_ir,
                                                 nondenovo=builder_init.non_denovo, logging=tlogger)

                for jnc in gene_obj.get_all_junctions():
                    jnc.to_rna_hdf5(out_f['%s/junctions' % gne_id], out_f[CONST_JUNCTIONS_DATASET_NAME],
                                    data_index=jnc_idx)
                    jnc_idx += 1

                gene_to_splicegraph(gene_obj, sgraph)

                del gene_obj
                del majiq_config.gene_tlb[gne_id]

            if majiq_config.gcnorm:
                factor, meanbins = gc_factor_calculation(gc_pairs, nbins=10)
                out_f.attrs['gc_values'] = (factor, meanbins)

            jj_list = out_f[CONST_JUNCTIONS_DATASET_NAME][()]
            indx = np.arange(jj_list.shape[0])[jj_list.sum(axis=1) >= majiq_config.MINREADS]
            tlogger.debug("[%s] Fitting NB function with constitutive events..." % gne_id)
            out_f.attrs['one_over_r'] = fit_nb(jj_list[indx, :], "%s/nbfit" % majiq_config.outDir,
                                               None, logger=tlogger)

            out_f.close()
            sgraph.close()
            majiq_utils.monitor('CHILD %s:: ENDLOOP' % chnk)

        except Exception:
            majiq_utils.monitor('CHILD %s:: EXCEPT' % chnk)
            traceback.print_exc()
            sys.stdout.flush()
            raise

        finally:
            majiq_io.close_rnaseq(samfl)

    db_f.close()
    tlogger.info("[%s] End" % chnk)


class Builder(BasicPipeline):

    def run(self):
        if not os.path.exists(self.conf):
            raise RuntimeError("Config file %s does not exist" % self.conf)
        sam_list = majiq_config.global_conf_ini(self.conf, self)
        self.builder(sam_list)

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

        if self.prebam:
            pool = mp.Pool(processes=self.nthreads, initializer=builder_init,
                           initargs=[None, None, sam_list, self.pcr_filename, self.gff_output, self.only_rna,
                                     self.non_denovo, get_build_temp_db_filename(majiq_config.outDir), list_of_genes,
                                     self.silent, self.debug],
                           maxtasksperchild=1)
            lchnksize = max(len(sam_list)/self.nchunks, 1)
            lchnksize = lchnksize if len(sam_list) % self.nchunks == 0 else lchnksize + 1
            values = list(zip(range(len(sam_list)), sam_list))
            pool.map_async(parsing_files, majiq_utils.chunks(values, lchnksize, extra=range(self.nthreads)))
            pool.close()
            pool.join()


        # Detect LSVs

        lock_array = [mp.Lock() for xx in sam_list]
        idx_array = [0] * len(sam_list)

        pool = mp.Pool(processes=self.nthreads, initializer=builder_init,
                       initargs=[idx_array, lock_array, sam_list, self.pcr_filename, self.gff_output,
                                 self.only_rna, self.non_denovo, get_build_temp_db_filename(majiq_config.outDir),
                                 None, self.silent, self.debug],
                       maxtasksperchild=1)

        lchnksize = max(len(list_of_genes)/self.nthreads, 1) + 1
        db_f = h5py.File(get_build_temp_db_filename(majiq_config.outDir))

        for exp_idx, sam_file in enumerate(sam_list):
            f = h5py.File(get_builder_majiq_filename(majiq_config.outDir, sam_file),
                          'w', compression='gzip', compression_opts=9)
            effective_readlen = (majiq_config.readLen - 16) + 1
            f.create_dataset(LSV_JUNCTIONS_DATASET_NAME,
                             (2, effective_readlen),
                             maxshape=(None, effective_readlen))


            # fill meta info
            f.attrs['sample_id'] = sam_file
            path = get_builder_temp_majiq_filename(majiq_config.outDir, sam_file)
            f.attrs['fitfunc'] = majiq_utils.get_fitfunc_from_rnafile(path)
            f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
            f.attrs['VERSION'] = VERSION

            f.close()

        # majiq_multi.queue_manager(input_h5dfp=None, output_h5dfp=out_files, lock_array=lock_array, result_queue=q,
        #                           num_chunks=self.nthreads, logger=self.logger)
        pool.map_async(merging_files, majiq_utils.chunks(list_of_genes, lchnksize, extra=range(self.nthreads)))
        pool.close()
        pool.join()

        ''' Closing HDF5 files'''
        db_f.close()
        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
