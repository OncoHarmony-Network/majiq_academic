import multiprocessing as mp
import os
import sys
import traceback
import h5py
import numpy as np
import majiq.src.io as majiq_io
import majiq.src.io_base as io_base
import majiq.grimoire.gene
import majiq.src.utils as majiq_utils
import majiq.src.multiproc as majiq_multi
from majiq.grimoire.exon import detect_exons
from majiq.src.normalize import gc_factor_calculation, gc_normalization
from majiq.src.analize import lsv_detection
from majiq.src.constants import *
from majiq.src.basic_pipeline import BasicPipeline, pipeline_run
from majiq.src.polyfitnb import fit_nb
from majiq.src.voila_wrapper import gene_to_splicegraph, init_splicegraph
from majiq.src.config import Config
import datetime


def build(args):
    pipeline_run(Builder(args))


def builder_init(lock_array, sam_list, pcr_filename, gff_output, only_rna,
                 non_denovo, silent, debug):

    builder_init.files_locks = lock_array
    builder_init.sam_list = sam_list
    builder_init.pcr_filename = pcr_filename
    builder_init.gff_output = gff_output
    builder_init.only_rna = only_rna
    builder_init.non_denovo = non_denovo
    builder_init.silent = silent
    builder_init.debug = debug


def merging_files(args_vals):

    try:
        # gne_id, chnk, loop_idx, total = args_vals
        majiq_config = Config()
        list_of_genes, chnk = args_vals
        logger = majiq_utils.get_logger("%s/%s.majiq.log" % (majiq_config.outDir, chnk),
                                        silent=majiq_config.silent, debug=majiq_config.debug)
        total = len(list_of_genes)
        for gne_idx, gne_id in enumerate(list_of_genes):
            loop_idx = gne_idx
            if loop_idx % 50 == 0:
                logger.info("[%s] Progress %s/%s" % (chnk, loop_idx, total))
            loop_id = '%s - %s' % (chnk, gne_id)
            logger.debug("[%s] Retrieving gene" % loop_id)
            junction_list = {}

            with h5py.File(get_build_temp_db_filename(majiq_config.outDir), 'r') as db_f:
                gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f, junction_list=junction_list, all_exp=True)

            splice_list = set()
            dict_of_junctions = {}

            rna_files = []
            vfunc_gc = []

            jset = set()
            for exp_idx in range(len(majiq_config.sam_list)):
                fname = get_builder_temp_majiq_filename(majiq_config.outDir, majiq_config.sam_list[exp_idx])
                with h5py.File(fname, 'r') as rfa:
                    try:
                        jset = jset.union(set(rfa["%s/junctions" % gne_id].keys()))
                        rna_files.append(fname)
                    except KeyError:
                        pass
                    if majiq_config.gcnorm:
                        vfunc_gc.append(gc_normalization(rfa.attrs['gc_values']))
                    else:
                        vfunc_gc = None

            njunc = len(jset)
            gene_obj.junc_matrix = np.zeros(shape=(njunc, len(majiq_config.sam_list), (majiq_config.readLen - 16 + 1)),
                                            dtype=np.uint32)
            if majiq_config.gcnorm:
                gene_obj.gc_content = np.zeros(shape=(njunc, len(majiq_config.sam_list), (majiq_config.readLen - 16 + 1)),
                                               dtype=np.float)

            for exp_idx, fname in enumerate(rna_files):
                with h5py.File(fname, 'r') as rnaf:
                    for jj_grp_id in rnaf["%s/junctions" % gne_id]:
                        jj_grp = rnaf["%s/junctions/%s" % (gne_id, jj_grp_id)]
                        junc = majiq.grimoire.gene.extract_junctions_hdf5(gene_obj, jj_grp, junction_list,
                                                                          annotated=jj_grp.attrs['annotated'],
                                                                          all_exp=True)
                        hdfidx = jj_grp.attrs['coverage_index']
                        gene_obj.junc_matrix[junc.get_index(), exp_idx, :] = rnaf[JUNCTIONS_DATASET_NAME][hdfidx, :]
                        if majiq_config.gcnorm:
                            gene_obj.gc_content[junc.get_index(), exp_idx, :] = rnaf[JUNCTIONS_GC_CONTENT][hdfidx, :]
                        if junc.intronic:
                            coord = junc.get_coordinates()
                            dict_of_junctions[coord[0]] = junc
                            dict_of_junctions[coord[1]] = junc
                        else:
                            splice_list.add((junc.start, '5prime', junc))
                            splice_list.add((junc.end, '3prime', junc))

            del junction_list
            majiq.grimoire.exon.detect_exons(gene_obj, list(splice_list), retrieve=True)
            del splice_list
            majiq.grimoire.gene.find_intron_retention(gene_obj, dict_of_junctions, majiq_config.non_denovo,
                                                      logging=logger)
            del dict_of_junctions
            if majiq_config.simplify:
                logger.debug('[%s] Simplifying gene' % loop_id)
                gene_obj.simplify()
            gene_obj.prepare_exons()

            logger.debug("[%s] Detecting LSV" % loop_id)
            nlsv = lsv_detection(gene_obj, gc_vfunc=vfunc_gc, lsv_list=majiq_config.sam_list,
                                 locks=builder_init.files_locks, rna_files=rna_files, logging=None)

            if nlsv:
                gene_to_splicegraph(gene_obj, builder_init.files_locks[-1])


            del majiq_config.gene_tlb[gne_id]
            del gene_obj


    except Exception:
        # majiq_utils.monitor('CHILD %s:: EXCEPT' % chnk)
        traceback.print_exc()
        sys.stdout.flush()
        raise

    finally:
        import logging
        logging.shutdown()




def parsing_files(args_vals):
    try:
        sam_file, chnk, loop_idx, total = args_vals
        majiq_config = Config()
        tlogger = majiq_utils.get_logger("%s/%s.majiq.log" % (majiq_config.outDir, chnk),
                                         silent=majiq_config.silent, debug=majiq_config.debug)

        tlogger.info("[%s] Starting new thread" % sam_file)
        majiq_utils.monitor('CHILD %s:: CREATION' % chnk)

        counter = [0] * 6
        loop_id = sam_file
        out_f = h5py.File(get_builder_temp_majiq_filename(majiq_config.outDir, sam_file), 'w')
        samfl = io_base.open_rnaseq("%s/%s.bam" % (majiq_config.sam_dir, sam_file))
        gc_pairs = {'GC': [], 'COV': []}
        jnc_idx = 0

        with h5py.File(get_build_temp_db_filename(majiq_config.outDir), 'r') as db_f:
            list_of_genes = db_f.keys()
            junc_mtrx = []
            gc_matrx = [] if majiq_config.gcnorm else None

            for gne_idx, gne_id in enumerate(list_of_genes):
                if gne_idx % 50 == 0:
                    tlogger.info("[%s] Progress %s/%s" % (loop_id, gne_idx, len(list_of_genes)))

                tlogger.debug("[%s] Retrieving gene" % gne_id)
                gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f)

                tlogger.debug("[%s] Reading BAM files" % gne_id)
                io_base.read_sam_or_bam(gene_obj, samfl, counter, h5py_file=db_f,
                                        nondenovo=majiq_config.non_denovo, info_msg=loop_id, logging=tlogger)

                if gene_obj.get_read_count() == 0:
                    continue

                out_f.create_group('%s/junctions' % gne_id)
                if majiq_config.gcnorm:
                    for ex in gene_obj.get_exon_list():
                        if ex.get_gc_content() > 0 and ex.get_coverage() > 0:
                            gc_pairs['GC'].append(ex.get_gc_content())
                            gc_pairs['COV'].append(ex.get_coverage())

                tlogger.debug("[%s] Detecting intron retention events" % gne_id)
                io_base.rnaseq_intron_retention(gene_obj, samfl, chnk,
                                                permissive=majiq_config.permissive_ir,
                                                nondenovo=majiq_config.non_denovo, logging=tlogger)

                for jnc in gene_obj.get_all_junctions():
                    jnc.to_rna_hdf5(out_f['%s/junctions' % gne_id], junc_mtrx,
                                    data_index=jnc_idx, gc_dataset=gc_matrx)

                    jnc_idx += 1

                del gene_obj
                del majiq_config.gene_tlb[gne_id]

        io_base.close_rnaseq(samfl)

        junc_mtrx = np.array(junc_mtrx)
        indx = np.arange(junc_mtrx.shape[0])[junc_mtrx.sum(axis=1) >= majiq_config.minreads]
        tlogger.debug("[%s] Fitting NB function with constitutive events..." % gne_id)
        out_f.attrs['one_over_r'] = fit_nb(junc_mtrx[indx, :], "%s/nbfit" % majiq_config.outDir,
                                           None, logger=tlogger)
        out_f.create_dataset(JUNCTIONS_DATASET_NAME, data=junc_mtrx, compression='gzip', compression_opts=9)
        if majiq_config.gcnorm:
            gc_matrx = np.array(gc_matrx)
            out_f.create_dataset(JUNCTIONS_GC_CONTENT, data=gc_matrx, compression='gzip', compression_opts=9,
                                 dtype=np.float)
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
        tlogger.info("[%s] End" % sam_file)
        tlogger.shutdown()


class Builder(BasicPipeline):

    def run(self):
        if self.simplify is not None and len(self.simplify) not in (0, 2):
            raise RuntimeError('Simplify requires 2 values type of junctions afected and E(PSI) threshold.')
        if not os.path.exists(self.conf):
            raise RuntimeError("Config file %s does not exist" % self.conf)
        majiq_config = Config(self.conf, self)
        self.builder(majiq_config)

    def builder(self, majiq_config):

        logger = majiq_utils.get_logger("%s/majiq.log" % majiq_config.outDir, silent=False,
                                        debug=self.debug)
        logger.info("Majiq Build v%s" % VERSION)
        logger.info("Command: %s" % " ".join(sys.argv))

        if self.prebam:

            p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
                           args=(majiq_io.read_gff, [self.transcripts, majiq_config.sam_list],
                                 '%s/tmp' % majiq_config.outDir, 'db', 0, False))

            logger.info("... waiting gff3 parsing")
            p.start()
            p.join()

            pool = mp.Pool(processes=self.nthreads, initializer=builder_init,
                           initargs=[None, majiq_config.sam_list, self.pcr_filename, self.gff_output, self.only_rna,
                                     self.non_denovo, self.silent, self.debug],
                           maxtasksperchild=1)

            lchnksize = max(int(len(majiq_config.sam_list)/self.nthreads), 1) + 1
            pool.map_async(parsing_files, majiq_utils.chunks(majiq_config.sam_list, lchnksize, extra=range(self.nthreads)))
            pool.close()
            pool.join()

        lock_array = [mp.Lock() for xx in majiq_config.sam_list]
        lock_array.append(mp.Lock())
        pool = mp.Pool(processes=self.nthreads, initializer=builder_init,
                       initargs=[lock_array, majiq_config.sam_list, self.pcr_filename, self.gff_output, self.only_rna,
                                 self.non_denovo, self.silent, self.debug],
                       maxtasksperchild=1)

        list_of_genes = None
        with h5py.File(get_build_temp_db_filename(majiq_config.outDir), 'r') as db_f:
                list_of_genes = list(db_f.keys())
        lchnksize = max(len(list_of_genes)/self.nthreads, 1) + 1
        init_splicegraph(get_builder_splicegraph_filename(majiq_config.outDir))
        for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            with h5py.File(get_builder_majiq_filename(majiq_config.outDir, sam_file), 'w') as f:
                effective_readlen = (majiq_config.readLen - 16) + 1
                f.create_dataset(JUNCTIONS_DATASET_NAME, (NRANDOM_JUNCTIONS, effective_readlen),
                                 maxshape=(None, effective_readlen))

                # fill meta info
                f.attrs['sample_id'] = sam_file
                path = get_builder_temp_majiq_filename(majiq_config.outDir, sam_file)
                f.attrs['fitfunc'] = majiq_utils.get_fitfunc_from_rnafile(path)
                f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
                f.attrs['data_index'] = 0
                f.attrs['num_lsvs'] = 0
                f.attrs['genome'] = majiq_config.genome
                f.attrs['VERSION'] = VERSION

        pool.map_async(merging_files, majiq_utils.chunks2(list_of_genes, lchnksize, range(self.nthreads)))
        pool.close()
        pool.join()

        for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            with h5py.File(get_builder_majiq_filename(majiq_config.outDir, sam_file), 'r+') as f:

                n_juncs = f.attrs['data_index']
                shp = f[JUNCTIONS_DATASET_NAME].shape
                f[JUNCTIONS_DATASET_NAME].resize((n_juncs, shp[1]))
                logger.info('%s LSVs found in %s' % (f.attrs['num_lsvs'], sam_file))

        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
