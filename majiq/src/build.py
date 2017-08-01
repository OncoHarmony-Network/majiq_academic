import multiprocessing as mp
import os
import sys
import traceback
import h5py

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
from majiq.src.multiproc import QueueMessage, process_conf, queue_manager, process_wrapper


def build(args):
    pipeline_run(Builder(args))


def merging_files(list_of_genes, chnk, majiq_config, process_conf, logger):
    try:
        total = len(list_of_genes)
        dummy = 0
        for gne_idx, gne_id in enumerate(list_of_genes):
            loop_idx = gne_idx
            if loop_idx % 50 == 0:
                logger.info("[%s] Progress %s/%s" % (chnk, loop_idx, total))
            loop_id = '%s - %s' % (chnk, gne_id)
            logger.debug("[%s] Retrieving gene" % loop_id)
            junction_list = {}

            with h5py.File(get_build_temp_db_filename(majiq_config.outDir), 'r') as db_f:
                gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f, junction_list=junction_list, all_exp=True)

            splice_list = {xx.get_coordinates(): xx for xx in gene_obj.get_all_junctions(filter=False)}
            dict_of_junctions = {}

            rna_files = []
            vfunc_gc = []

            jset = set()
            fitfunc_r = []
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
                    fitfunc_r.append(rfa.attrs['one_over_r'])

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
                            splice_list[(junc.start, junc.end)] = junc

            del junction_list
            majiq.grimoire.exon.detect_exons(gene_obj, splice_list, retrieve=True)
            del splice_list
            majiq.grimoire.gene.find_intron_retention(gene_obj, dict_of_junctions, majiq_config.non_denovo,
                                                      logging=logger)
            del dict_of_junctions
            if majiq_config.simplify:
                logger.debug('[%s] Simplifying gene' % loop_id)
                gene_obj.simplify()
            gene_obj.prepare_exons()

            logger.debug("[%s] Detecting LSV" % loop_id)

            # nlsv = lsv_detection(gene_obj, gc_vfunc=vfunc_gc, lsv_list=majiq_config.sam_list,
            #                      locks = builder_init.files_locks, rna_files = rna_files, logging = None)

            nlsv = lsv_detection(gene_obj, gc_vfunc=vfunc_gc, fitfunc_r=fitfunc_r, queue=process_conf.queue)

            dummy += nlsv
            if nlsv:
                gene_to_splicegraph(gene_obj, process_conf.lock[-1])
            # print(gene_obj.id, dummy)
            del majiq_config.gene_tlb[gne_id]
            del gene_obj

    except Exception:
        # majiq_utils.monitor('CHILD %s:: EXCEPT' % chnk)
        traceback.print_exc()
        sys.stdout.flush()
        raise

    finally:
        qm = QueueMessage(QUEUE_MESSAGE_END_WORKER, None, chnk)
        process_conf.queue.put(qm, block=True)
        process_conf.lock[chnk].acquire()
        process_conf.lock[chnk].release()
        import logging
        logging.shutdown()


# def parsing_files(args_vals):
    # try:
    #     sam_file_list, chnk = args_vals
    #     majiq_config = Config()
    #     tlogger = majiq_utils.get_logger("%s/%s.majiq.log" % (majiq_config.outDir, chnk),
    #                                      silent=majiq_config.silent, debug=majiq_config.debug)
    #
    #     majiq_utils.monitor('CHILD %s:: CREATION' % chnk)

def parsing_files(sam_file_list, chnk, majiq_config, process_conf, logger):

        db_f = h5py.File(get_build_temp_db_filename(majiq_config.outDir), 'r')
        list_of_genes = list(db_f.keys())
        for gne_idx, gne_id in enumerate(list_of_genes):
            logger.debug("[%s] Retrieving gene" % gne_id)
            majiq.grimoire.gene.retrieve_gene(gne_id, db_f)
        db_f.close()

        majiq_utils.monitor('CHILD %s:: LOADED GENES' % chnk)

        for sam_file in sam_file_list:
            logger.info("[%s] Starting new file" % sam_file)
            loop_id = sam_file
            out_f = h5py.File(get_builder_temp_majiq_filename(majiq_config.outDir, sam_file), 'w')
            samfl = io_base.open_rnaseq("%s/%s.bam" % (majiq_config.sam_dir, sam_file))
            gc_pairs = {'GC': [], 'COV': []}
            jnc_idx = 0

            junc_mtrx = []
            gc_matrx = [] if majiq_config.gcnorm else None

            for gne_idx, gne_id in enumerate(list_of_genes):

                if gne_idx % 50 == 0:
                    logger.info("[%s] Progress %s/%s" % (loop_id, gne_idx, len(list_of_genes)))

                gene_obj = majiq_config.gene_tlb[gne_id]
                gene_obj.reset_to_db()
                # tlogger.debug("[%s] Retrieving gene" % gne_id)
                # gene_obj = majiq.grimoire.gene.retrieve_gene(gne_id, db_f)

                logger.debug("[%s] Reading BAM files" % gne_id)
                io_base.read_sam_or_bam(gene_obj, samfl, h5py_file=db_f, info_msg=loop_id, logging=logger)

                if gene_obj.get_read_count() == 0:
                    continue

                out_f.create_group('%s/junctions' % gne_id)
                if majiq_config.gcnorm:
                    for ex in gene_obj.get_exon_list():
                        if ex.get_gc_content() > 0 and ex.get_coverage() > 0:
                            gc_pairs['GC'].append(ex.get_gc_content())
                            gc_pairs['COV'].append(ex.get_coverage())

                logger.debug("[%s] Detecting intron retention events" % gne_id)
                io_base.rnaseq_intron_retention(gene_obj, samfl, chnk, logging=logger)

                for jnc in gene_obj.get_all_junctions():
                    jnc.to_rna_hdf5(out_f['%s/junctions' % gne_id], junc_mtrx,
                                    data_index=jnc_idx, gc_dataset=gc_matrx)
                    jnc_idx += 1

                # del gene_obj
                # del majiq_config.gene_tlb[gne_id]

            io_base.close_rnaseq(samfl)

            junc_mtrx = np.array(junc_mtrx)
            indx = np.arange(junc_mtrx.shape[0])[junc_mtrx.sum(axis=1) >= majiq_config.minreads]
            logger.debug("[%s] Fitting NB function with constitutive events..." % sam_file)
            out_f.attrs['one_over_r'] = fit_nb(junc_mtrx[indx, :], "%s/nbfit" % majiq_config.outDir,
                                               None, logger=logger)
            out_f.create_dataset(JUNCTIONS_DATASET_NAME, data=junc_mtrx, compression='gzip', compression_opts=9)
            if majiq_config.gcnorm:
                gc_matrx = np.array(gc_matrx)
                out_f.create_dataset(JUNCTIONS_GC_CONTENT, data=gc_matrx, compression='gzip', compression_opts=9,
                                     dtype=np.float)
                factor, meanbins = gc_factor_calculation(gc_pairs, nbins=10)
                out_f.attrs['gc_values'] = (factor, meanbins)
            out_f.close()
        majiq_utils.monitor('CHILD %s:: ENDLOOP' % chnk)

    # except Exception:
    #     majiq_utils.monitor('CHILD %s:: EXCEPT' % chnk)
    #     traceback.print_exc()
    #     sys.stdout.flush()
    #     raise
    #
    # finally:
    #     tlogger.info("[%s] End" % sam_file)
    #     import logging
    #     logging.shutdown()


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

            if self.nthreads > 1:
                pool = mp.Pool(processes=self.nthreads, initializer=process_conf, initargs=[self, None, None, None])
                lchnksize = max(int(len(majiq_config.sam_list)/self.nthreads), 1) + 1
                pool.map_async(parsing_files, majiq_utils.chunks2(majiq_config.sam_list, lchnksize,
                                                                  extra=range(self.nthreads)))
                pool.close()
                pool.join()
            else:
                parsing_files(majiq_config.sam_list, 0, majiq_config=majiq_config, process_conf=self, logger=logger)

        with h5py.File(get_build_temp_db_filename(majiq_config.outDir), 'r') as db_f:
                list_of_genes = list(db_f.keys())

        lchnksize = max(len(list_of_genes)/self.nthreads, 1) + 1
        init_splicegraph(get_builder_splicegraph_filename(majiq_config.outDir))

        out_h5p_list = []
        for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            #TODO: CHANGE NAME BOOTS FOR MAJIQ
            f = h5py.File('%s/%s.boots.hdf5' % (majiq_config.outDir, sam_file), 'w')
            f.create_dataset('junctions', (5000,  majiq_config.m), maxshape=(None, majiq_config.m))
            f.create_dataset('junc_cov', (5000, 2), maxshape=(None, 2))

            # fill meta info
            f.attrs['sample_id'] = sam_file
            f.attrs['date'] = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
            f.attrs['VERSION'] = VERSION
            f.attrs['lsv_idx'] = 0
            f.attrs['num_lsvs'] = 0
            f.attrs['genome'] = majiq_config.genome
            out_h5p_list.append(f)

        self.queue = mp.Queue()
        self.lock = [mp.Lock() for xx in range(self.nthreads)]
        self.lock.append(mp.Lock())

        pool = mp.Pool(processes=self.nthreads, initializer=process_conf, initargs=[merging_files, self],
                       maxtasksperchild=1)
        [xx.acquire() for xx in self.lock[:-1]]

        pool.map_async(process_wrapper, majiq_utils.chunks2(list_of_genes, lchnksize, range(self.nthreads)))
        pool.close()

        queue_manager(None, out_h5p_list, self.lock, self.queue, num_chunks=self.nthreads, logger=logger)
        pool.join()
        for exp_idx, sam_file in enumerate(majiq_config.sam_list):
            f = out_h5p_list[exp_idx]
            n_juncs = f.attrs['lsv_idx']
            shp = f[JUNCTIONS_DATASET_NAME].shape
            f[JUNCTIONS_DATASET_NAME].resize((n_juncs, shp[1]))
            logger.info('%s LSVs found in %s' % (f.attrs['num_lsvs'], sam_file))
            f.close()

        logger.info("MAJIQ Builder is ended succesfully!")
        logger.info("Alakazam! Done.")
