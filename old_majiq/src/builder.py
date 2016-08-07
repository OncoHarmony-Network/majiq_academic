#!/usr/bin/python
import Queue
import multiprocessing as mp
import os
import sys
import traceback
from multiprocessing import Pool

import h5py

import majiq.src.config as majiq_config
import majiq.src.utils as utils
import old_majiq.grimoire.gene as majiq_gene
import old_majiq.src.analize as analize
import old_majiq.src.io as majiq_io
import old_majiq.src.io_utils
import old_majiq.src.voila_wrapper
from old_majiq.src.normalize import prepare_gc_content, gc_factor_calculation


def __builder_init(out_queue, lock_arr):
    majiq_builder.queue = out_queue
    majiq_builder.lock_arr = lock_arr


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
    gene_list = old_majiq.src.io_utils.load_bin_file(annot_file)

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
    analize.lsv_detection(gene_list, chnk, only_real_data=only_rna, out_queue=majiq_builder.queue,
                          logging=logging)


    prepare_gc_content(gene_list, temp_dir)

    # if pcr_validation:
    #     utils.get_validated_pcr_lsv(lsv, temp_dir)
    # if gff_output:
    #     majiq_lsv.extract_gff(lsv, temp_dir)
    old_majiq.src.voila_wrapper.generate_visualization_output(gene_list, temp_dir, majiq_builder.queue)
    if not logging is None:
        logging.info("[%s] Preparing output" % chnk)
    lsv = None
    const = None
    utils.send_output(lsv, const, temp_dir, majiq_builder.queue, chnk, majiq_builder.lock_arr[chnk])


def gather_files(out_dir, prefix='', gff_out=None, pcr_out=None, nthreads=1, logger=None):

    #GATHER
    logger.info("Gather outputs")
    if prefix != '':
        prefix = '%s.' % prefix

    if majiq_config.gcnorm:
        gc_factor_calculation(10)

    if nthreads > 1:
        nthr = min(nthreads, 4)
        pool = mp.Pool(processes=nthr)#, maxtasksperchild=1)

    for name, ind_list in majiq_config.tissue_repl.items():
        for idx, exp_idx in enumerate(ind_list):
            if nthreads > 1:
                pool.apply_async(utils.merge_and_create_majiq_file, [exp_idx, prefix])
            else:
                utils.merge_and_create_majiq_file(exp_idx, prefix)

    if nthreads > 1:
        pool.close()
        pool.join()

    if not gff_out is None:
        logger.info("Gather lsv and generate gff")
        fp = open('%s/%s' % (out_dir, gff_out), 'w+')
        for chnk in range(majiq_config.num_final_chunks):
            temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
            yfile = '%s/temp_gff.pkl' % temp_dir
            if not os.path.exists(yfile):
                continue
            gff_list = old_majiq.src.io_utils.load_bin_file(yfile)
            for gff in gff_list:
                fp.write("%s\n" % gff)
        fp.close()

    if not pcr_out is None:
        logger.info("Gather pcr results")
        fp = open('%s/pcr_match.tab' % majiq_config.outDir, 'w+')
        for chnk in range(majiq_config.num_final_chunks):
            temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
            yfile = '%s/pcr.pkl' % temp_dir
            if not os.path.exists(yfile):
                continue
            pcr_l = old_majiq.src.io_utils.load_bin_file(yfile)
            for pcr in pcr_l:
                fp.write("%s\n" % pcr)
        fp.close()


def __parallel_lsv_quant(samfiles_list, chnk, pcr_validation=False, gff_output=None, only_rna=False,
                         nondenovo=False, silent=False, debug=0):
    try:
        print "START child,", mp.current_process()._identity
        tlogger = utils.get_logger("%s/%s.old_majiq.log" % (majiq_config.outDir, mp.current_process()._identity[0]),
                                   silent=silent, debug=debug)
        majiq_builder(samfiles_list, chnk, pcr_validation=pcr_validation,
                      gff_output=gff_output, only_rna=only_rna, nondenovo=nondenovo, logging=tlogger)
        print "END child, ", mp.current_process().name
    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


def __parallel_gff3(transcripts, pcr_filename, nthreads, silent=False, debug=0):

    try:
        print "START child,", mp.current_process()._identity
        tlogger = utils.get_logger("%s/db.old_majiq.log" % majiq_config.outDir,
                                   silent=silent, debug=debug)
        majiq_io.read_gff(transcripts, pcr_filename, nthreads, logging=tlogger)
        print "END child, ", mp.current_process().name
    except Exception as e:
        traceback.print_exc()
        sys.stdout.flush()
        raise()


#########
# MAIN  #
#########


def main(params):

#    import yappi
#    yappi.start()

    if not os.path.exists(params.conf):
        raise RuntimeError("Config file %s does not exist" % params.conf)
    majiq_config.global_conf_ini(params.conf, params)

    logger = utils.get_logger("%s/old_majiq.log" % majiq_config.outDir, silent=params.silent, debug=params.debug)
    logger.info("")
    logger.info("Command: %s" % params)

    if not params.onlygather:
        p = mp.Process(target=__parallel_gff3, args=(params.transcripts, params.pcr_filename, params.nthreads))
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

        if len(sam_list) == 0:
            return
        if params.nthreads > 1:
            lock_array = [mp.Lock() for xx in range(majiq_config.num_final_chunks)]
            q = mp.Queue()
            pool = Pool(processes=params.nthreads, initializer=__builder_init,
                        initargs=[q, lock_array], maxtasksperchild=1)

        args_list = []
        for chnk in range(majiq_config.num_final_chunks):
            temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
            utils.create_if_not_exists(temp_dir)

            if params.nthreads == 1:
                majiq_builder(sam_list, chnk, pcr_validation=params.pcr_filename, gff_output=params.gff_output,
                              only_rna=params.only_rna, nondenovo=params.non_denovo, logging=logger)
            else:
                lock_array[chnk].acquire()
                pool.apply_async(__parallel_lsv_quant, [sam_list, chnk,
                                                      params.pcr_filename,
                                                      params.gff_output,
                                                      params.only_rna,
                                                      params.non_denovo,
                                                      params.silent,
                                                      params.debug
                                                      ])

        count = 0

        lsv_list = []
        junc_list = []
        splicegraph = []
        file_list = []
        junc_idx = [0] * majiq_config.num_experiments
        for exp_idx, exp in enumerate(majiq_config.exp_list):
            fname = "%s/%s.old_majiq.hdf5" % (majiq_config.outDir, majiq_config.exp_list[exp_idx])
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

        if params.nthreads > 1:
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

#    yappi.get_func_stats().print_all()

    majiq_config.print_numbers()
    logger.info("End of execution")
