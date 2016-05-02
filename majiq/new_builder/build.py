import os
import multiprocessing as mp

import majiq.new_builder.gene
from majiq.src.basic_pipeline import BasicPipeline, _pipeline_run
import majiq.src.config as majiq_config
import majiq.src.utils.utils as majiq_utils
import multiproc as majiq_multi
import majiq_io as majiq_io
from analize import lsv_detection
import pysam
from multiprocessing import Pool
import Queue
import sys
import traceback
import h5py
import types

def build(args):
    _pipeline_run(Builder(args))

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


def majiq_builder(list_of_genes):

    created = mp.Process()
    current = mp.current_process()
    print 'running:', current.name, current._identity
    print 'created:', created.name, created._identity

    tlogger = majiq_utils.get_logger("%s/%s.kk.majiq.log" % (majiq_config.outDir, mp.current_process()._identity[0]),
                                     silent=majiq_builder.silent, debug=majiq_builder.debug)

    tlogger.debug("[%s] Starting new chunk" % created._identity[1])
    chnk = created._identity[1]
    db_f = h5py.File(majiq_builder.dbfile)
    if isinstance(list_of_genes, types.StringTypes):
        list_of_genes = [list_of_genes]

    tlogger.info(list_of_genes)


    counter = [0] * 6
    for gne_id in list(list_of_genes):
        try:
            tlogger.info("[%s] Retrieving gene" % chnk)
            gene_obj = majiq.new_builder.gene.retrieve_gene(gne_id, db_f)

            tlogger.info("[%s] Reading BAM files" % chnk)
            samfile = [pysam.Samfile(xx, "rb") for xx in majiq_builder.sam_list]

            majiq_io.read_sam_or_bam(gene_obj, samfile, chnk, counter,
                                     nondenovo=majiq_builder.non_denovo, logging=tlogger)

            if gene_obj.get_read_count().sum() == 0:
                continue


            tlogger.info("[%s] Detecting intron retention events" % chnk)
            majiq_io.rnaseq_intron_retention(gene_obj, samfile, chnk,
                                             permissive=majiq_config.permissive_ir,
                                             nondenovo=majiq_builder.non_denovo, logging=tlogger)

            for ss in samfile:
                ss.close()

            tlogger.info("[%s] Detecting LSV" % chnk)
            lsv_detection(gene_obj, only_real_data=majiq_builder.only_rna, out_queue=majiq_builder.queue, logging=tlogger)

        except Exception as e:
            traceback.print_exc()
            sys.stdout.flush()
        finally:
            del majiq_config.gene_tlb[gne_id]



    db_f.close()
    tlogger.info("[%s] Waiting to be freed" % chnk)
    majiq_builder.queue.put([3, mp.current_process()._identity[0]], block=True)
    majiq_builder.lock_arr[chnk].acquire()
    majiq_builder.lock_arr[chnk].release()
    return


class Builder(BasicPipeline):
    def run(self):
        if not os.path.exists(self.conf):
            raise RuntimeError("Config file %s does not exist" % self.conf)
        sam_list = majiq_config.global_conf_ini(self.conf, self)
        self.builder(sam_list)


    def queue_manager(self, lock_array, result_queue, first_id=0):

        count = []
        lsv_list = []
        junc_list = []
        splicegraph = []
        file_list = []
        junc_idx = [0] * majiq_config.num_experiments
        lsv_idx = [0] * majiq_config.num_experiments

        for exp_idx, exp in enumerate(majiq_config.exp_list):
            fname = "%s/%s.majiq.hdf5" % (majiq_config.outDir, majiq_config.exp_list[exp_idx])
            f = h5py.File(fname, 'w', compression='gzip', compression_opts=9)
            file_list.append(f)

            fname_sg = "%s/%s.splicegraph.hdf5" % (majiq_config.outDir, majiq_config.exp_list[exp_idx])
            f_splicegraph = h5py.File(fname_sg, 'w', compression='gzip', compression_opts=9)

            #as_table = f.create_group('LSVs')
            lsv_list.append(f)
            effective_readlen = (majiq_config.readLen - 16) + 1
            as_table = f.create_dataset("/lsv_junctions",
                                        (majiq_config.nrandom_junctions, effective_readlen),
                                        maxshape=(None, effective_readlen))
            non_as_table = f.create_dataset("/const_junctions",
                                            (majiq_config.nrandom_junctions, effective_readlen),
                                            maxshape=(None, effective_readlen))

            junc_list.append(non_as_table)

            splicegraph.append(f_splicegraph)

        nthr_count = 0
        while True:
            try:
                val = result_queue.get(block=True, timeout=10)
                if val[0] == 0:
                    for exp_idx in majiq_config.tissue_repl[val[2]]:
                        lsv_idx[exp_idx] = val[1].to_hdf5(lsv_list[exp_idx], lsv_idx[exp_idx], exp_idx)

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
                    nthr_count += 1

            except Queue.Empty:
                if nthr_count < majiq_config.num_final_chunks:
                    continue
                break

    def builder(self, sam_list=[]):

        logger = majiq_utils.get_logger("%s/majiq.log" % majiq_config.outDir, silent=False,
                                        debug=self.debug)
        logger.info("")
        logger.info("Command: %s" % self)

        manager = mp.Manager()

        list_of_genes = manager.list()
        p = mp.Process(target=majiq_multi.parallel_lsv_child_calculation,
                       args=(majiq_io.read_gff, [self.transcripts, self.pcr_filename, self.nthreads, list_of_genes],
                             '%s/tmp' % majiq_config.outDir, 'db', 0, False))

        logger.info("... waiting gff3 parsing")
        p.start()
        p.join()

        #lock_array = [mp.Lock() for xx in range(majiq_config.num_final_chunks)]
        lock_array = {}
        q = mp.Queue()
        db_filename = "%s/tmp/db.hdf5" % majiq_config.outDir

        pool = mp.Pool(processes=self.nthreads, initializer=builder_init,
                       initargs=[q, lock_array, sam_list, self.pcr_filename,self.gff_output,
                                 self.only_rna, self.non_denovo, db_filename, self.silent, self.debug],
                       maxtasksperchild=1)

        lchnksize = max(len(list_of_genes)/self.nchunks, 1)

        for proc in pool._pool:
            pid = proc._identity[0]
            lock_array[pid] = mp.Lock()
            lock_array[pid].acquire()

        pool.map_async(majiq_builder, list_of_genes)#, chunksize=lchnksize)
        pool.close()
        self.queue_manager(lock_array, q)
        pool.join()





    #

    #
    #     for chnk in range(majiq_config.num_final_chunks):
    #         temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
    #         majiq_utils.create_if_not_exists(temp_dir)
    #         if self.nthreads == 1:
    #             majiq_builder(sam_list, chnk, pcr_validation=self.pcr_filename, gff_output=self.gff_output,
    #                           only_rna=self.only_rna, nondenovo=self.non_denovo, logging=logger)
    #         else:
    #             lock_array[chnk].acquire()
    #             pool.apply_async(majiq_multi.parallel_lsv_child_calculation, [majiq_builder, [sam_list, chnk,
    #                                                                                           self.pcr_filename,
    #                                                                                           self.gff_output,
    #                                                                                           self.only_rna,
    #                                                                                           self.non_denovo,
    #                                                                                           self.silent,
    #                                                                                           self.debug]
    #                                                                                          ])
    #
    #     count = 0
    #
    #     lsv_list = []
    #     junc_list = []
    #     splicegraph = []
    #     file_list = []
    #     junc_idx = [0] * majiq_config.num_experiments
    #     for exp_idx, exp in enumerate(majiq_config.exp_list):
    #         fname = "%s/%s.majiq.hdf5" % (majiq_config.outDir, majiq_config.exp_list[exp_idx])
    #         f = h5py.File(fname, 'w', compression='gzip', compression_opts=9)
    #         file_list.append(f)
    #
    #         fname_sg = "%s/%s.splicegraph.hdf5" % (majiq_config.outDir, majiq_config.exp_list[exp_idx])
    #         f_splicegraph = h5py.File(fname_sg, 'w', compression='gzip', compression_opts=9)
    #
    #         #as_table = f.create_group('LSVs')
    #         lsv_list.append(f)
    #         effective_readlen = (majiq_config.readLen - 16) + 1
    #         non_as_table = f.create_dataset("/const_junctions",
    #                                         (majiq_config.nrandom_junctions, effective_readlen),
    #                                         maxshape=(None, effective_readlen))
    #         junc_list.append(non_as_table)
    #
    #         splicegraph.append(f_splicegraph)
    #
    #     #logger.debug("... waiting childs")
    #     pool.close()
    #     while True:
    #         try:
    #             val = q.get(block=True, timeout=10)
    #             if val[0] == 0:
    #                 for exp_idx in majiq_config.tissue_repl[val[2]]:
    #                     val[1].to_hdf5(lsv_list[exp_idx])
    #             elif val[0] == 1:
    #                 for jdx, exp_idx in enumerate(majiq_config.tissue_repl[val[2]]):
    #                     if junc_idx[exp_idx] >= majiq_config.nrandom_junctions:
    #                         continue
    #                     junc_list[exp_idx][junc_idx[exp_idx], :] = val[1][jdx, :].toarray()
    #                     junc_idx[exp_idx] += 1
    #
    #             elif val[0] == 2:
    #                 val[1].to_hdf5(splicegraph[val[2]])
    #                 pass
    #             elif val[0] == 3:
    #                 lock_array[val[1]].release()
    #                 count += 1
    #
    #         except Queue.Empty:
    #             if count < majiq_config.num_final_chunks:
    #                 continue
    #             break
    #     pool.join()
    #
    #     for ff in file_list:
    #         ff.close()
    #
    #     for ff in splicegraph:
    #         ff.close()
    #
    #
    #     #GATHER
    # #    gather_files(majiq_config.outDir, '', params.gff_output, params.pcr_filename,
    # #                       nthreads=params.nthreads, logger=logger)
    #
    #     yappi.get_func_stats().print_all()

        #majiq_config.print_numbers()
        logger.info("End of execution")



# def majiq_builder(samfiles_list, chnk, pcr_validation=None, gff_output=None, create_tlb=True, only_rna=False,
#                   nondenovo=False, logging=None):
#
#     if not logging is None:
#         logging.debug("Building chunk %s" % chnk)
#
#     temp_dir = "%s/tmp/chunk_%s" % (majiq_config.outDir, chnk)
#     annot_file = '%s/annot_genes.pkl' % temp_dir
#     if not os.path.exists(annot_file):
#         majiq_builder.queue.put([3, chnk], block=True)
#         majiq_builder.lock_arr[chnk].acquire()
#         majiq_builder.lock_arr[chnk].release()
#         return
#
#     # temp_file = open(annot_file, 'rb')
#     gene_list = majiq.src.io_utils.load_bin_file(annot_file)
#
#     if create_tlb:
#         if not logging is None:
#             logging.debug("[%s] Recreatin Gene TLB" % chnk)
#         majiq_gene.recreate_gene_tlb(gene_list)
#
#     if not logging is None:
#         logging.info("[%s] Reading BAM files" % chnk)
#     majiq_io.read_sam_or_bam(samfiles_list, gene_list, chnk,
#                              nondenovo=nondenovo, logging=logging)
#     if not logging is None:
#         logging.info("[%s] Detecting intron retention events" % chnk)
#     majiq_io.rnaseq_intron_retention(samfiles_list, gene_list, chnk, permissive=majiq_config.permissive_ir,
#                                      nondenovo=nondenovo, logging=logging)
#     if not logging is None:
#         logging.info("[%s] Detecting LSV" % chnk)
#     const = analize.lsv_detection(gene_list, chnk, only_real_data=only_rna, out_queue=majiq_builder.queue,
#                                   logging=logging)
#
#
#     prepare_gc_content(gene_list, temp_dir)
#
#     # if pcr_validation:
#     #     utils.get_validated_pcr_lsv(lsv, temp_dir)
#     # if gff_output:
#     #     majiq_lsv.extract_gff(lsv, temp_dir)
#     majiq.src.voila_wrapper.generate_visualization_output(gene_list, temp_dir, majiq_builder.queue)
#     if not logging is None:
#         logging.info("[%s] Preparing output" % chnk)
#     lsv = None
#     utils.send_output(lsv, const, temp_dir, majiq_builder.queue, chnk, majiq_builder.lock_arr[chnk])


### TODO REMOVE FROM HERE

import argparse
def new_subparser():
    return argparse.ArgumentParser(add_help=False)

def main():
    """
    Main MAJIQ parser with all flags and subcommands
    """
    #REMINDER parser.add_parser(..... parents='[bla, ble]')
    parser = argparse.ArgumentParser(description="MAJIQ is a suite of tools for the analysis of Alternative "
                                                 "Splicing Events and Alternative Splicing Quantification.")

    #common flags (first ones are required)
    common = new_subparser()
    common.add_argument('--nthreads', default=4, type=int, help='Number of threads')
    common.add_argument('--nchunks', default=4, type=int, help='Number of chunks')
    common.add_argument('--tmp', default="/tmp/", help='Path to save the temporary files. [Default: %(default)s]')
    common.add_argument('--output', required=True, help='Path to save the pickle output to.')
    common.add_argument('--logger', default=None, help='Path for the logger. Default is output directory')
    common.add_argument('--silent', action='store_true', default=False, help='Silence the logger.')
    common.add_argument('--plotpath', default=None,
                        help='Path to save the plot to, if not provided will show on a matplotlib popup window')
    common.add_argument('--debug', type=int, default=0,
                        help="Activate this flag for debugging purposes, activates logger and jumps some "
                             "processing steps.")

    buildparser = new_subparser()
    buildparser.add_argument('transcripts', action="store", help='read file in SAM format')
    buildparser.add_argument('-conf', default=None, required=True, help='Provide study configuration file with all '
                                                                        'the execution information')
    buildparser.add_argument('--nogc', dest="gcnorm", action='store_false', default=True,
                             help='psianddelta GC content normalization [Default: GC content normalization activated]')
    buildparser.add_argument('--pcr', dest='pcr_filename', action="store", help='PCR bed file as gold_standard')
    buildparser.add_argument('--gff_output', dest='gff_output', default="lsvs.gff", action="store",
                             help='Filename where a gff with the lsv events will be generated')
    buildparser.add_argument('--min_denovo', default=2, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'denovo junction is real". '
                             '[Default: %(default)s]')
    buildparser.add_argument('--minreads', default=3, type=int,
                             help='Minimum number of reads threshold combining all positions in a LSV to consider that'
                                  'the LSV "exist in the data". '
                             '[Default: %(default)s]')
    buildparser.add_argument('--min_intronic_cov', default=1.5, type=float,
                             help='Minimum number of reads on average in intronic sites, only for intron retention.'
                                  'Default: %(default)s]')
    buildparser.add_argument('--num_chunks', default=-1, type=float,
                             help='Numbers of chunks the execution will be divided. That differs of nthread in the '
                                  'concurrency. Chunks is the total chunks of the execution, nthreads set how many of '
                                  'this chunks will be executed at the same time.')
    buildparser.add_argument('--minpos', default=2, type=int, help='Minimum number of start positions with at least 1 '
                                                                   'read in a LSV to consider that the LSV "exist in '
                                                                   'the data"')

    buildparser.add_argument('--only_rna', default=False, action='store_true', help='Use only rna detected junction in '
                                                                                    'order to detect LSV. If an exon '
                                                                                    'has only one junction with '
                                                                                    'coverage, it is not going to be '
                                                                                    'detected as an LSV')
    buildparser.add_argument('--non_denovo', default=False, action='store_true', help='Avoid denovo detection of '
                                                                                      'junction, splicesites and exons.'
                                                                                      ' This will speedup the execution'
                                                                                      ' but reduce the number of LSVs '
                                                                                      'detected')
    buildparser.add_argument('--only_gather', action='store_true', dest='onlygather', default=False)
    buildparser.add_argument('--permissive_ir', action='store_true', dest='permissive', default=False)


    #flags shared by calcpsi and deltapair

    subparsers = parser.add_subparsers(help='')

    parser_preprocess = subparsers.add_parser('build', help='Preprocess SAM/BAM files as preparation for the rest of '
                                                            'the tools (psi, deltapsi)', parents=[common, buildparser])
    parser_preprocess.set_defaults(func=build)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()


