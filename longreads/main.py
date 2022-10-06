import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import majiqv2, flairParser
from tool_comparison import ToolComparer
from graph import exon
import pprint
from config import get_args
import csv
import traceback
import multiprocessing
import time
import os, sys
from multiprocessing import Manager, Pool
from graph import GeneNotFoundInSplicegraph
import glob

#majiq_splicegraph_path = '/slowdata/lrdata/majiq/splicegraph.sql'
#majiq_gene_id="gene:ENSG00000109534"


#flair_gtf_path = '/slowdata/lrdata/flair/flair_filter_transcripts.gtf'
#flair_gene_id = 'ENSG00000109534.16'
#
# save_path = '/tmp/lr_o'
# os.makedirs(save_path, exist_ok=True)

args = get_args()

majiq_splicegraph_path = args.majiq_splicegraph_path

majiqParser = majiqv2.MajiqV2Reader(majiq_splicegraph_path)
flair_gtf_path = args.flair_gtf_path


os.makedirs(args.output_path, exist_ok=True)

error_file_path = os.path.join(args.output_path, 'comparison.errors.txt')
if os.path.exists(error_file_path):
    os.remove(error_file_path)

if args.gene_ids_file:
    majiq_gene_ids = []
    flair_gene_ids = []
    with open(args.gene_ids_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                majiq_gene_ids.append(line)
                flair_gene_ids.append(line)

elif args.gene_id:

    if args.gene_id:
        majiq_gene_ids = [args.gene_id]
    else:
        majiq_gene_ids = majiqParser.gene_ids

    if args.gene_id_flair:
        flair_gene_ids = [args.gene_id_flair]
    else:
        flair_gene_ids = majiq_gene_ids

else:
    majiq_gene_ids = []
    flair_gene_ids = []
    for gene_id in majiqParser.gene_ids:
        majiq_gene_ids.append(gene_id)
        flair_gene_ids.append(gene_id)



print('~~~Parsing Flair~~~')
flairreader = flairParser.FlairReader(flair_gtf_path)
print('~~~Done Parsing Flair~~~')

def compare_gene(_args):
    """
    Single thread process self_contained

    The general strategy involved each process writing to it's own file, and then merging them at the end of the run
    """
    gene_id, modules, output_path_format, q = _args

    tc = ToolComparer(args)
    if q:
        pid = multiprocessing.current_process().pid
        tsv_path = output_path_format.format('.' + str(pid))
    else:
        tsv_path = output_path_format.format('')

    try:


        majiqParser.parse_splicegraph(gene_id)

        annotated_starts = majiqParser.annotated_starts(gene_id)
        annotated_ends = majiqParser.annotated_ends(gene_id)
        annotated_exons_starts, annotated_exons_ends = majiqParser.annotated_exons(gene_id)
        all_exons_starts, all_exons_ends = majiqParser.all_exons(gene_id)
        annotated_exons_order = majiqParser.annotated_exons_order(gene_id)

        #full_flair_exons = tuple(x for x in flairreader.gene(gene_id, extent=None, ignore_starts_ends=False))
        #gene_partial_count = tc.add_partials(full_flair_exons, annotated_starts, annotated_ends)
        modules_list = [majiqParser.moduleExtent(i) for i in range(majiqParser.getNumModules())] if modules else [None]
        modules_list = flairreader.extend_modules(modules_list, flairreader.get_exons(gene_id)) if modules else [None]

        for module_idx, module in enumerate(modules_list):
            if args.verbose >= 1:
                print(gene_id, '-------------------------')
                if modules:
                    print('module IDX', module_idx)



            """
            for in-module, by default the exons we receive from majiq start/end are technically not part of the module
            UNLESS, they have different start/end. As a simple way to deal with this, for matching purposes, we will trim all
            exon coordinates at the start/end of the module to match the module coordinates
            
            """

            flair_exons = flairreader.get_exons(gene_id, majiq_module_extent=module, modules=modules)

            majiq_exons, majiq_denovo, majiq_has_reads = majiqParser.allpaths_data(
                modules=modules,
                module_idx=module_idx if modules else None,
                max_paths=args.max_paths,
                majiq_module_extent=module
            )


            counts = tc.add_data(majiq_exons, majiq_denovo, majiq_has_reads, flair_exons, annotated_starts, annotated_ends, annotated_exons_starts, annotated_exons_ends, annotated_exons_starts, annotated_exons_ends, annotated_exons_order)
            #counts['partial'] = gene_partial_count

            row = [gene_id]
            if modules:
                row.append(module_idx)
            for count in counts.values():
                row.append(count)

            with open(tsv_path, 'a') as tsv:
                writer = csv.writer(tsv, delimiter='\t')
                writer.writerow(row)

            if args.verbose >= 1:
                pprint.pprint(tc.counts)
    except KeyboardInterrupt:
        raise
    except GeneNotFoundInSplicegraph:
        pass
    except RecursionError:
        if args.verbose >= 1:
            print("Recursion too great, gene skipped!", gene_id)
        with open(error_file_path, 'a') as f:
            f.write(traceback.format_exc() + '\n')
    except:
        print("Some error with gene!", gene_id)
        if args.debug:
            print(traceback.format_exc())
            return
        with open(error_file_path, 'a') as f:
            f.write(traceback.format_exc() + '\n')
    finally:
        if q:
            q.put(None)

def compare_tools(all_gene_ids, modules=False):

    work_size = len(all_gene_ids)

    main_tsv_path = os.path.join(args.output_path, 'comparison.tsv')
    output_path_format = main_tsv_path + "{0}"

    tc = ToolComparer(args)
    fieldnames = ['gene_id']
    if modules:
        fieldnames.append('module_idx')

    for key in tc.extra_count_keys:
        fieldnames.append(key)

    if args.debug or args.threads:
        with open(main_tsv_path, 'w') as tsv:
            writer = csv.writer(tsv, delimiter='\t')
            writer.writerow(fieldnames)


    if args.debug_num_genes and args.debug_num_genes > 0:
        all_gene_ids = all_gene_ids[:args.debug_num_genes]

    skipped_from_flair = 0
    # for i, gene_id in enumerate(all_gene_ids):
    #
    #     if not flairreader.has_gene(gene_id):
    #         skipped_from_flair += 1
    #         with open(error_file_path, 'a') as f:
    #             f.write(f"gene_id not found for flair, skipping: {gene_id}\n")
    #         continue
    #
    #     if not majiqParser.has_gene(gene_id):
    #         with open(error_file_path, 'a') as f:
    #             f.write(f"gene_id not found for majiq, skipping: {gene_id}\n")
    #         continue




    if args.threads == 1:
        try:
            for i, gene_id in enumerate(all_gene_ids):
                # t1 = time.time()
                # gene_id, modules, tc, output_path_format, q
                compare_gene((gene_id, modules, output_path_format, None))
                # t2 = time.time()
                # print(t2-t1)
                print('Processing Genes [%d/%d]\r' % (i, work_size), end="")
        except KeyboardInterrupt:
            print('                                                  \r', end="")


        print('                                                  \r', end="")

    else:
        old_read_files = glob.glob(args.output_path + "/comparison.tsv.*")
        for _f in old_read_files:
            os.remove(_f)


        manager = Manager()
        q = manager.Queue()

        p = Pool(args.threads)



        # voila_index = p.map(self._heterogen_pool_add_index, zip(lsv_ids, range(work_size), repeat(work_size)))
        classifier_pool = p.map_async(compare_gene, ((x, modules, output_path_format, q) for x in all_gene_ids),)

        # monitor loop
        while True:

            if classifier_pool.ready():
                break
            else:
                size = q.qsize()
                print('Processing Genes [%d/%d]\r' % (size, work_size), end="")
                time.sleep(2)

        print('                                                  \r', end="")
        res = classifier_pool.get()



        #print('                                                  \r', end="")
        #print(f'{round((skipped_from_flair / (i+1)) * 100.0, 1)}% of genes were missing from flair and skipped')


    if not (args.threads == 1):
        read_files = glob.glob(args.output_path + "/comparison.tsv.*")

        with open(main_tsv_path, "wb") as outfile:
            outfile.write(('\t'.join(fieldnames) + '\n').encode())
            for f in read_files:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
                os.remove(f)

    if not modules and args.threads == 1:
        final_counts_path = os.path.join(args.output_path, 'comparison.totals.txt')
        with open(final_counts_path, 'w') as f:
            line = ''
            for majiq in ('T', 'F'):
                for flair in ('T', 'F'):
                    for annotated in ('T', 'F'):
                        line += f'{majiq}{flair}{annotated}\t'
            for key in tc.extra_count_keys:
                line += f'{key}\t'
            line += '\n'
            f.write(line)

            line = ''
            for count in tc.counts.values():
                line += f'{count}\t'
            line += '\n'
            f.write(line)



if __name__ == "__main__":

    print('~~~Running comparison~~~')
    compare_tools(majiq_gene_ids, modules=args.per_module)
    print('~~~Finished Running comparison~~~')
