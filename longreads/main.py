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



def compare_tools(modules=False):

    work_size = len(majiq_gene_ids)

    main_tsv_path = os.path.join(args.output_path, 'comparison.tsv')
    tc = ToolComparer(args)
    with open(main_tsv_path, 'w') as tsv:
        fieldnames = ['gene_id']
        if modules:
            fieldnames.append('module_idx')
        for majiq in ('T', 'F'):
            for flair in ('T', 'F'):
                for annotated in ('T', 'F'):
                    fieldnames.append(f'{majiq}{flair}{annotated}')
        for key in tc.extra_count_keys:
            fieldnames.append(key)
        writer = csv.writer(tsv, delimiter='\t')
        writer.writerow(fieldnames)

        skipped_from_flair = 0
        for i, (majiq_gene_id, flair_gene_id) in enumerate(zip(majiq_gene_ids, flair_gene_ids)):

            if args.debug_num_genes and i > args.debug_num_genes:
                break

            print('Processing Genes [%d/%d] (%s)\r' % (i, work_size, majiq_gene_id), end="")


            try:
                if not flairreader.has_gene(flair_gene_id):
                    skipped_from_flair += 1
                    with open(error_file_path, 'a') as f:
                        f.write(f"gene_id not found for flair, skipping: {flair_gene_id}\n")
                    continue

                if not majiqParser.has_gene(majiq_gene_id):
                    with open(error_file_path, 'a') as f:
                        f.write(f"gene_id not found for majiq, skipping: {majiq_gene_id}\n")
                    continue

                majiqParser.parse_splicegraph(majiq_gene_id)

                annotated_starts = majiqParser.annotated_starts(majiq_gene_id)
                annotated_ends = majiqParser.annotated_ends(majiq_gene_id)
                full_flair_exons = tuple(x[0] for x in flairreader.gene(flair_gene_id, extent=None, ignore_starts_ends=False))
                gene_partial_count = tc.add_partials(full_flair_exons, annotated_starts, annotated_ends)


                for module_idx in range(majiqParser.getNumModules() if modules else 1):
                    if args.verbose >= 1:
                        print(majiq_gene_id, '-------------------------')
                        if modules:
                            print('module IDX', module_idx)




                    majiq_module_extent = majiqParser.moduleExtent(module_idx) if modules else None

                    """
                    for in-module, by default the exons we receive from majiq start/end are technically not part of the module
                    UNLESS, they have different start/end. As a simple way to deal with this, for matching purposes, we will trim all
                    exon coordinates at the start/end of the module to match the module coordinates
                    
                    """



                    flair_exons = set()
                    ord_flair_exons = tuple(x[0] for x in flairreader.gene(flair_gene_id, extent=majiq_module_extent, ignore_starts_ends=True))

                    for transcript in ord_flair_exons:
                        if modules:
                            flair_exons.add(tuple(exon(max(majiq_module_extent[0], e.start) if e.start != -1 else -1, min(majiq_module_extent[1], e.end) if e.end != -1 else -1) for e in transcript))
                        else:
                            flair_exons.add(tuple(exon(e.start, e.end) for e in transcript))

                    majiq_exons = set()
                    majiq_denovo = {}
                    majiq_has_reads = {}

                    num_paths = 0
                    for (ord_majiq_transcript, majiq_meta, denovo, has_reads) in majiqParser.getAllPaths(module_idx=module_idx if modules else None):
                        num_paths += 1
                        if args.max_paths == 0 or num_paths > args.max_paths:
                            raise RecursionError()
                        if modules:
                            set_key = tuple(exon(max(majiq_module_extent[0], e.start) if e.start != -1 else -1, min(majiq_module_extent[1], e.end) if e.end != -1 else -1) for e in ord_majiq_transcript)
                        else:
                            set_key = tuple(exon(e.start, e.end) for e in ord_majiq_transcript)
                        majiq_exons.add(set_key)
                        majiq_denovo[set_key] = denovo
                        majiq_has_reads[set_key] = has_reads


                    counts = tc.add_data(majiq_exons, majiq_denovo, majiq_has_reads, flair_exons)
                    counts['partial'] = gene_partial_count

                    row = [majiq_gene_id]
                    if modules:
                        row.append(module_idx)
                    for count in counts.values():
                        row.append(count)
                    writer.writerow(row)

                    if args.verbose >= 1:
                        pprint.pprint(tc.counts)
            except KeyboardInterrupt:
                print('                                                  \r', end="")
                break
            except RecursionError:
                print("Recursion too great, gene skipped!", majiq_gene_id, flair_gene_id)
                with open(error_file_path, 'a') as f:
                    f.write(traceback.format_exc() + '\n')
            except:
                print("Some error with gene!", majiq_gene_id, flair_gene_id)
                if args.debug:
                    print(traceback.format_exc())
                    break
                with open(error_file_path, 'a') as f:
                    f.write(traceback.format_exc() + '\n')

        print('                                                  \r', end="")
        print(f'{round((skipped_from_flair / i) * 100.0, 1)}% of genes were missing from flair and skipped' )


    if not modules:
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
    compare_tools(modules=args.per_module)
    print('~~~Finished Running comparison~~~')
