import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import majiqv2, flairParser
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


IN_ANNOTATION = True
IN_FLAIR = True
IN_MAJIQ = True

from collections import namedtuple
countComp = namedtuple('ComparisonCount', 'majiq flair annotated')

class ToolComparer:

    def __init__(self):
        """
        Here we gather counts for each permutation of tools used + "gene annotation", of which we consider majiq-non-denovo
        to be an authoritative source.
        """

        self.counts = {}
        for majiq in (True, False):
            for flair in (True, False):
                for annotated in (True, False):
                    self.counts[countComp(majiq, flair, annotated)] = 0

    def incCountPrint(self, tmpcounts, transcript, key):
        tmpcounts[key] += 1
        if args.verbose >= 2:
            print("PATH", key, transcript)

    def compare_fuzzy(self, set1, set2, fuzziness):
        """
        Return "only in set1", "only in set2" and "in both sets" by fuzzy matching
        To be considered a match, the length of the element must be the same, and also each inner integer value
        must be within N (fuzziness) absolute value of the other set
        note in the fuzzy match case, the value of SET2 will be used in the return
        """
        only_in_set1 = set()
        in_both_sets = set()

        for set1elem in set1:
            for set2elem in set2:
                if set1elem == set2elem:
                    in_both_sets.add(set2elem)
                    break
                if len(set1elem) == len(set2elem):
                    for coords1, coords2 in zip(set1elem, set2elem):
                        if abs(coords1[0] - coords2[0]) <= fuzziness and abs(coords1[1] - coords2[1]) <= fuzziness:
                            in_both_sets.add(set2elem)
                            break
                    else:
                        continue
                    break
            else:
                only_in_set1.add(set1elem)

        only_in_set2 = set2.difference(in_both_sets)

        return only_in_set1, only_in_set2, in_both_sets

    def compare_exact(self, set1, set2):
        """
        Return "only in set1", "only in set2" and "in both sets"  by exact matching the elements of sets
        """
        only_in_set1 = set1.difference(set2)
        only_in_set2 = set2.difference(set1)
        in_both_sets = set1.intersection(set2)
        return only_in_set1, only_in_set2, in_both_sets

    def add_data(self, majiq_result, majiq_denovo, majiq_has_reads, flair_result):

        if args.fuzziness == 0:
            only_in_flair, only_in_majiq, in_flair_and_majiq = self.compare_exact(flair_result, majiq_result)
        else:
            only_in_flair, only_in_majiq, in_flair_and_majiq = self.compare_fuzzy(flair_result, majiq_result, args.fuzziness)

        tmpcounts = {}
        for majiq in (True, False):
            for flair in (True, False):
                for annotated in (True, False):
                    tmpcounts[countComp(majiq, flair, annotated)] = 0

        for transcript in in_flair_and_majiq:
            if majiq_denovo[transcript]:
                self.incCountPrint(tmpcounts, transcript, countComp(True, True, False))
            elif not majiq_has_reads[transcript]:
                self.incCountPrint(tmpcounts, transcript, countComp(False, True, True))
            else:
                self.incCountPrint(tmpcounts, transcript, countComp(True, True, True))
        for transcript in only_in_majiq:
            if majiq_denovo[transcript]:
                self.incCountPrint(tmpcounts, transcript, countComp(True, False, False))
            elif not majiq_has_reads[transcript]:
                self.incCountPrint(tmpcounts, transcript, countComp(False, False, True))
            else:
                self.incCountPrint(tmpcounts, transcript, countComp(True, False, True))
        for transcript in only_in_flair:
            self.incCountPrint(tmpcounts, transcript, countComp(False, True, False))


        for k, v in tmpcounts.items():
            self.counts[k] += v

        return tmpcounts


def compare_tools(modules=False):

    work_size = len(majiq_gene_ids)

    main_tsv_path = os.path.join(args.output_path, 'comparison.tsv')
    tc = ToolComparer()
    with open(main_tsv_path, 'w') as tsv:
        fieldnames = ['gene_id']
        if modules:
            fieldnames.append('module_idx')
        for majiq in ('T', 'F'):
            for flair in ('T', 'F'):
                for annotated in ('T', 'F'):
                    fieldnames.append(f'{majiq}{flair}{annotated}')
        writer = csv.writer(tsv, delimiter='\t')
        writer.writerow(fieldnames)

        for i, (majiq_gene_id, flair_gene_id) in enumerate(zip(majiq_gene_ids, flair_gene_ids)):

            print('Processing Genes [%d/%d]\r' % (i, work_size), end="")

            try:
                if not flairreader.has_gene(flair_gene_id):
                    with open(error_file_path, 'a') as f:
                        f.write(f"gene_id not found for flair, skipping: {flair_gene_id}")
                    continue

                if not majiqParser.has_gene(majiq_gene_id):
                    with open(error_file_path, 'a') as f:
                        f.write(f"gene_id not found for majiq, skipping: {majiq_gene_id}")
                    continue

                majiqParser.parse_splicegraph(majiq_gene_id)

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
                    ord_flair_exons = tuple(x[0] for x in flairreader.gene(flair_gene_id, extent=majiq_module_extent))
                    for transcript in ord_flair_exons:
                        if modules:
                            flair_exons.add(tuple(exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in transcript))
                        else:
                            flair_exons.add(tuple(exon(e.start, e.end) for e in transcript))

                    majiq_exons = set()
                    majiq_denovo = {}
                    majiq_has_reads = {}
                    for (ord_majiq_transcript, majiq_meta, denovo, has_reads) in majiqParser.getAllPaths(module_idx=module_idx if modules else None):
                        if modules:
                            set_key = tuple(exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in ord_majiq_transcript)
                        else:
                            set_key = tuple(exon(e.start, e.end) for e in ord_majiq_transcript)
                        majiq_exons.add(set_key)
                        majiq_denovo[set_key] = denovo
                        majiq_has_reads[set_key] = has_reads

                    counts = tc.add_data(majiq_exons, majiq_denovo, majiq_has_reads, flair_exons)

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
            except:
                print("Some error with gene!", majiq_gene_id, flair_gene_id)
                with open(error_file_path, 'a') as f:
                    f.write(traceback.format_exc())

        print('                                                  \r', end="")


    if not modules:
        final_counts_path = os.path.join(args.output_path, 'comparison.totals.txt')
        with open(final_counts_path, 'w') as f:
            line = ''
            for majiq in ('T', 'F'):
                for flair in ('T', 'F'):
                    for annotated in ('T', 'F'):
                        line += f'{majiq}{flair}{annotated}\t'
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