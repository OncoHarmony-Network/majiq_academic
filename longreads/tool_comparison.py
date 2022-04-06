import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import majiqv2, flairParser
from graph import exon
import pprint

majiq_splicegraph_path = '/slowdata/lrdata/majiq/splicegraph.sql'
majiq_gene_id="ENSG00000110492.17"
parser = majiqv2.MajiqV2Reader(majiq_splicegraph_path)
parser.parse_splicegraph(majiq_gene_id)

flair_gtf_path = '/slowdata/lrdata/flair/flair_filter_transcripts.gtf'
flair_gene_id = 'ENSG00000110492.17'
#
# save_path = '/tmp/lr_o'
# os.makedirs(save_path, exist_ok=True)



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

    def add_module_data(self, majiq_result, majiq_denovo, majiq_has_reads, flair_result):
        
        only_in_flair = flair_result.difference(majiq_result)
        only_in_majiq = majiq_result.difference(flair_result)
        in_flair_and_majiq = flair_result.intersection(majiq_result)

        for transcript in in_flair_and_majiq:
            if majiq_denovo[transcript]:
                self.counts[countComp(True, True, False)] += 1
            elif not majiq_has_reads[transcript]:
                self.counts[countComp(False, True, True)] += 1
            else:
                self.counts[countComp(True, True, True)] += 1
        for transcript in only_in_majiq:
            if majiq_denovo[transcript]:
                self.counts[countComp(True, False, False)] += 1
            elif not majiq_has_reads[transcript]:
                self.counts[countComp(False, False, True)] += 1
            else:
                self.counts[countComp(True, False, True)] += 1
        for transcript in only_in_flair:
            self.counts[countComp(False, True, False)] += 1


def compare_tools():
    for module_idx in range(parser.getNumModules()):

        majiq_module_extent = parser.moduleExtent(module_idx)

        """
        for in-module, by default the exons we receive from majiq start/end are technically not part of the module
        UNLESS, they have different start/end. As a simple way to deal with this, for matching purposes, we will trim all
        exon coordinates at the start/end of the module to match the module coordinates
        
        """

        tc = ToolComparer()
        flair_exons = set()
        ord_flair_exons = tuple(x[0] for x in flairreader.gene(flair_gene_id, extent=majiq_module_extent))
        for transcript in ord_flair_exons:
            flair_exons.add(tuple(exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in transcript))

        majiq_exons = set()
        majiq_denovo = {}
        majiq_has_reads = {}
        for (ord_majiq_transcript, majiq_meta, denovo, has_reads) in parser.getAllPaths(module_idx=module_idx):
            set_key = tuple(exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in ord_majiq_transcript)
            majiq_exons.add(set_key)
            majiq_denovo[set_key] = denovo
            majiq_has_reads[set_key] = has_reads


        tc.add_module_data(majiq_exons, majiq_denovo, majiq_has_reads, flair_exons)

        pprint.pprint(tc.counts)

def compare_all_tools():

    tc = ToolComparer()
    flair_exons = set()
    ord_flair_exons = tuple(x[0] for x in flairreader.gene(flair_gene_id, extent=majiq_module_extent))
    for transcript in ord_flair_exons:
        flair_exons.add(tuple(exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in transcript))

    majiq_exons = set()
    majiq_denovo = {}
    majiq_has_reads = {}
    for (ord_majiq_transcript, majiq_meta, denovo, has_reads) in parser.getAllPaths(module_idx=module_idx):
        set_key = tuple(exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in ord_majiq_transcript)
        majiq_exons.add(set_key)
        majiq_denovo[set_key] = denovo
        majiq_has_reads[set_key] = has_reads


    tc.add_module_data(majiq_exons, majiq_denovo, majiq_has_reads, flair_exons)

    pprint.pprint(tc.counts)

        #
        #
        # for transcript in ord_majiq_exons:
        #     majiq_exons.add(tuple(exon(max(majiq_module_extent[0], e.start), min(majiq_module_extent[1], e.end)) for e in transcript))
        #
        # print(ord_majiq_exons)
        #print(ord_majiq_transcripts)



if __name__ == "__main__":
    compare_tools()
