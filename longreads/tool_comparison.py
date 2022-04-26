import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import majiqv2, flairParser
from graph import exon
import pprint
from config import get_args
import csv
import traceback


from collections import namedtuple
countComp = namedtuple('ComparisonCount', 'majiq flair annotated')

class ToolComparer:

    def __init__(self, args):
        """
        Here we gather counts for each permutation of tools used + "gene annotation", of which we consider majiq-non-denovo
        to be an authoritative source.
        """

        self.extra_count_keys = ['partial']
        self.counts = self._makeCountsObj()
        self.args = args


    def _makeCountsObj(self):
        tmpcounts = {}
        for majiq in (True, False):
            for flair in (True, False):
                for annotated in (True, False):
                    tmpcounts[countComp(majiq, flair, annotated)] = 0
        for key in self.extra_count_keys:
            tmpcounts[key] = 0
        return tmpcounts

    def incCountPrint(self, tmpcounts, transcript, key):
        tmpcounts[key] += 1
        if self.args.verbose >= 2:
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

    def _is_partial_isoform(self, majiq_transcripts, flair_transcript, fuzziness=0):
        for majiq_transcript in majiq_transcripts:
            for flair_exon in flair_transcript:
                for majiq_exon in majiq_transcript:
                    if abs(flair_exon.start - majiq_exon.start) <= fuzziness and abs(flair_exon.end - majiq_exon.end) <= fuzziness:
                        break
                else:
                    # for this flair exon, we went through all of the majiq exons, but could not find a match, we need
                    # to try the next majiq transcript
                    break

            else:
                # we got through all of the flair exons, and each one found a match in this majiq transcript
                return True

            # we broke somewhere, indicating we need to repeat the process
            continue

        return False

    def add_data(self, majiq_result, majiq_denovo, majiq_has_reads, flair_result):

        if self.args.fuzziness == 0:
            only_in_flair, only_in_majiq, in_flair_and_majiq = self.compare_exact(flair_result, majiq_result)
        else:
            only_in_flair, only_in_majiq, in_flair_and_majiq = self.compare_fuzzy(flair_result, majiq_result, self.args.fuzziness)

        tmpcounts = self._makeCountsObj()

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
            if self._is_partial_isoform(only_in_majiq, transcript, self.args.fuzziness):
                self.incCountPrint(tmpcounts, transcript, 'partial')
            else:
                self.incCountPrint(tmpcounts, transcript, countComp(False, True, False))

        for k, v in tmpcounts.items():
            self.counts[k] += v

        return tmpcounts

