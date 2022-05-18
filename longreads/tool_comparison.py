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

        self.extra_count_keys = ['TTT',
                                 'TTF', 
                                 'TFT', 
                                 'TFF',
                                 'FTT',
                                 'FTF',
                                 'FFT',
                                 'majiq_combination', 
                                 'majiq_novel', 
                                 'flair_combination', 
                                 'flair_novel',
                                 'flair_denovo_partial',
                                 'partial'
                                 ]
        self.counts = self._makeCountsObj()
        self.args = args


    def _makeCountsObj(self):
        tmpcounts = {}
        # for majiq in (True, False):
        #     for flair in (True, False):
        #         for annotated in (True, False):
        #             tmpcounts[countComp(majiq, flair, annotated)] = 0
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
        only_in_set2 = set()
        in_both_sets = set()
        superset = set1.union(set2)

        def compare(set1elem, set2elem):
            if set1elem == set2elem:
                return True
            if len(set1elem) == len(set2elem):
                for coords1, coords2 in zip(set1elem, set2elem):
                    startCondition = coords1[0] == -1 or coords2[0] == -1 or (abs(coords1[0] - coords2[0]) <= fuzziness)
                    endCondition = coords1[1] == -1 or coords2[1] == -1 or (abs(coords1[1] - coords2[1]) <= fuzziness)
                    if not startCondition or not endCondition:
                        break
                else:
                    return set2elem
                        #print('true by cond', coords1[0] == -1, coords2[0] == -1, abs(coords1[0] - coords2[0]) <= fuzziness, coords1[1] == -1, coords2[1] == -1, abs(coords1[1] - coords2[1]) <= fuzziness)
                        #return True
            return False

        for transcript in superset:
            result1 = False
            for t in set1:
                if compare(transcript, t):
                    if result1:
                        if self.args.verbose:
                            print("Warning, duplicate matches for transcript", transcript)
                    result1 = True
            result2 = False
            for t in set2:
                if compare(transcript, t):
                    if result2:
                        if self.args.verbose:
                            print("Warning, duplicate matches for transcript", transcript)
                    result2 = t

            if result1 and result2:
                in_both_sets.add(result2)
            elif result1:
                only_in_set1.add(transcript)
            else:
                only_in_set2.add(transcript)

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

    def add_partials(self, flair_result, annotated_starts, annotated_ends):

        total_partials = 0
        for transcript in flair_result:

            if (transcript[0].start not in annotated_starts) or (transcript[-1].end not in annotated_ends):
                self.counts['partial'] += 1
                total_partials += 1

        return total_partials

    # def add_flair_partials(self, transcript, annotated_starts, annotated_ends):

    #     total_partials = 0
    #     if (transcript[0].start not in annotated_starts) or (transcript[-1].end not in annotated_ends):
    #         self.counts['partial'] += 1
    #         total_partials += 1
    #     return total_partials


    def set_flair_unknown_ends(self, flair_result):
        _flair_result = set()
        for transcript in flair_result:
            _flair_result.add((
                exon(-1, transcript[0].end),
                *transcript[1:-1],
                exon(transcript[-1].start, -1)
            ))
        return _flair_result
    
    def current_coordinate(self, current_transcript):
        _all_coordinate = set()
        for exon in current_transcript:
            _all_coordinate.add(exon.start)
            _all_coordinate.add(exon.end)
        return _all_coordinate
    
    def all_annotated(self, in_flair_and_majiq, only_in_majiq, majiq_denovo):
        _annotated_coordinate = set()
        for transcript in in_flair_and_majiq.union(only_in_majiq):
            if not majiq_denovo[transcript]:
                for exon in transcript:       
                    _annotated_coordinate.add(exon.start)
                    _annotated_coordinate.add(exon.end)

        return _annotated_coordinate
                

    def add_data(self, majiq_result, majiq_denovo, majiq_has_reads, flair_result, annotated_starts, annotated_ends):
        """

        """
        tmpcounts = self._makeCountsObj()

        # before removing start / end information, we use it to check for partial isoforms



        # if self.args.fuzziness == 0:
        #     only_in_flair, only_in_majiq, in_flair_and_majiq = self.compare_exact(flair_result, majiq_result)
        # else:
        # print('F', flair_result)
        # print('M', majiq_result)
        only_in_flair, only_in_majiq, in_flair_and_majiq = self.compare_fuzzy(flair_result, majiq_result, self.args.fuzziness)
        # print(only_in_flair)
        # print(in_flair_and_majiq)
  
        _annotated_coordinate = self.all_annotated(in_flair_and_majiq, only_in_majiq, majiq_denovo)
        flair_partials = self.add_partials(only_in_flair, annotated_starts, annotated_ends)

        tmpcounts['flair_denovo_partial'] += flair_partials

        for transcript in in_flair_and_majiq:
            if majiq_denovo[transcript]:
                self.incCountPrint(tmpcounts, transcript, 'TTF')
            else:
                if majiq_has_reads[transcript]:
                    self.incCountPrint(tmpcounts, transcript, 'TTT')
                else:
                    self.incCountPrint(tmpcounts, transcript, 'FTT')


        for transcript in only_in_majiq:
            if majiq_denovo[transcript]:
                self.incCountPrint(tmpcounts, transcript, 'TFF')
                if not _annotated_coordinate.difference(self.current_coordinate(transcript)):
                    self.incCountPrint(tmpcounts, transcript, 'majiq_combination')            
                else:
                    self.incCountPrint(tmpcounts, transcript, 'majiq_novel')
            else:
                if not majiq_has_reads[transcript]:
                    self.incCountPrint(tmpcounts, transcript, 'FFT')
                else:
                    self.incCountPrint(tmpcounts, transcript, 'TFT')

        for transcript in only_in_flair:
            self.incCountPrint(tmpcounts, transcript, 'FTF')
            if not _annotated_coordinate.difference(self.current_coordinate(transcript)):
                self.incCountPrint(tmpcounts, transcript, 'flair_combination')
                if (transcript[0].start not in annotated_starts) or (transcript[-1].end not in annotated_ends):
                    self.incCountPrint(tmpcounts, transcript, 'flair_combination_partial')
            else:
                self.incCountPrint(tmpcounts, transcript, 'flair_novel')
                if (transcript[0].start not in annotated_starts) or (transcript[-1].end not in annotated_ends):
                    self.incCountPrint(tmpcounts, transcript, 'flair_combination_novel')
        
        for k, v in tmpcounts.items():
            self.counts[k] += v

        return tmpcounts

