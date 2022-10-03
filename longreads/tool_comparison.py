import os

import matplotlib.pyplot as plt
import matplotlib.patches as patches

import majiqv2, flairParser
from graph import exon, junction
import pprint
from config import get_args
import csv
import traceback
from itertools import product, combinations


from collections import namedtuple
countComp = namedtuple('ComparisonCount', 'majiq flair annotated')

class ToolComparer:

    def __init__(self, args):
        """
        Here we gather counts for each permutation of tools used + "gene annotation", of which we consider majiq-non-denovo
        to be an authoritative source.

        majiq_combination
        majiq_novel
        flair_only_combination
        flair_only_novel
        flair_combination_novel
        flair_only_partial
        flair_combination_partial
        flair_novel_partial
        flair_partial_combination_novel

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
                                 'flair_only_combination',
                                 'flair_novel',
                                 'flair_combination_novel',
                                 'flair_combination_partial',
                                 'flair_novel_partial',
                                 'flair_partial_combination_novel',
                                 'partial',
                                 'flair_novel_alt3',
                                 'flair_novel_alt5',
                                 'flair_novel_exon',
                                 'flair_novel_intron',
                                 'flair_novel_alt3_alt5',
                                 'flair_novel_alt3_exon',
                                 'flair_novel_alt3_intron',
                                 'flair_novel_alt5_exon',
                                 'flair_novel_alt5_intron',
                                 'flair_novel_intron_exon',
                                 'flair_novel_alt3_alt5_exon',
                                 'flair_novel_alt3_alt5_intron',
                                 'flair_novel_alt3_intron_exon',
                                 'flair_novel_alt5_intron_exon',
                                 'flair_novel_alt3_alt5_intron_exon',
                                 'flair_FTF_unknown'
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

    

    def compare_fuzzy(self, set1, set2, fuzziness_5, fuzziness_3):
        """
        Return "only in set1", "only in set2" and "in both sets" by fuzzy matching
        To be considered a match, the length of the element must be the same, and also each inner integer value
        must be within N (fuzziness) absolute value of the other set
        note in the fuzzy match case, the value of SET1 will be used in the return
        """
        only_in_set1 = set()
        only_in_set2 = set()
        in_both_sets = set()
        superset = set1.union(set2)

        def fuzzy_distance(flair_transcript, majiq_transcript):
            
            total_distance = 0
            print(len(majiq_transcript), len(flair_transcript))
            print("M_before ",majiq_transcript)
            print("F_before ",flair_transcript)

            if len(flair_transcript) <= len(majiq_transcript):
                for i in range(len(majiq_transcript)-1):
                    junc_majiq = junction(majiq_transcript[i].end, majiq_transcript[i + 1].start)
                    if junc_majiq.start > abs(flair_transcript[0].start) and junc_majiq.end < flair_transcript[0].end:
                        return 0
                    elif junc_majiq.start > flair_transcript[-1].start and junc_majiq.end < abs(flair_transcript[-1].end):
                        return 0

                for k in range(len(majiq_transcript) - len(flair_transcript)):
                    majiq_transcript = majiq_transcript[k:k+len(flair_transcript)]
                    
                    for coords1, coords2 in zip(flair_transcript, majiq_transcript):
                        startCondition = coords1[0] <= -2 or coords2[0] <= -2 or (abs(coords1[0] - coords2[0]) <= fuzziness_5)
                        endCondition = coords1[1] <= -2 or coords2[1] <= -2 or (abs(coords1[1] - coords2[1]) <= fuzziness_3)
                        print("start ",startCondition)
                        print("end ", endCondition)
                        if not startCondition or not endCondition:
                            break
                    else:
                        print("M ",majiq_transcript)
                        print("F ",flair_transcript)
                        print("M_start: ", abs(majiq_transcript[k].start))
                        print("M_end: ", abs(majiq_transcript[k].end))
                        print("F_start: ", abs(flair_transcript[k].start))
                        print("F_end: ", abs(flair_transcript[k].end))
                        dist5 = abs(abs(majiq_transcript[k].start) - abs(flair_transcript[k].start))
                        print("dist5 ",dist5)
                        dist3 = abs(abs(majiq_transcript[k].end) - abs(flair_transcript[k].end))
                        print("dist3 ",dist3)

                        total_distance += dist5 + dist3
                    # if (dist5 > fuzziness_5) or (dist3 > fuzziness_3):
                    #     return False

                #print("TOTAL DISTANCE ",total_distance)

            return total_distance


        def compare(set1elem, set2elem):
        
            if len(set1elem) == 1:
                # we can not compare a one-exon flair transcript with majiq in our current paradigm, because both
                # the start and end coordinate will be TSS/TES ; which majiq does not measure automatic mark as match
                return True, 0
            if set1elem == set2elem:
                return True, 0
            if len(set1elem) <= len(set2elem):
                # set1elem: Flair, set2elem: MAJIQ
                for i in range(len(set2elem)-1):
                    junc_majiq = junction(set2elem[i].end, set2elem[i + 1].start)
                    if junc_majiq.start > abs(set1elem[0].start) and junc_majiq.end < set1elem[0].end:
                        return False, 0
                    elif junc_majiq.start > set1elem[-1].start and junc_majiq.end < abs(set1elem[-1].end):
                        return False, 0
                for k in range(len(set2elem) - len(set1elem)+1):
                    slide_set2 = set2elem[k:k+len(set1elem)]

                    for coords1, coords2 in zip(set1elem, slide_set2):
                        startCondition = coords1[0] <= -2 or coords2[0] <= -2 or (abs(coords1[0] - coords2[0]) <= fuzziness_5)
                        endCondition = coords1[1] <= -2 or coords2[1] <= -2 or (abs(coords1[1] - coords2[1]) <= fuzziness_3)
                        #print("start ",startCondition)
                        #print("end ", endCondition)
                        if not startCondition or not endCondition:
                            break
                    else:
                        total_distance = fuzzy_distance(set1elem, set2elem)
                        # print("why ",total_distance)
                        # print("why2 ",set2elem)
                        return set2elem, total_distance

                        #print('true by cond', coords1[0] == -1, coords2[0] == -1, abs(coords1[0] - coords2[0]) <= fuzziness, coords1[1] == -1, coords2[1] == -1, abs(coords1[1] - coords2[1]) <= fuzziness)
                        #return True
            return False, 0

        
        only_in_set2 = set2.copy()
        for f_transcript in set1:
            same = False
            for m_transcript in set2:
                matched, total_distance = compare(f_transcript, m_transcript)
                # print("matched ",matched)
                # print("dist", total_distance)
                # print("m_trans ",m_transcript)
                if matched:
                    same = True
                    in_both_sets.add((f_transcript, m_transcript, total_distance))
    
                    if m_transcript in only_in_set2:
                        only_in_set2.remove(m_transcript)
                 
            if not same:
                only_in_set1.add(f_transcript)

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
                self.counts['flair_novel_partial'] += 1
                total_partials += 1

        return total_partials

    def set_flair_unknown_ends(self, flair_result):
        _flair_result = set()
        for transcript in flair_result:
            _flair_result.add((
                exon(-transcript[0].start, transcript[0].end),
                *transcript[1:-1],
                exon(transcript[-1].start, -transcript[-1].end)
            ))
        return _flair_result
    
    def current_coordinate(self, current_transcript):
        _all_coordinate = set()
        for exon in current_transcript:
            if exon.start > 0:
                _all_coordinate.add(exon.start)
            if exon.end > 0:
                _all_coordinate.add(exon.end)
        return _all_coordinate

    def get_junctions(self, transcript):
        junctions = set()
        for i in range(len(transcript)-1):
            junctions.add(junction(transcript[i].end, transcript[i+1].start))
            # print("See the process",junctions)
        return junctions
    
    def all_annotated(self, in_flair_and_majiq, only_in_majiq, majiq_denovo):
        """
        This finds in annotated as dictated by non denovo paths only, generally superseded be directly using annotated exons from the splicegraph
        """
        _annotated_coordinate = set()
        for transcript in in_flair_and_majiq.union(only_in_majiq):
            if not majiq_denovo[transcript]:
                for exon in transcript:       
                    _annotated_coordinate.add(exon.start)
                    _annotated_coordinate.add(exon.end)

        return _annotated_coordinate

    def substring_FTF(self, **kwargs):
        """
        kwargs: partial, novel, combination
        """
        if all(kwargs.get(x, False) for x in ('partial', 'novel', 'combination')):
            return 'flair_partial_combination_novel'
        elif all(kwargs.get(x, False) for x in ('partial', 'novel')):
            return 'flair_novel_partial'
        elif all(kwargs.get(x, False) for x in ('partial', 'combination')):
            return 'flair_combination_partial'
        elif all(kwargs.get(x, False) for x in ('novel', 'combination')):
            return 'flair_combination_novel'
        elif kwargs.get('novel', False):
            return 'flair_novel'
        elif kwargs.get('combination', False):
            return 'flair_only_combination'
        else:
            #print("unexpected", kwargs)
            assert False

    def substring_FTF_novel(self, **kwargs):
        """
        kwargs: novel_alt3, novel_alt5, novel_intron, novel_exon
        """
        novel_substring = 'flair_novel'
        novel_substring += "_alt3" if kwargs.get('novel_alt3',False) else ""
        novel_substring += "_alt5" if kwargs.get('novel_alt5',False) else ""
        novel_substring += "_intron" if kwargs.get('novel_intron',False) else ""
        novel_substring += "_exon" if kwargs.get('novel_exon',False) else ""

        assert novel_substring != 'flair_novel' 

        return novel_substring 


    def add_data(self, majiq_result, majiq_denovo, majiq_has_reads, flair_result, annotated_starts, annotated_ends, all_exons_starts, all_exons_ends, annotated_exons_starts, annotated_exons_ends, annotated_exons_order):
        """

        """
        tmpcounts = self._makeCountsObj()
        known_junctions = set()
        annotated_exon_coords = set(annotated_exons_starts).union(set(annotated_exons_ends))

        # before removing start / end information, we use it to check for partial isoforms

        only_in_flair, only_in_majiq, in_flair_and_majiq = self.compare_fuzzy(flair_result, majiq_result, self.args.fuzziness5, self.args.fuzziness3)
        # print("Only_flair",only_in_flair)
        # print("Only_majiq",only_in_majiq)
        for f_transcript, m_transcript, total_distance in in_flair_and_majiq:
        # for f_transcript, m_transcript in in_flair_and_majiq:
            #print("total :",total_distance)
            known_junctions = known_junctions.union(self.get_junctions(m_transcript))
            #known_junctions = known_junctions.union(self.get_junctions(f_transcript))
            if majiq_denovo[m_transcript]:
                self.incCountPrint(tmpcounts, m_transcript, 'TTF')
            else:
                if majiq_has_reads[m_transcript]:
                    self.incCountPrint(tmpcounts, m_transcript, 'TTT')
                else:
                    self.incCountPrint(tmpcounts, m_transcript, 'FTT')
            if (-f_transcript[0].start not in annotated_starts) or (-f_transcript[-1].end not in annotated_ends):
                self.incCountPrint(tmpcounts, f_transcript, 'partial')


        for transcript in only_in_majiq:
            known_junctions = known_junctions.union(self.get_junctions(transcript))
            if majiq_denovo[transcript]:
                self.incCountPrint(tmpcounts, transcript, 'TFF')
                if self.current_coordinate(transcript).issubset(annotated_exon_coords):
                    self.incCountPrint(tmpcounts, transcript, 'majiq_combination')            
                else:
                    self.incCountPrint(tmpcounts, transcript, 'majiq_novel')
            else:
                if not majiq_has_reads[transcript]:
                    self.incCountPrint(tmpcounts, transcript, 'FFT')
                else:
                    self.incCountPrint(tmpcounts, transcript, 'TFT')

        # fun debugging help things
        # print('F', flair_result)
        # print('M', majiq_result)
        # print("A_ex", annotated_exon_coords)
        # print("Known_j", known_junctions)
        # print(' | f_o', only_in_flair)
        # print('flair transcript number:', len(only_in_flair))
        # print('m_o', only_in_majiq)
        # print('f_m_b', in_flair_and_majiq)

        for transcript in only_in_flair:                                     
            self.incCountPrint(tmpcounts, transcript, 'FTF')
            partial, novel, combination = False, False, False
            novel_alt3, novel_alt5, novel_intron, novel_exon = False, False, False, False

            flair_new_exon = set()

            # check for flair exons in between annotated exons
            for flair_exon in transcript[1:-1]:
                for i in range(len(annotated_exons_order)-1):
                    E1, E2 = annotated_exons_order[i], annotated_exons_order[i+1]
                    if flair_exon.start > E1.end and flair_exon.end < E2.start:
                        # print("E1 start ",E1.start)
                        # print("E2 end ",E2.end)
                        flair_new_exon.add(flair_exon.start)
                        flair_new_exon.add(flair_exon.end)
                        #print("new exon: ", flair_new_exon)
                        novel_exon = True


            # check for flair exons existing before/after any annotated exons
            for flair_exon in transcript:
                #print(flair_exon)
                if all(abs(flair_exon.end) < e.start for e in annotated_exons_order) or \
                   all(abs(flair_exon.start) > e.end for e in annotated_exons_order):
                    flair_new_exon.add(flair_exon.start)
                    flair_new_exon.add(flair_exon.end)
                    novel = True
                    novel_exon = True


            for i in range(len(transcript)-1):
                junc = junction(transcript[i].end, transcript[i + 1].start)
                # print(known_junctions)
                if {junc.start, junc.end}.issubset(annotated_exon_coords) and junc not in known_junctions:
                    combination = True
                elif junc not in known_junctions:
                    novel = True
                    if junc.end not in annotated_exon_coords and junc.end not in flair_new_exon:
                        novel_alt5 = True
                    if junc.start not in annotated_exon_coords and junc.start not in flair_new_exon:
                        novel_alt3 = True

            for majiq_exon_1, majiq_exon_2 in combinations(zip(all_exons_starts, all_exons_ends), 2):
                # print("ann_start ", annotated_exons_starts)
                # print("ann_end ", annotated_exons_ends)
                # print("1 ",majiq_exon_1[1])
                # print("2 ",majiq_exon_2[0])
                junc_start = majiq_exon_1[1]
                junc_end = majiq_exon_2[0]
                for flair_exon in transcript:
                    # print("flair ", flair_exon)
                    if abs(flair_exon.start) <= junc_start and abs(flair_exon.end) > junc_end:
                        # print(abs(flair_exon.start), junc_start, abs(flair_exon.end), junc_end)
                        # print("flair start ",flair_exon.start)
                        # print("flair end ",flair_exon.end)
                        # print("junc start ",junction_.start)
                        # print("junc end ",junction_.end)
                        novel_intron = True
                        novel = True
                    else:
                        continue
                    break
            
            assert not (not novel and novel_exon)

            if novel:
                novel_name = self.substring_FTF_novel(novel_alt3=novel_alt3, novel_alt5=novel_alt5, novel_intron=novel_intron, novel_exon=novel_exon)
                self.incCountPrint(tmpcounts, transcript, novel_name)      

            if (-transcript[0].start not in annotated_starts) or (-transcript[-1].end not in annotated_ends):
                partial = True


            try:
                name = self.substring_FTF(partial=partial, novel=novel, combination=combination)
            except:
                #print("Error with substring_FTF", partial, novel, combination, transcript)
                #raise
                name = 'flair_FTF_unknown'

            self.incCountPrint(tmpcounts, transcript, name)

        
        for k, v in tmpcounts.items():
            self.counts[k] += v

        return tmpcounts

