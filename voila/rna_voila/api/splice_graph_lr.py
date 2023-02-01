
import pickle
from intervaltree import Interval, IntervalTree
from statistics import median
from math import ceil

"""
       majiq + annotation
       majiq denovo
       lr + annotation
       lr denovo
       majiq + lr
       majiq + lr + denovo 
       """
combined_colors = {
    'sla': '#332288',
    'sl': '#88CCEE',
    'l': '#AA4499',
    'la': '#CC6677',
    's': '#44AA99',
    'sa': '#882255',
    'ao': 'grey'
}

lrdb = None

class SpliceGraphLR:
    def __init__(self, filename):
        global lrdb
        if lrdb is None:
            with open(filename, 'rb') as f:
                self.lrdb = pickle.load(f)
            lrdb = self.lrdb
        else:
            self.lrdb = lrdb

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def gene(self, gene_id, annotated_exons):
        ret = []


        for transcript in self.lrdb.get(gene_id, []):
            d = {
                'experiment': transcript['experiment'],
                'exons': [
                ],
                'junctions': [
                    {
                        "annotated": 1,
                        "color": combined_colors['l'],
                        "end": j[1],
                        "has_reads": 1,
                        "is_constitutive": 0,
                        "is_simplified": 0,
                        "start": j[0],
                    }
                    for j in transcript['junctions']
                ],
                'junction_reads': {
                    transcript['experiment']: {j[0]: {j[1]: r} for j, r in zip(transcript['junctions'], transcript['junction_reads'])}
                },
                'intron_retention': [
                    {
                        "annotated": 1,
                        "color": combined_colors['l'],
                        "end": j[1],
                        "has_reads": 1,
                        "is_constitutive": 0,
                        "is_simplified": 0,
                        "start": j[0],
                    }
                    for j in transcript['intron_retention']
                ],
                'intron_retention_reads': {
                    transcript['experiment']: {j[0]: {j[1]: r} for j, r in zip(transcript['intron_retention'], transcript['intron_retention_reads'])}
                }
            }
            annotated_exons = IntervalTree.from_tuples(annotated_exons)
            for lr_exon in transcript['exons']:
                matching_annotated = annotated_exons.overlap(lr_exon[0], lr_exon[1])
                ex_d = {'color': 'grey'}
                if matching_annotated:
                    matching_annotated = matching_annotated.pop()

                    ex_d['start'] = lr_exon[0]
                    ex_d['end'] = lr_exon[1]
                    ex_d['annotated'] = 1
                    ex_d['annotated_start'] = matching_annotated[0]
                    ex_d['annotated_end'] = matching_annotated[1]
                    ex_d['ext_color'] = 'orange'
                    annotated_exons.remove(matching_annotated)
                else:
                    ex_d['start'] = lr_exon[0]
                    ex_d['end'] = lr_exon[1]
                    ex_d['annotated'] = 0
                    ex_d['annotated_start'] = lr_exon[0]
                    ex_d['annotated_end'] = lr_exon[1]

                d['exons'].append(ex_d)
            for annot_exon in annotated_exons:
                d['exons'].append({
                    'start': annot_exon.begin,
                    'end': annot_exon.end,
                    'annotated': 1,
                    'annotated_start': annot_exon.begin,
                    'annotated_end': annot_exon.end,
                    'color': 'hidden'
                })
            d['exons'].sort(key=lambda x: x['start'])


            ret.append(d)
        return ret

    def _overlap_categories(self, gene_id, shortread, subkey):

        sr_junctions = set((j['start'], j['end']) for j in shortread[subkey] if j['annotated'] == 0)
        annot_sr_junctions = set((j['start'], j['end']) for j in shortread[subkey] if j['annotated'] == 1 and j['has_reads'] == 1)
        annot_only_junctions = set((j['start'], j['end']) for j in shortread[subkey] if j['annotated'] == 1 and j['has_reads'] == 0)
        annot_junctions = set((j['start'], j['end']) for j in shortread[subkey] if j['annotated'] == 1)
        lr_junctions = set()
        for transcript in self.lrdb.get(gene_id, []):
            for j in transcript[subkey]:
                lr_junctions.add((j[0], j[1]))

        j_sla = annot_sr_junctions & lr_junctions
        j_l = lr_junctions - annot_junctions
        j_sl = sr_junctions & j_l
        j_la = (lr_junctions & annot_junctions) - j_sla
        j_s = sr_junctions
        j_sa = annot_sr_junctions - j_sla
        j_ao = annot_only_junctions - j_la

        #assert bool(j_sla & j_l & j_sl & j_la & j_s & j_sa & j_ao) is False  # there should be no overlaps
        return j_sla, j_l, j_sl, j_la, j_s, j_sa, j_ao

    def _debugprint(self, j_sla, j_l, j_sl, j_la, j_s, j_sa, j_ao):
        import pprint
        print('--------------------------------')
        pprint.pprint({
            'short long annotation': j_sla,
            'long only': j_l,
            'short and long': j_sl,
            'long and annotation': j_la,
            'short only': j_s,
            'short and annotation': j_sa,
            'annotation only': j_ao
        })
        print('--------------------------------')

    def _combine_summary(self, gene_id, shortread, subkey):

        if subkey == 'junctions':
            readssubkey = 'junction_reads'
        else:
            readssubkey = 'intron_retention_reads'

        sr_reads = {exp:v for exp, v in shortread[readssubkey].items()} #  if exp.endswith('Combined')
        lr_reads = {v['experiment']: {(j[0], j[1],): r for j, r in zip(v[subkey], v[readssubkey])} for v in self.lrdb.get(gene_id, [])}
        j_sla, j_l, j_sl, j_la, j_s, j_sa, j_ao = self._overlap_categories(gene_id, shortread, subkey)
        #self._debugprint(j_sla, j_l, j_sl, j_la, j_s, j_sa, j_ao)

        shortread[subkey] = []
        shortread[readssubkey] = {'combined':{}}


        for _type, juncset in zip(('sla', 'l', 'sl', 'la', 's', 'sa', 'ao'), (j_sla, j_l, j_sl, j_la, j_s, j_sa, j_ao)):
            for junc in juncset:
                shortread[subkey].append({
                    "annotated": 1,
                    "color": combined_colors[_type],
                    "end": junc[1],
                    #"gene_id": "ENSMUSG00000031134",
                    "has_reads": 1,
                    "is_constitutive": 0,
                    "is_simplified": 0,
                    "start": junc[0],
                })

                _sr_reads = []
                _lr_reads = []

                if _type in ('sla', 'sl', 's', 'sa'):
                    for v in sr_reads.values():
                        reads = v.get(junc[0], {}).get(junc[1], 0)
                        if reads:
                            _sr_reads.append(reads)
                if _type in ('sla', 'l', 'sl', 'la'):
                    for v in lr_reads.values():
                        reads = v.get((junc[0], junc[1],), 0)
                        if reads:
                            _lr_reads.append(reads)

                _sr_reads = ceil(median(_sr_reads)) if _sr_reads else 0
                _lr_reads = ceil(median(_lr_reads)) if _lr_reads else 0
                # all_reads = _sr_reads + _lr_reads
                # final_reads = ceil(median(all_reads)) if all_reads else 0
                if junc[0] not in shortread[readssubkey]['combined']:
                    shortread[readssubkey]['combined'][junc[0]] = {junc[1]: (_sr_reads, _lr_reads,)}
                else:
                    if junc[1] in shortread[readssubkey]['combined'][junc[0]]:
                        print('hmmm?', junc)
                        # pass
                        assert False
                        #shortread[readssubkey]['combined'][junc[0]][junc[1]] += final_reads
                    else:
                        shortread[readssubkey]['combined'][junc[0]][junc[1]] = (_sr_reads, _lr_reads,)

        return shortread

    def _overlaps(self, s1, e1, s2, e2):
        return s1 <= e2 and s2 <= e1

    def _combine_exons(self, gene_id, shortread):
        """
        This will extend the exons with extensions Note: in the rare case that both SR and LR extend the same exon by
        different amounts, we don't really have a good way to show this yet
        """

        #print(shortread['exons'])
        for sr_exon in shortread['exons']:

            for transcript in self.lrdb.get(gene_id, []):
                #print(transcript)
                for lr_exon in transcript['exons']:
                    #print(lr_exon, sr_exon, self._overlaps(lr_exon[0], lr_exon[1], sr_exon['annotated_start'], sr_exon['annotated_end']))
                    if self._overlaps(lr_exon[0], lr_exon[1], sr_exon['annotated_start'], sr_exon['annotated_end']):
                        if lr_exon[0] < sr_exon['start']:
                            sr_exon['start'] = lr_exon[0]
                            sr_exon['ext_color'] = 'orange'
                        if lr_exon[1] > sr_exon['end']:
                            sr_exon['end'] = lr_exon[1]
                            sr_exon['ext_color'] = 'orange'

        return shortread


    def combined_gene(self, gene_id, shortread):

        shortread = self._combine_summary(gene_id, shortread, 'junctions')
        shortread = self._combine_summary(gene_id, shortread, 'intron_retention')
        shortread = self._combine_exons(gene_id, shortread)

        return shortread