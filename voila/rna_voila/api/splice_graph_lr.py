
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
    'sl': '#44AA99',
    'l': '#EA07C5',
    'la': '#DDCC77',
    's': '#117733',
    'sa': '#882255',
    'ao': '#808080'
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

    def _get_transcripts(self, gene_id):
        if gene_id in self.lrdb:
            return self.lrdb[gene_id]['transcripts']
        return []

    def gene(self, gene_id, annotated_exons):
        ret = []

        # I'm not sure why it happens that we get exons like this but it seems to once in a while...
        annotated_exons = [x for x in annotated_exons if x[0] != x[1]]

        for transcript in self._get_transcripts(gene_id):
            d = {
                'id': transcript['id'],
                'reads': transcript['transcript_reads'],
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
                        "presence": "l"
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
                        "presence": "l"
                    }
                    for j in transcript['intron_retention']
                ],
                'intron_retention_reads': {
                    transcript['experiment']: {j[0]: {j[1]: r} for j, r in zip(transcript['intron_retention'], transcript['intron_retention_reads'])}
                },
                'exon_reads': {
                    transcript['experiment']: {j[0]: {j[1]: r} for j, r in zip(transcript['exons'], transcript['exon_reads'])}
                }
            }



            _annotated_exons = IntervalTree.from_tuples(annotated_exons)

            for lr_exon in transcript['exons']:
                matching_annotated = _annotated_exons.overlap(lr_exon[0], lr_exon[1])
                ex_d = {'color': combined_colors['ao'], 'presence': 'la'}
                if matching_annotated:
                    matching_annotated = matching_annotated.pop()

                    ex_d['start'] = lr_exon[0]
                    ex_d['end'] = lr_exon[1]
                    ex_d['annotated'] = 1
                    ex_d['annotated_start'] = matching_annotated[0]
                    ex_d['annotated_end'] = matching_annotated[1]
                    ex_d['ext_color'] = combined_colors['l']
                    _annotated_exons.remove(matching_annotated)
                else:
                    ex_d['start'] = lr_exon[0]
                    ex_d['end'] = lr_exon[1]
                    ex_d['annotated'] = 0
                    ex_d['annotated_start'] = lr_exon[0]
                    ex_d['annotated_end'] = lr_exon[1]
                    ex_d['color'] = combined_colors['l']
                    ex_d['presence'] = 'l'

                d['exons'].append(ex_d)

            for annot_exon in _annotated_exons:
                d['exons'].append({
                    'start': annot_exon.begin,
                    'end': annot_exon.end,
                    'annotated': 1,
                    'annotated_start': annot_exon.begin,
                    'annotated_end': annot_exon.end,
                    'color': 'hidden',
                    'presence': 'ao'
                })
            d['exons'].sort(key=lambda x: x['start'])

            exon_number = 1
            for exon in (d['exons'] if transcript['strand'] == "+" else reversed(d['exons'])):
                if exon['presence'] == 'l':
                    exon['number'] = ''
                else:
                    exon['number'] = exon_number
                    exon_number += 1

            d['start'] = d['exons'][0]['start']
            d['end'] = d['exons'][-1]['end']
            d['strand'] = transcript['strand']


            ret.append(d)
        return ret

    def _overlap_categories(self, gene_id, shortread, subkey):

        sr_junctions = set((j['start'], j['end']) for j in shortread[subkey] if j['annotated'] == 0)
        annot_sr_junctions = set((j['start'], j['end']) for j in shortread[subkey] if j['annotated'] == 1 and j['has_reads'] == 1)
        annot_only_junctions = set((j['start'], j['end']) for j in shortread[subkey] if j['annotated'] == 1 and j['has_reads'] == 0)
        annot_junctions = set((j['start'], j['end']) for j in shortread[subkey] if j['annotated'] == 1)
        lr_junctions = set()
        for transcript in self._get_transcripts(gene_id):
            for j in transcript[subkey]:
                lr_junctions.add((j[0], j[1]))

        """
        This is the main logic used for determine which of the six categories the junction or intron falls into
        """
        j_sla = annot_sr_junctions & lr_junctions
        j_l = lr_junctions - annot_junctions - sr_junctions
        j_sl = sr_junctions & (lr_junctions - annot_junctions)
        j_la = (lr_junctions & annot_junctions) - j_sla
        j_s = sr_junctions - j_sla - j_l - j_sl
        j_sa = annot_sr_junctions - j_sla
        j_ao = annot_only_junctions - j_la

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

    def _even_strsize(self, _sr_reads, _lr_reads):
        """
        If there are less characters in the lr or sr, we even them so that the separator appears centered.
        Currently, in voila, the SVG read counts font is sans-sarif, not monospace. Therefore, this is only
        an approximate solution which uses two non-breaking spaces to achieve the desired effect.
        """
        sr_letter_count = len(str(_sr_reads))
        lr_letter_count = len(str(_lr_reads))
        max_lettercount = max(sr_letter_count, lr_letter_count)
        _sr_reads = '\u00A0\u00A0' * (max_lettercount - sr_letter_count) + str(_sr_reads)
        _lr_reads = str(_lr_reads) + '\u00A0\u00A0' * (max_lettercount - lr_letter_count)

        return _sr_reads, _lr_reads

    def _combine_summary(self, gene_id, shortread, subkey):

        if subkey == 'junctions':
            readssubkey = 'junction_reads'
        else:
            readssubkey = 'intron_retention_reads'

        sr_reads = {exp:v for exp, v in shortread[readssubkey].items()} #  if exp.endswith('Combined')
        lr_reads = {v['experiment']: {(j[0], j[1],): r for j, r in zip(v[subkey], v[readssubkey])} for v in self._get_transcripts(gene_id)}
        j_sla, j_l, j_sl, j_la, j_s, j_sa, j_ao = self._overlap_categories(gene_id, shortread, subkey)
        #self._debugprint(j_sla, j_l, j_sl, j_la, j_s, j_sa, j_ao)

        shortread[subkey] = []
        shortread[readssubkey] = {'combined':{}}


        for _type, juncset in zip(('sla', 'l', 'sl', 'la', 's', 'sa', 'ao'), (j_sla, j_l, j_sl, j_la, j_s, j_sa, j_ao)):
            for junc in juncset:
                shortread[subkey].append({
                    "annotated": 1,
                    "color": combined_colors[_type],
                    "presence": _type,
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
                    reads = self._format_reads(_sr_reads, _lr_reads)
                    if reads:
                        shortread[readssubkey]['combined'][junc[0]] = {junc[1]: reads}
                else:
                    if junc[1] in shortread[readssubkey]['combined'][junc[0]]:
                        print('error: there seem to be overlapping junctions between the six categories!', junc)
                        # pass
                        assert False
                        #shortread[readssubkey]['combined'][junc[0]][junc[1]] += final_reads
                    else:
                        reads = self._format_reads(_sr_reads, _lr_reads)
                        if reads:
                            shortread[readssubkey]['combined'][junc[0]][junc[1]] = reads

        return shortread

    def _format_reads(self, _sr_reads, _lr_reads):
        if _sr_reads or _lr_reads:
            #_sr_reads, _lr_reads = self._even_strsize(_sr_reads, _lr_reads)
            #return f"{_sr_reads}╦{_lr_reads}"
            return (_sr_reads, _lr_reads,)
        return None

    def _overlaps(self, s1, e1, s2, e2):
        return s1 <= e2 and s2 <= e1

    def _combine_exons(self, gene_id, shortread):
        """
        This will extend the exons with extensions Note: in the rare case that both SR and LR extend the same exon by
        different amounts, we don't really have a good way to show this yet
        """


        #print(shortread['exons'])
        for transcript in self._get_transcripts(gene_id):

            for lr_exon in transcript['exons']:

                found_lr_exon = False
                for sr_exon in shortread['exons']:
                #print(transcript)

                    #print(lr_exon, sr_exon, self._overlaps(lr_exon[0], lr_exon[1], sr_exon['annotated_start'], sr_exon['annotated_end']))
                    if self._overlaps(lr_exon[0], lr_exon[1], sr_exon['annotated_start'], sr_exon['annotated_end']):

                        found_lr_exon = True
                        if lr_exon[0] < sr_exon['start']:
                            sr_exon['start'] = lr_exon[0]
                            sr_exon['ext_color'] = combined_colors['l']
                        if lr_exon[1] > sr_exon['end']:
                            sr_exon['end'] = lr_exon[1]
                            sr_exon['ext_color'] = combined_colors['l']
                        if not sr_exon['color']:  # dubious
                            sr_exon['color'] = combined_colors['l']

                if not found_lr_exon:
                    shortread['exons'] = list(shortread['exons'])
                    shortread['exons'].append({
                        'start': lr_exon[0], 'end': lr_exon[1], 'annotated': 0, 'color': combined_colors['l'], 'annotated_start': lr_exon[0], 'annotated_end': lr_exon[1]
                    })
                    shortread['exons'] = sorted(shortread['exons'], key=lambda d: d['start'])
                    shortread['start'] = min(shortread['start'], lr_exon[0])
                    shortread['end'] = max(shortread['end'], lr_exon[1])


        return shortread

    def combined_gene(self, gene_id, shortread):

        #long_reads_only = not shortread['junction_reads'] and not shortread['intron_retention_reads']

        shortread = self._combine_summary(gene_id, shortread, 'junctions')
        shortread = self._combine_summary(gene_id, shortread, 'intron_retention')
        shortread = self._combine_exons(gene_id, shortread)



        return shortread

    def lsvs(self, gene_id):
        if gene_id in self.lrdb:
            return self.lrdb[gene_id]['lsvs']
        return {}

    def lsv(self, gene_id, lsv_id):
        return self.lsvs(gene_id)[lsv_id]

    def has_lsv(self, gene_id, lsv_id):
        return lsv_id in self.lsvs(gene_id)

    def combined_junctions(self, gene_id, lsv_junctions):
        """
        Figure out LR junction LSVs that match with short read labeled LSVs
        This is the old version of on demand PSI calculation that will probably not be used anymore
        """
        lr_equiv_reads = []
        total_reads = 0

        all_lr_transcripts = self._get_transcripts(gene_id)


        for seek_junc in lsv_junctions:
            for transcript in all_lr_transcripts:
                for j, r in zip(transcript['junctions'], transcript['junction_reads']):
                    if j[0] == seek_junc[0] and j[1] == seek_junc[1]:
                        total_reads += r
                        lr_equiv_reads.append(r)
                        break
                else:
                    continue
                break
            else:
                lr_equiv_reads.append(0)

        assert len(lr_equiv_reads) == len(lsv_junctions)
        if total_reads > 0:
            lr_equiv_reads = [r / total_reads for r in lr_equiv_reads]
        return lr_equiv_reads
