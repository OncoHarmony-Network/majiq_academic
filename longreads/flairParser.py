import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from graph import exon, junction
from graph import module as _module
from gtfparse import read_gtf


class FlairReader:

    def __init__(self, gtf_path):
        self.modules = {}
        self.df = read_gtf(gtf_path)

    def has_gene(self, gene_id):
        return len(self.df[self.df['gene_id'] == gene_id].index) != 0


    def get_exons(self, gene_id, majiq_module_extent=None, modules=False):
        flair_exons = set()
        ord_flair_exons = tuple(x for x in self.gene(gene_id, extent=majiq_module_extent, ignore_starts_ends=True))

        for transcript in ord_flair_exons:
            if modules:
                flair_exons.add(tuple(exon(max(majiq_module_extent[0], e.start) if e.start > 0 else e.start, min(majiq_module_extent[1], e.end) if e.end > 0 else e.end) for e in transcript))
            else:
                flair_exons.add(tuple(exon(e.start, e.end) for e in transcript))

        return flair_exons

    def extend_modules(self, module_extents, flair_transcripts):

        add_modules = []
        remove_modules = []
        for transcript in flair_transcripts:
            for i in range(len(transcript)-1):
                junc = junction(transcript[i].end, transcript[i + 1].start)

                start_in_module = None
                end_in_module = None
                for module in module_extents:
                    if junc.start >= module.start and junc.start <= module.end:
                        start_in_module = module
                    if junc.end >= module.start and junc.end <= module.end:
                        end_in_module = module

                if start_in_module and end_in_module and start_in_module != end_in_module:
                    remove_modules.append(start_in_module)
                    remove_modules.append(end_in_module)
                    add_modules.append(_module(start_in_module.start, end_in_module.end))

                if start_in_module and not end_in_module:
                    # extend to the next exon
                    remove_modules.append(start_in_module)
                    add_modules.append(_module(start_in_module.start, junc.end))

                if not start_in_module and end_in_module:
                    # extend to the previous exon
                    remove_modules.append(end_in_module)
                    add_modules.append(_module(junc.start, end_in_module.end))

        for module in remove_modules:
            if module in module_extents:
                module_extents.remove(module)
        for module in add_modules:
            if module not in module_extents:
                module_extents.append(module)

        return module_extents

    def gene(self, gene_id, extent=None, ignore_starts_ends=False):

        df_gene = self.df[self.df['gene_id'] == gene_id]

        found_transcripts = set()
        transcript_exons = []



        def append_next(_transcript_exons):
            if _transcript_exons:
                #_transcript_exons = sorted(_transcript_exons, key=lambda e: e.start)
                if ignore_starts_ends:
                    _transcript_exons[0] = exon(-_transcript_exons[0].start, _transcript_exons[0].end)
                    _transcript_exons[-1] = exon(_transcript_exons[-1].start, -_transcript_exons[-1].end)
                _transcript_exons = tuple(_transcript_exons)
                if not _transcript_exons in found_transcripts:
                    found_transcripts.add(_transcript_exons)
                    return tuple(_transcript_exons)

        for index, row in df_gene.iterrows():

            if row.feature == 'transcript':
                ret = append_next(transcript_exons)
                if ret:
                    yield ret

                transcript_exons = []
                #transcript_meta = {x: getattr(row, x) for x in ('strand', 'gene_id', 'transcript_id')}
                continue

            elif row.feature == 'exon':
                _exon = exon(int(row.start), int(row.end))


                if extent:
                    if (_exon.end >= extent[0] and _exon.end <= extent[1]) or \
                       (_exon.start >= extent[0] and _exon.start <= extent[1]):
                        transcript_exons.append(_exon)
                else:
                    transcript_exons.append(_exon)

        ret = append_next(transcript_exons)  # catch the last append
        if ret:
            yield ret


if __name__ == "__main__":

    gtf_path = '/slowdata/lrdata/flair/flair_filter_transcripts.gtf'

    flairreader = FlairReader(gtf_path)

    for exons, meta in flairreader.gene('ENSG00000138326.18'):
        print(exons)
        print(meta)
        break