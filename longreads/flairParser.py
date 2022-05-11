from pygtftk.gtf_interface import GTF
from graph import exon
from gtfparse import read_gtf


class FlairReader:

    def __init__(self, gtf_path):
        self.modules = {}
        self.df = read_gtf(gtf_path)

    def has_gene(self, gene_id):
        return len(self.df[self.df['gene_id'] == gene_id].index) != 0


    def gene(self, gene_id, extent=None, ignore_starts_ends=False):

        df_gene = self.df[self.df['gene_id'] == gene_id]

        found_transcripts = set()
        transcript_exons = []
        transcript_meta = {}



        def append_next(_transcript_exons):
            if _transcript_exons:
                #_transcript_exons = sorted(_transcript_exons, key=lambda e: e.start)
                if ignore_starts_ends:
                    _transcript_exons[0] = exon(-1, _transcript_exons[0].end)
                    _transcript_exons[-1] = exon(_transcript_exons[-1].start, -1)
                _transcript_exons = tuple(_transcript_exons)
                if not _transcript_exons in found_transcripts:
                    found_transcripts.add(_transcript_exons)
                    return tuple(_transcript_exons), transcript_meta

        for index, row in df_gene.iterrows():

            if row.feature == 'transcript':
                ret = append_next(transcript_exons)
                if ret:
                    yield ret

                transcript_exons = []
                transcript_meta = {x: getattr(row, x) for x in ('strand', 'gene_id', 'transcript_id')}
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