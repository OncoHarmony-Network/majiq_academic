from pygtftk.gtf_interface import GTF
from graph import exon
from gtfparse import read_gtf


class FlairReader:

    def __init__(self, gtf_path):
        self.modules = {}
        self.df = read_gtf(gtf_path)

    def has_gene(self, gene_id):
        return len(self.df[self.df['gene_id'] == gene_id].index) != 0

    def gene(self, gene_id, extent=None):

        df_gene = self.df[self.df['gene_id'] == gene_id]

        found_transcripts = set()
        transcript_exons = set()
        transcript_meta = {}

        for index, row in df_gene.iterrows():

            if row.feature == 'transcript':
                if transcript_exons:
                    transcript_exons = frozenset(transcript_exons)
                    if not transcript_exons in found_transcripts:
                        found_transcripts.add(transcript_exons)
                        yield tuple(transcript_exons), transcript_meta

                transcript_exons = []
                transcript_meta = {x: getattr(row, x) for x in ('strand', 'gene_id', 'transcript_id')}
                continue

            elif row.feature == 'exon':
                _exon = exon(int(row.start), int(row.end))
                if extent:
                    if (_exon.end > extent[0] and _exon.end < extent[1]) or \
                       (_exon.start > extent[0] and _exon.start < extent[1]):
                        transcript_exons.append(_exon)
                else:
                    transcript_exons.append(_exon)

        if transcript_exons:
            transcript_exons = frozenset(transcript_exons)
            if not transcript_exons in found_transcripts:
                found_transcripts.add(transcript_exons)
                yield tuple(transcript_exons), transcript_meta


if __name__ == "__main__":

    gtf_path = '/slowdata/lrdata/flair/flair_filter_transcripts.gtf'

    flairreader = FlairReader(gtf_path)

    for exons, meta in flairreader.gene('ENSG00000138326.18'):
        print(exons)
        print(meta)
        break