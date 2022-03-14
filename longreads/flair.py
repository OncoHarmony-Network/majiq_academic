from pygtftk.gtf_interface import GTF
from graph import exon


class FlairReader:

    def __init__(self):
        self.modules = {}

    @classmethod
    def parse_gtf(cls, gtf_path, gene_id, extent=None):
        gtf = GTF(gtf_path)

        keys = ["start", "end", "gene_id", "transcript_id", "strand", "exon_number", "feature"]
        ki = {k: i for i, k in enumerate(keys)}

        found_transcripts = set()
        transcript_exons = set()
        transcript_meta = {}

        for row in gtf.extract_data(keys):
            if row[ki['gene_id']] != gene_id:
                continue

            if row[ki['feature']] == 'transcript':
                if transcript_exons:
                    transcript_exons = frozenset(transcript_exons)
                    if not transcript_exons in found_transcripts:
                        found_transcripts.add(transcript_exons)
                        yield tuple(transcript_exons), transcript_meta

                transcript_exons = []
                transcript_meta = {x: row[ki[x]] for x in ('strand', 'gene_id', 'transcript_id')}
                continue

            elif row[ki['feature']] == 'exon':
                _exon = exon(int(row[ki["start"]]), int(row[ki["end"]]))
                if extent:
                    if _exon.start > extent[0] or _exon.end < extent[1]:
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
    for exons, meta in FlairReader.parse_gtf(gtf_path, 'ENSG00000138326.18'):
        print(exons)
        print(meta)
        break