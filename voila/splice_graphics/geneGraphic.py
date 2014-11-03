import json

__author__ = 'abarrera'


class GeneGraphic(object):

    __eq__ = lambda self, other: (self.chrom == other.chrom and self.strand == other.strand
                                  and self.start < other.end and self.end > other.start)
    __ne__ = lambda self, other: (self.chrom != other.chrom or self.strand != other.strand
                                  or self.start >= other.end or self.end <= other.start)
    __lt__ = lambda self, other: (self.chrom < other.chrom
                                  or (self.chrom == other.chrom
                                      and (self.end < other.start
                                           or (self.end > other.start and self.start < other.end
                                               and self.strand == '+' and other.strand == '-'))))
    __gt__ = lambda self, other: (self.chrom > other.chrom
                                  or (self.chrom == other.chrom
                                      and (self.start > other.end
                                           or (self.end > other.start and self.start < other.end
                                               and self.strand == '-' and other.strand == '+'))))

    def __init__(self, name, strand, exons, junctions, chrom=None, id=None):
        self.id = id
        self.name = name
        self.strand = strand
        self.exons = exons
        self.junctions = junctions
        self.chrom = chrom
        self.start = exons[0].coords[0]
        self.end = exons[-1].coords[1]

    def get_id(self):
        return self.id

    def get_name(self):
        return self.name

    def get_strand(self):
        return self.strand

    def get_exons(self):
        return self.exons

    def get_junctions(self):
        return self.junctions

    def get_chrom(self):
        return self.chrom

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_coords(self):
        # return [self.start, self.end]
        # if self.strand == '+':
        return [self.exons[0].coords[0], self.exons[-1].coords[1]]
        # if self.strand == '-':
        #     return [self.exons[-1].coords[1], self.exons[0].coords[0]]

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def to_bed12(self):
        #  "%(a)s, %(a)s" % {'a':'test'}
        #fields = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
        bed_fields = [
            self.chrom,
            self.start,
            self.end,
            'lsv',
            0,
            self.strand,
            self.start,
            self.end,
            0,
            len(self.exons),
            ','.join([str(abs(e.get_coords()[0] - e.get_coords()[1])) for e in self.exons]), #[::int(self.strand + '1')]
            ','.join([str(abs(e.get_coords()[0] - self.start)) for e in self.exons]),
        ]
        return "\t".join([str(b) for b in bed_fields])