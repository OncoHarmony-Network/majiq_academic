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
    # __eq__ = lambda self, other: (self.strand == other.strand
    #                               and self.exons[0].coords[0] < other.exons[-1].coords[1] and self.exons[-1].coords[1] > other.exons[0].coords[0])
    # __ne__ = lambda self, other: (self.strand != other.strand
    #                               or self.exons[0].coords[0] >= other.exons[-1].coords[1] or self.exons[-1].coords[1] <= other.exons[0].coords[0])
    # __lt__ = lambda self, other: ((self.exons[-1].coords[1] < other.exons[0].coords[0]
    #                                        or (self.exons[-1].coords[1] > other.exons[0].coords[0] and self.exons[0].coords[0] < other.exons[-1].coords[1]
    #                                            and self.strand == '+' and other.strand == '-')))
    # __gt__ = lambda self, other: ((self.exons[0].coords[0] > other.exons[-1].coords[1]
    #                                        or (self.exons[-1].coords[1] > other.exons[0].coords[0] and self.exons[0].coords[0] < other.exons[-1].coords[1]
    #                                            and self.strand == '-' and other.strand == '+')))

    def __init__(self, name, strand, exons, junctions, chrom=None):
        self.name = name
        self.strand = strand
        self.exons = exons
        self.junctions = junctions
        self.chrom = chrom
        self.start = exons[0].coords[0]
        self.end = exons[-1].coords[1]

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
    #
    # def to_bed12(self):
    #     #  "%(a)s, %(a)s" % {'a':'test'}
    #     bed_fields = {
    #         'chrom': self.chrom,
    #         'chromStart': self.start,
    #         'chromEnd': self.end,
    #         'name': None,
    #         'score': 0,
    #         'strand': self.strand,
    #         'thickStart': self.start,
    #         'thickEnd': self.end,
    #         'itemRgb': None,
    #         'blockCount': len(self.exons),
    #         'blockSizes'
    #     }
    #     pass
        # chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        # chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        # chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
        # The 9 additional optional BED fields are:
        #
        # name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        # score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
        # shade
        # score in range  	 166	167-277	278-388	389-499	500-611	612-722	723-833	834-944	 945
        # strand - Defines the strand - either '+' or '-'.
        # thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
        # thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
        # itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
        # blockCount - The number of blocks (exons) in the BED line.
        # blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
        # blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.