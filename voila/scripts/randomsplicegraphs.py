import random
import uuid
from itertools import izip

from voila.splice_graphics import ExonGraphic, JunctionGraphic, GeneGraphic


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)] * n)


class RandomSpliceGraph(object):
    def __init__(self, exons=None, junctions=None, experiments=None, gene_id=None):
        """
        Generate random exons and random junctions for creating a GeneGraphic object.
        :param exons: number of exons
        :param junctions: number junctions
        :param experiments: number of experiments
        :param gene_id: gene graphics gene id
        """
        self.gene_id = gene_id
        if not gene_id:
            self.gene_id = uuid.uuid4().hex

        self.experiments = experiments
        if not experiments:
            self.experiments = random.randrange(0, 10)

        self.junctions = junctions
        if not junctions:
            self.junctions = random.randrange(1, 20)

        self.exons = exons
        if not exons:
            self.exons = random.randrange(2, 20)

    def exons_generator(self, number, experiments):
        """
        Generate random exons.
        :param number: number of exons
        :param experiments: number of experiments
        :return: list of exons
        """
        starts_ends = grouped(sorted(random.sample(range(1000), number * 2)), 2)
        types = grouped((random.randrange(0, 3) for _ in range(number * experiments)), experiments)
        return [ExonGraphic(None, None, start=start, end=end, exon_type_list=exon_type_list)
                for (start, end), exon_type_list in zip(starts_ends, types)]

    def junction_generator(self, exons, number, experiments):
        """
        Generate random junctions.
        :param exons: list of exons
        :param number: number of junctions
        :param experiments: number of experiments
        :return: list of junctions
        """

        exon_pairs = (sorted(exon_pair) for exon_pair in
                      grouped((random.randrange(0, len(exons)) for _ in range(number * 2)), 2))
        types = grouped((random.randrange(0, 3) for _ in range(number * experiments)), experiments)
        reads = grouped((random.randrange(0, 500) for _ in range(number * experiments)), experiments)
        clean_reads = grouped((random.randrange(0, 100) for _ in range(number * experiments)), experiments)

        junctions = []
        for (start_exon, end_exon), junc, reads, clean_reads in zip(exon_pairs, types, reads, clean_reads):
            coords = sorted((random.randrange(exons[start_exon].start, exons[start_exon].end),
                             random.randrange(exons[end_exon].start, exons[end_exon].end)))

            junction = JunctionGraphic(start=coords[0],
                                       end=coords[1],
                                       junction_type_list=junc,
                                       reads_list=reads,
                                       clean_reads_list=clean_reads)
            junctions.append(junction)

        return junctions

    def get_gene_graphic(self):
        exons = self.exons_generator(number=self.exons, experiments=self.experiments)
        junctions = self.junction_generator(number=self.junctions, exons=exons, experiments=self.experiments)
        return GeneGraphic(self.gene_id, exons=exons, junctions=junctions)
