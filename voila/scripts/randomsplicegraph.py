import random
import uuid
from itertools import izip

from voila.splice_graphics import ExonGraphic, JunctionGraphic, GeneGraphic


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)] * n)


class RandomSpliceGraph(object):
    def __init__(self, experiments=None, exons=None, junctions=None, gene_id=None):
        """
        Generate random exons and random junctions for creating a GeneGraphic object. The developer admits this isn't
        the cleanest or nicest piece of software.
        :param exons: number of exons
        :param junctions: number junctions
        :param experiments: number of experiments
        :param gene_id: gene graphics gene id
        """
        self.gene_id = gene_id
        if not gene_id:
            self.gene_id = uuid.uuid4().hex

        self.junctions = junctions
        if junctions is None:
            self.junctions = random.randrange(0, 10)

        self.exons = exons
        if exons is None:
            self.exons = random.randrange(1, 30)

        self.experiments = experiments
        if experiments is None:
            self.experiments = 2

        self.gene_max_length = 2500
        self.exon_length_multiplier = 2
        self.junction_start_stop_iterations = 100

    def get_gene_graphic(self):
        exons = self.exons_generator()
        junctions = self.junctions_generator(exons=exons)
        return GeneGraphic(self.gene_id, exons=exons, junctions=junctions)

    def exons_generator(self):
        """
        Generate random exons.
        :return: list of exons
        """
        # get a list of random unique integers and group them into twos
        starts_ends = grouped(
            sorted(random.sample(range(self.gene_max_length / self.exon_length_multiplier), self.exons * 2)),
            2
        )

        # list of exon types
        exons_types = grouped((random.randrange(0, 3) for _ in range(self.exons * self.experiments)), self.experiments)

        return [ExonGraphic(a3=None,
                            a5=None,
                            start=start * self.exon_length_multiplier,
                            end=end * self.exon_length_multiplier,
                            exon_type_list=exon_type_list)
                for (start, end), exon_type_list in zip(starts_ends, exons_types)]

    def junctions_generator(self, exons):
        """
        Generate random junctions.
        :param exons: list of exons
        :return: list of junctions
        """
        junctions = []

        # randomly choose two exons
        exon_pairs = (
            sorted(exon_pair)
            for exon_pair in grouped((random.randrange(0, len(exons))
                                      for _ in range(self.junctions * 2)), 2)
        )

        # list of junction types
        junc_types = grouped((random.randrange(0, 3) for _ in range(self.junctions * self.experiments)),
                             self.experiments)

        # list of junction reads
        reads = grouped((random.randrange(0, 2000) for _ in range(self.junctions * self.experiments)), self.experiments)

        for (start_exon, end_exon), t, r in zip(exon_pairs, junc_types, reads):
            start = 0
            stop = 0
            count = 0

            # because you can start and end on the same exon, randomly choose a start and stop that aren't equal
            while start == stop:
                start, stop = sorted((random.randrange(exons[start_exon].start, exons[start_exon].end),
                                      random.randrange(exons[end_exon].start, exons[end_exon].end)))
                count += 1
                if count > self.junction_start_stop_iterations:
                    raise Exception('Yeah, we can\'t decide on a start and stop for this junction...')

            junction = JunctionGraphic(start=start,
                                       end=stop,
                                       junction_type_list=t,
                                       reads_list=r)
            junctions.append(junction)

        return junctions
