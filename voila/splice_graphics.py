import json
import os
from multiprocessing import Manager, Pool
from multiprocessing.process import Process
from multiprocessing.queues import JoinableQueue

import h5py

from voila import constants
from voila.constants import PROCESS_COUNT
from voila.hdf5 import HDF5
from voila.utils.voilaLog import voilaLog


class GeneGraphic(HDF5):
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

    def __init__(self, id, name=None, strand=None, exons=list(), junctions=list(), chrom=None):
        super(GeneGraphic, self).__init__()

        self.id = id
        self.name = name
        self.strand = strand
        self.exons = exons
        self.junctions = junctions
        self.chrom = chrom
        try:
            self.start = exons[0].get_coords()[0]
            self.end = exons[-1].get_coords()[1]
        except:
            self.start = None
            self.end = None

    def merge(self, other):
        # quick sanity check
        self_dict = self.__dict__.copy()
        other_dict = other.__dict__.copy()

        for attr in ['junctions', 'exons', 'end']:
            del self_dict[attr]
            del other_dict[attr]

        if self_dict != other_dict:
            voilaLog().warning('Attemping to merge two genes that aren\'t technically equal.', self_dict, other_dict)

        # adjust end point
        self.end = max(other.end, self.end)

        # concat exons and junctions
        self.exons += other.exons
        self.junctions += other.junctions

        # merge exons
        self.merge_overlapping_exons()

        # remove duplicate junctions
        self.remove_duplicate_junctions()

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
        return [self.start, self.end]

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def to_bed12(self):
        #  "%(a)s, %(a)s" % {'a':'test'}
        # fields = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
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
            ','.join([str(abs(e.get_coords()[0] - e.get_coords()[1])) for e in self.exons]),
            ','.join([str(abs(e.get_coords()[0] - self.start)) for e in self.exons]),
        ]
        return "\t".join([str(b) for b in bed_fields])

    def to_hdf5(self, h, use_id=True):
        if use_id:
            h = h.create_group(self.id)

        super(GeneGraphic, self).to_hdf5(h, use_id)

    def cls_list(self):
        return {'exons':
                    {'class': ExonGraphic, 'args': (None, None, None, None)},
                'junctions':
                    {'class': JunctionGraphic, 'args': ((), None, None)}
                }

    def merge_overlapping_exons(self):
        exons = sorted(self.exons, key=lambda exon: exon.coords)
        merged_exons = [exons[0]]
        for exon in exons[1:]:
            if not merged_exons[-1].merge(exon):
                merged_exons.append(exon)

        self.exons = merged_exons

    def remove_duplicate_junctions(self):
        self.junctions = sorted(set(self.junctions), key=lambda junction: junction.coords)

    def get_missing_exons(self, master_gene):

        # Create splice graph for this gene as compared to master splice graph
        master_exons_iter = iter(master_gene.exons)
        sample_exons_iter = iter(self.exons)

        # Figure out which exons are missing from the sample and make sure that the sample has the updated exons
        m = master_exons_iter.next()
        s = sample_exons_iter.next()

        missing = True
        output_exons = []

        try:

            while True:
                if m.overlaps(s):
                    missing = False
                    s = sample_exons_iter.next()
                else:
                    m = m.copy()

                    if missing:
                        m.type_exon = constants.EXON_TYPE_DB

                    output_exons.append(m)
                    missing = True
                    m = master_exons_iter.next()

        except StopIteration:
            output_exons.append(m)
            try:
                while True:
                    output_exons.append(master_exons_iter.next())
            except StopIteration:
                pass

        self.exons = output_exons

    def get_missing_junctions(self, master_gene):

        for junction in master_gene.junctions:
            if not junction in self.junctions:
                jg = junction.copy()
                jg.type_junction = constants.JUNCTION_TYPE_DB
                self.junctions.append(jg)
        self.junctions.sort(key=lambda junction: junction.coords)

    @classmethod
    def create_master(cls, splice_graphs, gene_id):
        master = None
        for h5_file in splice_graphs:
            voilaLog().info('loading ' + h5_file)
            with h5py.File(h5_file) as h:
                gene = cls.easy_from_hdf5(h, gene_id)
                try:
                    master.merge(gene)
                except AttributeError:
                    master = gene
        return master

    @classmethod
    def easy_from_hdf5(cls, h, gene_id):
        return GeneGraphic(None).from_hdf5(h[gene_id])

    def __str__(self):
        return str(self.__dict__)


class ExonException(Exception):
    pass


class ExonGraphic(HDF5):
    def __init__(self, a3, a5, coords, type_exon, coords_extra=list(), intron_retention=None, lsv_type=0,
                 alt_starts=list(), alt_ends=list()):
        """
        ExonGraphic constructor

        :param a3: list of junction indexes for alternative 3'acceptor sites.
        :param a5: list of junction indexes for alternative 5' donor sites.
        :param coords: exon coords.
        :param type_exon: 0: Annotated+found in RNASeq data; 1: Only found in RNASeq data; 2: Only annotated.
        :param coords_extra: coordinates of parts of the exon detected de novo.
        :param intron_retention: boolean to indicate that the intron after the exon has been identified as retained.
        :param lsv_type: signature of the lsv (i.e. t|1e1.1|1e1.2). [Currently not used]
        :param alt_starts: list of coordinates for alternative start sites of the exon (no junction "landing" on it)
        :param alt_ends: list of coordinates for alternative end sites of the exon (no junction "departing" from it)
        """
        super(ExonGraphic, self).__init__()
        self.a3 = a3
        self.a5 = a5
        self.coords = coords
        self.type_exon = type_exon
        self.coords_extra = coords_extra

        self.intron_retention = intron_retention
        self.lsv_type = lsv_type

        self.alt_starts = alt_starts
        self.alt_ends = alt_ends

    def merge(self, other):

        if not self.overlaps(other):
            return False

        other_dict = other.__dict__.copy()

        self.coords = [min(self.coords[0], other_dict['coords'][0]), max(self.coords[1], other_dict['coords'][1])]
        del other_dict['coords']

        if other_dict['type_exon'] < 4 and self.type_exon < 4:
            self.type_exon = min(other_dict['type_exon'], self.type_exon)
        elif other_dict['type_exon'] != self.type_exon:
            voilaLog().warning('Attempting to merge a missing end exon with a normal exon. ' + str((self, other)))

        del other_dict['type_exon']

        for attr in other_dict:
            try:
                for item in other_dict[attr]:
                    if item not in self.__dict__[attr]:
                        self.__dict__[attr].append(item)
            except TypeError:
                # if object attribute isn't iterable...
                pass

        return True

    def overlaps(self, other):
        assert self.coords[0] <= other.coords[0]
        return other.coords[0] <= self.coords[1]

    def get_a3_list(self):
        return self.a3

    def get_a5_list(self):
        return self.a5

    def get_coords(self):
        # Mask unkonwn start or ends
        def mask_unknown(coords):
            if coords[0] is None:
                coords[0] = -1
            if coords[1] is None:
                coords[1] = -1
            return coords

        return mask_unknown(list(self.coords))

    def get_type(self):
        return self.type_exon

    def get_coords_extra(self):
        return self.coords_extra

    def get_intron_retention(self):
        return self.intron_retention

    def get_lsv_type(self):
        return self.lsv_type

    def set_intron_retention(self, intron_retention):
        self.intron_retention = intron_retention

    def get_alt_starts(self):
        return self.alt_starts

    def get_alt_ends(self):
        return self.alt_ends

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def __str__(self):
        return str((self.coords, self.type_exon))

    def __repr__(self):
        return self.__str__()

    def __copy__(self):
        eg = type(self)(None, None, None, None)
        eg.__dict__.update(self.__dict__)
        return eg

    def copy(self):
        return self.__copy__()


class JunctionGraphic(HDF5):
    def __init__(self, coords, type_junction, nreads, clean_nreads=0, transcripts=list(), ir=0):
        super(JunctionGraphic, self).__init__()
        self.coords = list(coords)
        self.type_junction = type_junction
        self.num_reads = nreads
        self.num_clean_reads = clean_nreads
        self.transcripts = transcripts
        self.ir = ir

    def get_coords(self):
        return self.coords

    def get_type(self):
        return self.type_junction

    def get_num_reads(self):
        return self.num_reads

    def get_num_clean_reads(self):
        return self.num_clean_reads

    def set_clean_reads(self, cr):
        self.num_clean_reads = cr

    def get_transcripts(self):
        return self.transcripts

    def set_transcripts(self, t):
        self.transcripts = t

    def get_ir(self):
        return self.ir

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def __hash__(self):
        return int(str(self.coords[0]) + str(self.coords[1]))

    def __eq__(self, other):
        return self.coords == other.coords

    def __str__(self):
        return str([self.coords, self.type_junction])

    def __repr__(self):
        return self.__str__()

    def __copy__(self):
        jg = type(self)((), None, None)
        jg.__dict__.update(self.__dict__)
        return jg

    def copy(self):
        return self.__copy__()


def splice_graph_from_hdf5(hdf5_filename):
    def worker():
        with h5py.File(hdf5_filename, 'r') as h:
            while True:
                id = queue.get()
                manager_dict[id] = GeneGraphic(None).from_hdf5(h[id])
                queue.task_done()

    def producer():
        with h5py.File(hdf5_filename, 'r') as h:
            for x in h:
                queue.put(x)

    log = voilaLog()

    if not os.path.isfile(hdf5_filename):
        log.error('unable to load file: {0}'.format(hdf5_filename))
        return

    log.info('Loading {0}.'.format(hdf5_filename))

    queue = JoinableQueue()
    manager_dict = Manager().dict()

    producer_proc = Process(target=producer)
    producer_proc.daemon = True
    producer_proc.start()

    pool = Pool(PROCESS_COUNT, worker)

    producer_proc.join()
    queue.join()

    pool.close()
    queue.close()

    return manager_dict.values()


class LsvGraphic(GeneGraphic):
    def __init__(self, type_lsv, coords, id, name=None, strand=None, exons=list(), junctions=list(), chrom=None):
        super(LsvGraphic, self).__init__(id, name, strand, exons, junctions, chrom)
        self.type = type_lsv
        self.coords = coords

    def get_type(self):
        return self.type

    def to_hdf5(self, h, use_id=True):
        super(LsvGraphic, self).to_hdf5(h, False)

    @classmethod
    def easy_from_hdf5(cls, h, lsv_id=None):
        return cls(None, None, None).from_hdf5(h)
