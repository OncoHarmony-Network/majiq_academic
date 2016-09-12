import json
from multiprocessing import Manager, Pool
from multiprocessing.process import Process
from multiprocessing.queues import JoinableQueue

import h5py

from voila.hdf5 import HDF5


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
        return [self.exons[0].get_coords()[0], self.exons[-1].get_coords()[1]]
        # if self.strand == '-':
        #     return [self.exons[-1].coords[1], self.exons[0].coords[0]]

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
            # [::int(self.strand + '1')]
            ','.join([str(abs(e.get_coords()[0] - self.start)) for e in self.exons]),
        ]
        return "\t".join([str(b) for b in bed_fields])

    def cls_list(self):
        return {'exons':
                    {'class': ExonGraphic, 'args': (None, None, None, None)},
                'junctions':
                    {'class': JunctionGraphic, 'args': ((), None, None)}
                }


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

    # def get_size(self):
    #     return self.

    def to_JSON(self, encoder=json.JSONEncoder):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)


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


class LsvGraphic(GeneGraphic):
    def __init__(self, type_lsv, coords, id, name=None, strand=None, exons=list(), junctions=list(), chrom=None):
        super(LsvGraphic, self).__init__(id, name, strand, exons, junctions, chrom)
        self.type = type_lsv
        self.coords = coords

    def get_type(self):
        return self.type


def splice_graph_from_hdf5(hdf5_filename, logger):
    def worker():
        with h5py.File(hdf5_filename, 'r', swmr=True) as h:
            while True:
                index = queue.get()
                manager_dict[index] = GeneGraphic(None).from_hdf5(h[index])
                queue.task_done()

    def producer():
        with h5py.File(hdf5_filename, 'r', swmr=True) as h:
            for x in h:
                queue.put(x)

    if logger:
        logger.info('Loading {0}.'.format(hdf5_filename))

    queue = JoinableQueue()
    manager_dict = Manager().dict()

    producer_proc = Process(target=producer)
    producer_proc.daemon = True
    producer_proc.start()

    pool = Pool(None, worker)

    producer_proc.join()
    queue.join()

    pool.close()
    queue.close()

    ggs = [None] * len(manager_dict)
    for key in manager_dict.keys():
        ggs[int(key)] = manager_dict[key]

    return ggs
