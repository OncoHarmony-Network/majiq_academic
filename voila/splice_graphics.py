import csv
import io
import json

import h5py

from voila import constants
from voila.hdf5 import HDF5
from voila.utils.voilaLog import voilaLog


class GeneGraphic(HDF5):
    def __init__(self, id, name=None, strand=None, exons=list(), junctions=list(), chrom=None):
        """
        Gene data.
        :param id: Gene ID
        :param name: Gene name
        :param strand: Gene strand, either '-' or '+'
        :param exons: List of this gene's exons
        :param junctions: List of this gene's junctions
        :param chrom: String name of this gene's chromosome
        """

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
        """
        Merge another gene into this one.
        :param other: Other gene object
        :return: None
        """

        # for a sanity check, we'll remove the known changing elements and verify the rest are equal
        self_dict = self.__dict__.copy()
        other_dict = other.__dict__.copy()
        for attr in ['junctions', 'exons', 'end']:
            del self_dict[attr]
            del other_dict[attr]

        # log warning if gene's some seem equal
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

    def get_coords(self):
        """
        Construct and return this gene's coordinates.
        :return: list of coordinates
        """
        return [self.start, self.end]

    def to_JSON(self, encoder=json.JSONEncoder):
        """
        Generate and return JSON representation for a gene object.
        :param encoder: encoder used in the json.dumps call
        :return: JSON string
        """
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def to_bed12(self):
        """
        Bed 12 representation of this gene object.
        :return: bed 12 string
        """
        fieldnames = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
                      'blockCount', 'blockSizes', 'blockStarts']
        csvfile = io.BytesIO()
        writer = csv.DictWriter(csvfile, fieldnames, delimiter='\t')
        writer.writerow({
            'chrom': self.chrom,
            'chromStart': self.start,
            'chromEnd': self.end,
            'name': 'lsv',
            'score': 0,
            'strand': self.strand,
            'thickStart': self.start,
            'thickEnd': self.end,
            'itemRgb': 0,
            'blockCount': len(self.exons),
            'blockSizes': ','.join([str(abs(e.coords[0] - e.coords[1])) for e in self.exons]),
            'blockStarts': ','.join([str(abs(e.coords[0] - self.start)) for e in self.exons])
        })
        return csvfile.getvalue()

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
        """
        For all the exons for this gene, merge any that overlap.
        :return: None
        """
        exons = sorted(self.exons, key=lambda exon: exon.coords)
        merged_exons = [exons[0]]
        for exon in exons[1:]:
            if not merged_exons[-1].merge(exon):
                merged_exons.append(exon)

        self.exons = merged_exons

    def remove_duplicate_junctions(self):
        """
        Remove any duplicate junctions for this gene.
        :return: None
        """
        self.junctions = sorted(set(self.junctions), key=lambda junction: junction.coords)

    def get_missing_exons(self, master_gene):
        """
        Get any exons that this gene is missing from a "master gene".
        :param master_gene: A gene object containing a "master" list of exons and junctions
        :return: None
        """
        # Create splice graph for this gene as compared to master splice graph
        master_exons_iter = iter(master_gene.exons)
        sample_exons_iter = iter(self.exons)

        # Figure out which exons are missing from the sample and make sure that the sample has the updated exons
        m = master_exons_iter.next()
        s = sample_exons_iter.next()

        missing = True
        exons = []

        try:

            while True:
                if m.overlaps(s):
                    missing = False
                    s = sample_exons_iter.next()
                else:
                    m = m.copy()

                    if missing:
                        m.type_exon = constants.EXON_TYPE_DB

                    exons.append(m)
                    missing = True
                    m = master_exons_iter.next()

        except StopIteration:
            exons.append(m)
            try:
                while True:
                    exons.append(master_exons_iter.next())
            except StopIteration:
                pass

        self.exons = exons

    def get_missing_junctions(self, master_gene):
        """
        Get any missing junctions from a "master gene".
        :param master_gene: A gene object containing a "master" list of exons and junctions
        :return: None
        """
        for junction in master_gene.junctions:
            if not junction in self.junctions:
                jg = junction.copy()
                jg.type_junction = constants.JUNCTION_TYPE_DB
                self.junctions.append(jg)
        self.junctions.sort(key=lambda junction: junction.coords)


    @classmethod
    def create_master(cls, splice_graphs, gene_id):
        """
        Creates a "master gene" from a list of splice graph files and the gene id.
        :param splice_graphs: List of splice graph files
        :param gene_id: The gene id
        :return: Gene object
        """
        master = None
        for h5_file in splice_graphs:
            voilaLog().info('loading ' + h5_file)
            with h5py.File(h5_file) as h:
                gene = cls.easy_from_hdf5(h[gene_id])
                try:
                    master.merge(gene)
                except AttributeError:
                    master = gene
        return master

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None).from_hdf5(h)


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
        """
        Merge another exon with this one.
        :param other: Exon object
        :return: Boolean. True if exon is merged with this one.
        """

        # if the other exon doesn't overlap, then return False
        if not self.overlaps(other):
            return False

        other_dict = other.__dict__.copy()

        self.coords = [min(self.coords[0], other_dict['coords'][0]), max(self.coords[1], other_dict['coords'][1])]
        del other_dict['coords']

        self.type_exon = min(other_dict['type_exon'], self.type_exon)
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
        """
        Checks if another exon overlaps with this one.
        :param other: Exon object
        :return: Boolean.  True if other exon overlaps with this one.
        """
        assert self.coords[0] <= other.coords[0], 'Exons have to be in order to check if they overlap.'
        return other.coords[0] <= self.coords[1]

    def to_JSON(self, encoder=json.JSONEncoder):
        """
        JSON representation of this class.
        :param encoder: JSON encoder passed to json.dumps
        :return: JSON string
        """
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def __copy__(self):
        eg = type(self)(None, None, None, None)
        eg.__dict__.update(self.__dict__)
        return eg

    def copy(self):
        """
        Attribute for attribute copy of this class.
        :return: Copied class object
        """
        return self.__copy__()


class JunctionGraphic(HDF5):
    def __init__(self, coords, type_junction, nreads, clean_nreads=0, transcripts=list(), ir=0):
        """
        Junction data.
        :param coords: junction start and end points
        :param type_junction: junction type.  See constants file.
        :param nreads: number of reads
        :param clean_nreads: number of clean reads
        :param transcripts: list of transtripts
        :param ir: intron retention
        """
        super(JunctionGraphic, self).__init__()
        self.coords = list(coords)
        self.type_junction = type_junction
        self.num_reads = nreads
        self.num_clean_reads = clean_nreads
        self.transcripts = transcripts
        self.ir = ir

    def to_JSON(self, encoder=json.JSONEncoder):
        """
        JSON representation of this class.
        :param encoder: json encoder passed to json.dumps
        :return: JSON string
        """
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def copy(self):
        """
        Copy of this class
        :return: Copied object
        """
        return self.__copy__()

    def __hash__(self):
        return int(str(self.coords[0]) + str(self.coords[1]))

    def __eq__(self, other):
        return self.coords == other.coords

    def __copy__(self):
        jg = type(self)((), None, None)
        jg.__dict__.update(self.__dict__)
        return jg


class LsvGraphic(GeneGraphic):
    def __init__(self, type_lsv, coords, id, name=None, strand=None, exons=list(), junctions=list(), chrom=None):
        """
        LSV Data.
        :param type_lsv: LSV type.  See constants file.
        :param coords: Start and end of LSV
        :param id: LSV id
        :param name: LSV name
        :param strand: LSV Strand. Either '-' or '+'
        :param exons: List of exons associated with this LSV
        :param junctions: List of junctions associtated with this LSV
        :param chrom: This LSV's Chromosome
        """
        super(LsvGraphic, self).__init__(id, name, strand, exons, junctions, chrom)
        self.type = type_lsv
        self.coords = coords

    def to_hdf5(self, h, use_id=True):
        super(LsvGraphic, self).to_hdf5(h, False)

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None, None).from_hdf5(h)
