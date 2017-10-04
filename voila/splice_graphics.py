import csv
import io
import sys
from os.path import join

import h5py
import numpy

from voila import constants
from voila.hdf5 import HDF5, ExonTypeDataSet, JunctionTypeDataSet, ReadsDataSet
from voila.utils.exceptions import NoExonsInGene, ExperimentIndexError, ExperimentUnknowField


class Experiment(HDF5):
    def field_by_experiment(self):
        """
        List of fields that are actually a list of data for a set of experiments.
        :return: list
        """
        return ()

    def get(self, field, experiment):
        """
        Get data from a field for a specifice experiment index.
        :param field: class attribute name
        :param experiment: experiment index
        :return: dict
        """
        if field in self.field_by_experiment():
            experiment_field = self.__dict__[field]

            if experiment_field is None:
                raise Exception(experiment_field, field, self.__class__.__name__)

            try:
                value = experiment_field[experiment]
                if isinstance(value, (numpy.int64, numpy.uint64)):
                    value = int(value)
                elif isinstance(value, numpy.float64):
                    value = float(value)
                return value
            except IndexError:
                raise ExperimentIndexError()

        else:
            raise ExperimentUnknowField(field)

    def get_experiment(self, experiment):
        """
        Get data for a sp
        :param experiment:
        :return:
        """
        d = self.__dict__.copy()
        for field in self.__dict__:
            if field in self.field_by_experiment():
                d[field] = self.get(field, experiment)
        return d


class GeneGraphic(HDF5):
    def __init__(self, gene_id, name=None, strand=None, exons=(), junctions=(), chromosome=None):
        """
        Gene data.
        :param gene_id: Gene ID
        :param name: Gene name
        :param strand: Gene strand, either '-' or '+'
        :param exons: List of this gene's exons
        :param junctions: List of this gene's junctions
        :param chromosome: String name of this gene's chromosome
        """

        # sanity check
        exon_lengths = [len(exon.exon_type) for exon in exons]
        junction_lengths = [len(junction.junction_type) for junction in junctions]
        assert all(item == exon_lengths[0] for item in
                   exon_lengths + junction_lengths), 'There is an uneven number of experiments.'

        super(GeneGraphic, self).__init__()
        self.gene_id = gene_id
        self.name = name
        self.strand = strand
        self.exons = sorted(exons)
        self.junctions = sorted(junctions)
        self.chromosome = chromosome

    def __len__(self):
        return self.end() - self.start()

    def __sizeof__(self):
        return sum([sys.getsizeof(self.__dict__[key]) for key in self.__dict__])

    def __eq__(self, other):
        return (self.chromosome == other.chrom and self.strand == other.strand
                and self.start < other.end and self.end > other.start)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return (self.chromosome < other.chromosome
                or (self.chromosome == other.chromosome
                    and (self.end < other.start
                         or (self.end > other.start and self.start < other.end
                             and self.strand == '+' and other.strand == '-'))))

    def __gt__(self, other):
        return (self.chromosome > other.chrom
                or (self.chromosome == other.chrom
                    and (self.start > other.end
                         or (self.end > other.start and self.start < other.end
                             and self.strand == '-' and other.strand == '+'))))

    def __copy__(self):
        gg = type(self)(None, None, (), ())
        gg.__dict__.update(self.__dict__)
        return gg

    def copy(self):
        """
        Attribute for attribute copy of this class.
        :return: Copied class object
        """
        return self.__copy__()

    def start(self):
        """
        Start of gene.
        :return: Integer
        """
        try:
            return self.exons[0].start
        except IndexError:
            raise NoExonsInGene(self.gene_id)

    def end(self):
        """
        End of gene.
        :return: Integer
        """
        try:
            return self.exons[-1].end
        except IndexError:
            raise NoExonsInGene(self.gene_id)

    def combine(self, experiment_index, gene_dict=None):
        """
        Used to create the "Combined" option in a splice graph when there's more than one experiment.
        :param experiment_index:
        :param gene_dict:
        :return:
        """
        if not gene_dict:
            return self.get_experiment(experiment_index)

        assert self.gene_id == gene_dict['gene_id'], 'cannot combine two different genes'

        gene_dict['junctions'] = [s.combine(experiment_index, d) for s, d in
                                  zip(self.junctions, gene_dict['junctions'])]

        gene_dict['exons'] = [s.combine(experiment_index, d) for s, d in
                              zip(self.exons, gene_dict['exons'])]

        return gene_dict

    def get_coordinates(self):
        """
        Construct and return this gene's coordinates.
        :return: list of coordinates
        """
        return self.start(), self.end()

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
            'chrom': self.chromosome,
            'chromStart': self.start,
            'chromEnd': self.end,
            'name': 'lsv',
            'score': 0,
            'strand': self.strand,
            'thickStart': self.start,
            'thickEnd': self.end,
            'itemRgb': 0,
            'blockCount': len(self.exons),
            'blockSizes': ','.join([str(abs(e.start - e.end)) for e in self.exons]),
            'blockStarts': ','.join([str(abs(e.start - self.start())) for e in self.exons])
        })
        return csvfile.getvalue()

    def to_hdf5(self, h, use_id=True):
        if use_id:
            try:
                h = h.create_group(join('/genes', self.gene_id))
            except ValueError:
                h = h[join('/genes', self.gene_id)]

        super(GeneGraphic, self).to_hdf5(h, use_id)

    def cls_dict(self):
        return {'exons': ExonGraphic, 'junctions': JunctionGraphic}

    def get_experiment(self, experiment):
        d = self.__dict__.copy()
        experiment = int(experiment)
        d.update({
            'exons': [e.get_experiment(experiment) for e in self.exons],
            'junctions': [j.get_experiment(experiment) for j in self.junctions],
            'start': self.start(),
            'end': self.end(),
            'length': len(self)
        })
        return d

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
            with h5py.File(h5_file) as h:
                try:
                    gene = cls.easy_from_hdf5(h[gene_id])
                    try:
                        master.merge(gene)
                    except AttributeError:
                        master = gene
                except KeyError:
                    pass
        return master

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None).from_hdf5(h)


class ExonGraphic(Experiment):
    def __init__(self, a3, a5, start, end, exon_type_list, coords_extra=list(), intron_retention=None, lsv_type=0,
                 alt_starts=list(), alt_ends=list()):
        """
        ExonGraphic constructor

        :param start: Start coordinate for exon
        :param end: End coordinate for exon
        :param a3: list of junction indexes for alternative 3'acceptor sites.
        :param a5: list of junction indexes for alternative 5' donor sites.
        :param exon_type_list: 0: Annotated+found in RNASeq data; 1: Only found in RNASeq data; 2: Only annotated.
        :param coords_extra: coordinates of parts of the exon detected de novo.
        :param intron_retention: boolean to indicate that the intron after the exon has been identified as retained.
        :param lsv_type: signature of the lsv (i.e. t|1e1.1|1e1.2). [Currently not used]
        :param alt_starts: list of coordinates for alternative start sites of the exon (no junction "landing" on it)
        :param alt_ends: list of coordinates for alternative end sites of the exon (no junction "departing" from it)
        """
        super(ExonGraphic, self).__init__()
        self.end = end
        self.start = start
        self.a3 = a3
        self.a5 = a5
        self.exon_type = exon_type_list
        self.coords_extra = coords_extra

        self.intron_retention = intron_retention
        self.lsv_type = lsv_type

        self.alt_starts = alt_starts
        self.alt_ends = alt_ends

    def __sizeof__(self):
        return sum([sys.getsizeof(self.__dict__[key] for key in self.__dict__)])

    def __lt__(self, other):
        return self.start < other.start

    def __contains__(self, item):
        return self.start <= item <= self.end

    def __len__(self):
        return self.end - self.start

    def coords(self):
        return self.start, self.end

    def field_by_experiment(self):
        return ['exon_type']

    def combine(self, experiment_index, exon_dict):
        exon_dict['exon_type'] = min(self.exon_type[experiment_index], exon_dict['exon_type'])
        return exon_dict

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

        self.start = min(self.start, other_dict['start'])
        self.end = max(self.end, other_dict['end'])
        del other_dict['start']
        del other_dict['end']

        # if both exons are 'normal'...
        if self.exon_type <= 3 and other_dict['type_exon'] <= 3:
            self.exon_type = min(other_dict['type_exon'], self.exon_type)
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
        assert self.start <= other.start, 'Exons have to be in order to check if they overlap.'
        return other.start <= self.end

    def exclude(self):
        return ['exon_type']

    def to_hdf5(self, h, use_id=True):
        super(ExonGraphic, self).to_hdf5(h, use_id)

        ExonTypeDataSet(h, self.exon_type).encode()

    def from_hdf5(self, h):
        self.exon_type = ExonTypeDataSet(h).decode()
        return super(ExonGraphic, self).from_hdf5(h)

    def __copy__(self):
        eg = type(self)(None, None, None, None, None)
        eg.__dict__.update(self.__dict__)
        return eg

    def copy(self):
        """
        Attribute for attribute copy of this class.
        :return: Copied class object
        """
        return self.__copy__()

    def to_dict(self):
        return self.__dict__.copy()

    def get_experiment(self, experiment):
        d = super(ExonGraphic, self).get_experiment(experiment)
        d['length'] = len(self)
        return d

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None, None, None, None).from_hdf5(h)


class JunctionGraphic(Experiment):
    def __init__(self, start, end, junction_type_list, reads_list, transcripts=(),
                 intron_retention=constants.NONE_IR_TYPE):
        """
        Junction data.
        :param start: Start coordinate for junction
        :param end: End coordinate for junction
        :param junction_type_list: junction type.  See constants file.
        :param reads_list: number of reads
        :param transcripts: list of transtripts
        :param intron_retention: intron retention
        """

        # sanity check
        assert reads_list is None or (len(junction_type_list) == len(reads_list)), 'There is an uneven number of experiments.'

        super(JunctionGraphic, self).__init__()
        self.start = start
        self.end = end
        self.junction_type = junction_type_list
        self.reads = reads_list
        self.clean_reads = None
        self.transcripts = transcripts
        self.intron_retention = intron_retention

    def __hash__(self):
        return int(str(self.start) + str(self.end))

    def __eq__(self, other):
        return (self.start, self.end) == (other.start, other.end)

    def __copy__(self):
        jg = type(self)(None, None, (), ())
        jg.__dict__.update(self.__dict__)
        return jg

    def __str__(self):
        return str(((self.start, self.end), self.junction_type))

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return (self.start, self.end) < (other.start, other.end)

    def coords(self):
        return self.start, self.end

    def field_by_experiment(self):
        return ['junction_type', 'reads']

    def copy(self):
        """
        Copy of this class
        :return: Copied object
        """
        return self.__copy__()

    def combine(self, experiment_index, junc_dict):
        junc_dict['reads'] += self.reads[experiment_index]
        junc_dict['junction_type'] = min(self.junction_type[experiment_index], junc_dict['junction_type'])
        if 'clean_reads' in junc_dict and junc_dict['clean_reads']:
            junc_dict['clean_reads'] += self.clean_reads[experiment_index]
        return junc_dict

    def exclude(self):
        return ['junction_type', 'reads']

    def to_hdf5(self, h, use_id=True):
        super(JunctionGraphic, self).to_hdf5(h, use_id)
        JunctionTypeDataSet(h, self.junction_type).encode()
        ReadsDataSet(h, self.reads).encode()

    def from_hdf5(self, h):
        self.junction_type = JunctionTypeDataSet(h).decode()
        self.reads = ReadsDataSet(h).decode()
        return super(JunctionGraphic, self).from_hdf5(h)

    def to_dict(self):
        return self.__dict__.copy()

    def junction_id(self, chromosome, strand):
        return '{0}:{1}:{2}-{3}'.format(chromosome, strand, self.start, self.end)

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None, (), ()).from_hdf5(h)


class LsvGraphic(HDF5):
    def __init__(self, lsv_type, start, end, lsv_id, name=None, strand=None, exons=(), junctions=(), chromosome=None):
        """
        LSV Data.
        :param lsv_type: LSV type.  See constants file.
        :param lsv_id: LSV id
        :param name: LSV name
        :param strand: LSV Strand. Either '-' or '+'
        :param exons: List of exons associated with this LSV
        :param junctions: List of junctions associtated with this LSV
        :param chromosome: This LSV's Chromosome
        """
        super(LsvGraphic, self).__init__()
        self.end = end
        self.start = start
        self.lsv_type = lsv_type
        self.lsv_id = lsv_id
        self.name = name
        self.strand = strand
        self.exons = exons
        self.junctions = junctions
        self.chromosome = chromosome

    def to_dict(self, ignore_keys=()):
        d = self.__dict__.copy()

        for key in ignore_keys:
            del d[key]

        d.update({
            'exons': [exon.to_dict() for exon in self.exons],
            'junctions': [junction.to_dict() for junction in self.junctions]
        })
        return d

    def cls_dict(self):
        return {'exons': ExonGraphic, 'junctions': JunctionGraphic}

    def junction_ids(self):
        return (junction.junction_id(self.chromosome, self.strand) for junction in self.junctions)

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None, None, None).from_hdf5(h)
