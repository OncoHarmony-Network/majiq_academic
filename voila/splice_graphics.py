import csv
import io
import json
import os
import sys
from multiprocessing import Manager
from multiprocessing import Process, Pool
from multiprocessing.queues import JoinableQueue
from os.path import join, isfile

import h5py

from voila import constants
from voila.constants import PROCESS_COUNT, EXPERIMENTS_NAMES
from voila.hdf5 import HDF5, ExonTypeDataSet, JunctionTypeDataSet, ReadsDataSet
from voila.utils.voilaLog import voilaLog


class ExperimentIndexError(IndexError):
    def __init__(self):
        super(ExperimentIndexError, self).__init__(
            'Attempted to access an out of range experiment.')


class UnknownField(Exception):
    def __init__(self, field):
        super(UnknownField, self).__init__(
            '"{0}" might not contain data that has been sorted into a list by experiment.'.format(field))


class Experiment(HDF5):
    def field_by_experiment(self):
        return ()

    def get(self, field, experiment):
        if field in self.field_by_experiment():
            try:
                return self.__dict__[field][experiment]
            except IndexError:
                raise ExperimentIndexError()
        else:
            raise UnknownField(field)

    def get_dict(self, experiment):
        d = {}
        for key in self.__dict__:
            if key in self.field_by_experiment():
                if self.__dict__[key]:
                    try:
                        d[key] = self.__dict__[key][experiment]
                    except IndexError:
                        raise ExperimentIndexError()
            else:
                d[key] = self.__dict__[key]
        return d


class NoExons(Exception):
    def __init__(self):
        super(NoExons, self).__init__('There no exons in this gene.')


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

    def start(self):
        """
        Start of gene.
        :return: Integer
        """
        try:
            return self.exons[0].start
        except IndexError:
            raise NoExons()

    def end(self):
        """
        End of gene.
        :return: Integer
        """
        try:
            return self.exons[-1].end
        except IndexError:
            raise NoExons()

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
            voilaLog().warning('Attemping to merge two genes that might not be equal.')

        # concat exons and junctions
        self.exons += other.exons
        self.junctions += other.junctions

        # merge exons
        self.merge_overlapping_exons()

        # remove duplicate junctions
        self.remove_duplicate_junctions()

    def get_coordinates(self):
        """
        Construct and return this gene's coordinates.
        :return: list of coordinates
        """
        return [self.start, self.end]

    def to_json(self, encoder=json.JSONEncoder):
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
            h = h.create_group(join('/genes', self.gene_id))

        super(GeneGraphic, self).to_hdf5(h, use_id)

    def cls_list(self):
        return {'exons': ExonGraphic, 'junctions': JunctionGraphic}

    def merge_overlapping_exons(self):
        """
        For all the exons for this gene, merge any that overlap.
        :return: None
        """
        exons = sorted(self.exons)
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

                    if missing and m.type_exon <= constants.EXON_TYPE_DB_OTHER_RNASEQ:
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
            if junction not in self.junctions:
                jg = junction.copy()
                jg.type_junction = constants.JUNCTION_TYPE_DB
                self.junctions.append(jg)
        self.junctions.sort()

    def get_dict(self, experiment):
        d = self.__dict__.copy()
        d['exons'] = [e.get_dict(experiment) for e in self.exons]
        d['junctions'] = [j.get_dict(experiment) for j in self.junctions]
        d['start'] = self.start()
        d['end'] = self.end()
        d['length'] = len(self)
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

    def field_by_experiment(self):
        return ['exon_type']

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

    def to_json(self, encoder=json.JSONEncoder):
        """
        JSON representation of this class.
        :param encoder: JSON encoder passed to json.dumps
        :return: JSON string
        """
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, cls=encoder)

    def exclude(self):
        return ['exon_type']

    def to_hdf5(self, h, use_id=True):
        super(ExonGraphic, self).to_hdf5(h, use_id)

        ExonTypeDataSet(h, self.exon_type).encode()

    def from_hdf5(self, h):
        self.exon_type = ExonTypeDataSet(h).decode()
        return super(ExonGraphic, self).from_hdf5(h)

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

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None, None, None, None).from_hdf5(h)

    def __str__(self):
        return str(((self.start, self.end), self.exon_type))

    def __repr__(self):
        return self.__str__()

    def __sizeof__(self):
        return sum([sys.getsizeof(self.__dict__[key] for key in self.__dict__)])

    def __lt__(self, other):
        return self.start < other.start

    def __contains__(self, item):
        return self.start <= item <= self.end


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
        assert len(junction_type_list) == len(reads_list), 'There is an uneven number of experiments.'

        super(JunctionGraphic, self).__init__()
        self.start = start
        self.end = end
        self.junction_type = junction_type_list
        self.reads = reads_list
        self.clean_reads = None
        self.transcripts = transcripts
        self.intron_retention = intron_retention

    def field_by_experiment(self):
        return ['junction_type', 'reads', 'clean_reads']

    def to_json(self, encoder=json.JSONEncoder):
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

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None, (), ()).from_hdf5(h)

    def __hash__(self):
        return int(str(self.start) + str(self.end))

    def __eq__(self, other):
        return (self.start, self.end) == (other.start, other.end)

    def __copy__(self):
        jg = type(self)((), None, None)
        jg.__dict__.update(self.__dict__)
        return jg

    def __str__(self):
        return str(((self.start, self.end), self.junction_type))

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        return (self.start, self.end) < (other.start, other.end)


class LsvGraphic(GeneGraphic):
    def __init__(self, lsv_type, start, end, gene_id, name=None, strand=None, exons=(), junctions=(),
                 chromosome=None):
        """
        LSV Data.
        :param lsv_type: LSV type.  See constants file.
        :param gene_id: LSV id
        :param name: LSV name
        :param strand: LSV Strand. Either '-' or '+'
        :param exons: List of exons associated with this LSV
        :param junctions: List of junctions associtated with this LSV
        :param chromosome: This LSV's Chromosome
        """
        super(LsvGraphic, self).__init__(gene_id, name, strand, exons, junctions, chromosome)
        self.end = end
        self.start = start
        self.lsv_type = lsv_type

    def start(self):
        """
        LSV has defined start, which is different then its super class where it's generated. This method is here to
        protect us from accidentally calling the super class 'start()' method.
        :return: integer
        """
        return self.start

    def end(self):
        """
        LSV has a defined end, which is different then its super where it's generated. This method is here to protect
        us from accidentally calling the super class 'end()' method.
        :return:
        """
        return self.end

    def to_hdf5(self, h, use_id=True):
        super(LsvGraphic, self).to_hdf5(h, False)

    @classmethod
    def easy_from_hdf5(cls, h):
        return cls(None, None, None, None).from_hdf5(h)


class Splicegraph(object):
    def __init__(self, splice_graph_file_name, mode):
        self.file_name = splice_graph_file_name
        self.hdf5 = None
        self.mode = mode

    def __enter__(self):
        self.hdf5 = h5py.File(self.file_name, self.mode)
        return self

    def __exit__(self, type, value, traceback):
        self.hdf5.close()

    def add_experiments_names(self, experiments_names):
        self.hdf5[EXPERIMENTS_NAMES] = experiments_names

    def erase_splice_graph_file(self):
        os.remove(self.file_name)
        self.__enter__()

    def add_gene(self, gene):
        gene.to_hdf5(self.hdf5)

    def get_genes_list(self):
        def worker():
            with h5py.File(self.file_name, 'r') as h:
                while True:
                    id = queue.get()
                    manager_dict[id] = GeneGraphic.easy_from_hdf5(h['/genes'][id])
                    queue.task_done()

        def producer():
            with h5py.File(self.file_name, 'r') as h:
                for x in h['/genes']:
                    queue.put(x)

        log = voilaLog()

        if not isfile(self.file_name):
            log.error('unable to load file: {0}'.format(self.file_name))
            return

        log.info('Loading {0}.'.format(self.file_name))

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

    def get_experiments_list(self):
        return self.hdf5[EXPERIMENTS_NAMES].value

    def get_gene_ids(self):
        return self.hdf5['genes'].keys()

    def get_gene(self, gene_id):
        return GeneGraphic.easy_from_hdf5(self.hdf5['genes'][gene_id])
