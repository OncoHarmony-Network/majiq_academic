import math
import os
from abc import ABC, abstractmethod
from itertools import zip_longest

from sqlalchemy import create_engine, event, exists, and_
from sqlalchemy.orm import sessionmaker

from voila import constants
from voila.api import splice_graph_model as model
from voila.utils.voila_log import voila_log

Session = sessionmaker()


class SpliceGraphSQL:
    def __init__(self, filename, delete=False):
        self.filename = filename
        if delete is True:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

        engine = create_engine('sqlite:///{0}'.format(filename))
        engine.execution_options(stream_results=True)
        event.listen(engine, 'connect', self._fk_pragma_on_connect)
        model.Base.metadata.create_all(engine)
        Session.configure(bind=engine)
        self.session = Session()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @staticmethod
    def _fk_pragma_on_connect(dbapi_con, con_record):
        dbapi_con.execute('pragma foreign_keys=ON')

    def close(self):
        self.commit()
        self.session.close_all()

    def commit(self):
        self.session.commit()

    def add_experiment_names(self, experiment_names):
        session = self.session
        session.add_all([model.Experiment(name=name) for name in experiment_names])

    def get_experiment_names(self):
        return (e for e, in self.session.query(model.Experiment.name).all())

    def check_version(self):
        voila_log().warning('need to implement check version.')

    def get_experiments(self):
        return tuple(e[0] for e in self.session.query(model.Experiment.name).all())

    ######################################################
    # These are Voila specific methods... they'll need to be refactored out of this class.
    ######################################################


    def get_page_count(self, args):

        log = voila_log()

        log.debug('Start page count')

        gene_count = self.session.query(model.Gene).count()
        if args.limit < gene_count:
            gene_count = args.limit

        # if hasattr(args, 'voila_file'):
        #     with VoilaHDF5(args.voila_file, 'r') as v:
        #         for gene_id in self.get_gene_ids(args):
        #             try:
        #                 if any(v.get_lsvs(args, gene_id)):
        #                     gene_count += 1
        #             except GeneIdNotFoundInVoilaFile:
        #                 pass
        #
        # else:
        #     log.debug('Gene limit is set to {0}'.format(args.limit))
        #     for _ in self.get_gene_ids(args):
        #         gene_count += 1
        #         if gene_count == args.limit:
        #             break

        return int(math.ceil(gene_count / float(constants.MAX_GENES)))

    def get_paginated_genes(self, args):
        def grouper(iterable, n, fillvalue=None):
            "Collect data into fixed-length chunks or blocks"
            # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
            args = [iter(iterable)] * n
            return zip_longest(*args, fillvalue=fillvalue)

        for page in grouper(self.get_gene_ids(args), constants.MAX_GENES):
            yield page


            # log = voila_log()
            # log.debug('Getting paginated genes')
            #
            # gene_list = []
            # gene_count = 0
            #
            # for gene_id in self.get_gene_ids(args):
            #     log.debug('Found {0}'.format(gene_id))
            #     gene_list.append(gene_id)
            #     gene_count += 1
            #
            #     if gene_count == args.limit:
            #         break
            #
            #     if len(gene_list) == constants.MAX_GENES:
            #         yield gene_list
            #         gene_list = []
            #
            # if gene_list:
            #     yield gene_list

    def get_gene_ids(self, args=None):
        if args and args.gene_ids:
            return args.gene_ids

        if args and hasattr(args, 'lsv_ids') and args.lsv_ids:
            return (lsv_id.split(':')[0] for lsv_id in args.lsv_ids)

        for gene_id, in self.session.query(model.Gene.id).all():
            yield gene_id


class SpliceGraphType(ABC):
    def __bool__(self):
        return self.exists

    def __iter__(self):
        return self.get.__iter__()

    @abstractmethod
    def add(self, *args, **kwargs):
        pass

    @property
    @abstractmethod
    def get(self, *args):
        pass

    @property
    @abstractmethod
    def exists(self):
        pass


class Exons(SpliceGraphSQL):
    class _Exon(SpliceGraphType):
        def __init__(self, session, gene_id, start, end):
            self.session = session
            self.gene_id = gene_id
            self.start = start
            self.end = end

        def add(self, **kwargs):
            coords_extra = kwargs.pop('coords_extra', [])
            alt_ends = kwargs.pop('alt_ends', [])
            alt_starts = kwargs.pop('alt_starts', [])

            exon = model.Exon(gene_id=self.gene_id, start=self.start, end=self.end, **kwargs)

            for ce_start, ce_end in coords_extra:
                exon.coords_extra.append(model.CoordsExtra(start=int(ce_start), end=int(ce_end)))
            for alt_end in alt_ends:
                exon.alt_ends.append(model.AltEnds(coordinate=int(alt_end)))
            for alt_start in alt_starts:
                exon.alt_starts.append(model.AltStarts(coordinate=int(alt_start)))

            self.session.add(exon)

            return exon

        @property
        def get(self):
            return self.session.query(model.Exon).get((self.gene_id, self.start, self.end))

        @property
        def exists(self):
            return self.session.query(
                exists().where(and_(model.Exon.gene_id == self.gene_id, model.Exon.start == int(self.start),
                                    model.Exon.end == int(self.end)))).scalar()

    def exon(self, gene_id, start, end):
        return self._Exon(self.session, gene_id, start, end)

    @property
    def exons(self):
        return self.session.query(model.Exon).all()


class Junctions(SpliceGraphSQL):
    class _Junction(SpliceGraphType):
        def __init__(self, session, gene_id, start, end):
            self.session = session
            self.gene_id = gene_id
            self.start = int(start)
            self.end = int(end)

        def add(self, **kwargs):
            reads = kwargs.pop('reads', [])

            junc = model.Junction(gene_id=self.gene_id, start=self.start, end=self.end, **kwargs)

            for r, e in reads:
                junc.reads.append(model.Reads(reads=int(r), experiment_id=e))

            self.session.add(junc)

            return junc

        @property
        def exists(self):
            return self.session.query(
                exists().where(and_(model.Junction.gene_id == self.gene_id, model.Junction.start == self.start,
                                    model.Junction.end == self.end))).scalar()

        @property
        def get(self):
            return self.session.query(model.Junction).get((self.gene_id, self.start, self.end))

        def update_reads(self, reads, experiment):
            r = model.Reads(junction_gene_id=self.gene_id, junction_start=self.start, junction_end=self.end,
                            experiment_name=experiment, reads=int(reads))
            self.session.add(r)

    def junction(self, gene_id, start, end):
        return self._Junction(self.session, gene_id, start, end)

    @property
    def junctions(self):
        return self.session.query(model.Junction).all()


class Genes(SpliceGraphSQL):
    class _Gene(SpliceGraphType):
        def __init__(self, session, gene_id):
            self.session = session
            self.gene_id = gene_id

        def add(self, **kwargs):
            g = model.Gene(id=self.gene_id, **kwargs)
            self.session.add(g)
            return g

        @property
        def get(self):
            return self.session.query(model.Gene).get(self.gene_id)

        @property
        def exists(self):
            return self.session.query(exists().where(model.Gene.id == self.gene_id)).scalar()

    @property
    def genes(self):
        return self.session.query(model.Gene).all()

    def gene(self, gene_id):
        return self._Gene(self.session, gene_id)
