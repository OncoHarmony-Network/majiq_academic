import errno
import json
import os
import types
import uuid
from abc import abstractmethod, ABC
from distutils.dir_util import copy_tree

import numpy
from jinja2 import Environment, FileSystemLoader, StrictUndefined

from voila import constants
from voila.api.view_splice_graph import ViewSpliceGraph
from voila.constants import EXEC_DIR
from voila.utils import utils_voila
from voila.utils.voila_log import voila_log
from voila.utils.voila_pool import VoilaPool
from voila.vlsv import Het, HetGroup


class Html(ABC):
    def __init__(self, args):
        self.args = args
        self.voila_links = {}
        self._database_name = None

    def database_name(self):
        if self._database_name is None:
            self._database_name = 'voila_{}'.format(uuid.uuid4().hex)

        return self._database_name

    @staticmethod
    def get_summaries_subfolder(args):
        summaries_subfolder = os.path.join(args.output, constants.SUMMARIES_SUBFOLDER)
        utils_voila.create_if_not_exists(summaries_subfolder)
        return summaries_subfolder

    @staticmethod
    def get_voila_links(lsv_dict, page_name):
        for gene_id in lsv_dict.keys():
            yield gene_id, '{0}#{1}'.format(os.path.join(constants.SUMMARIES_SUBFOLDER, page_name), gene_id)

    @classmethod
    def get_page_name(cls, args, index):
        # todo: this shouldn't require the voila file name to work. Maybe just the suffix.
        try:
            output_html = cls.get_output_html(args, args.voila_file[0])
        except AttributeError:
            output_html = cls.get_output_html(args, args.splice_graph)
        return '{0}_{1}'.format(index, output_html)

    @classmethod
    def get_next_page(cls, args, index, page_count):
        if index + 1 == page_count:
            return None
        else:
            return cls.get_page_name(args, index + 1)

    @classmethod
    def get_prev_page(cls, args, index):
        if index - 1 < 0:
            return None
        else:
            return cls.get_page_name(args, index - 1)

    @staticmethod
    def chunkify(lst, n):
        for i in range(n):
            yield lst[i::n]

    @classmethod
    def get_env(cls):
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, (numpy.int64, numpy.bool_, numpy.uint8, numpy.uint32)):
                    return obj.item()
                elif isinstance(obj, numpy.ndarray):
                    return obj.tolist()
                elif isinstance(obj, numpy.bytes_):
                    return obj.decode('utf-8')
                elif isinstance(obj, (Het, HetGroup)):
                    return dict(obj)
                elif isinstance(obj, types.GeneratorType):
                    return tuple(obj)
                return json.JSONEncoder.default(self, obj)

        def to_json(value):
            return json.dumps(value, cls=NumpyEncoder)

        def to_dict(value):
            return dict(value)

        env = Environment(extensions=["jinja2.ext.do"], loader=FileSystemLoader(cls.get_template_dir()),
                          undefined=StrictUndefined)
        env.filters.update({
            'to_json': to_json,
            'to_dict': to_dict,
        })

        return env

    @staticmethod
    def get_template_dir():
        return os.path.join(EXEC_DIR, "templates/")

    def copy_static(self, index=True):
        """
        Copy static files to output directory.
        :param index:
        :param args: command line arguments
        :return: None
        """
        args = self.args
        log = voila_log()
        log.info("Copy static files from Voila sources")
        static_dir = os.path.join(self.get_template_dir(), 'static')
        if index:
            copy_tree(static_dir, os.path.join(args.output, 'static'))
        copy_tree(static_dir, os.path.join(args.output, constants.SUMMARIES_SUBFOLDER, 'static'))

    @staticmethod
    def get_output_html(args, file_name=None):
        """
        Get output html file name.
        :param args: command line arguments
        :param file_name: input file name
        :return:
        """

        if file_name:
            return '{0}.html'.format(os.path.splitext(os.path.split(file_name)[1])[0])
        else:
            return '{0}.html'.format(args.type_analysis.replace("-", "_"))

    @staticmethod
    def table_marks_set(size):
        """
        Calculate the number of elements to show in LSV tables.

        :param size: total number of LSVs.
        :return: set of total number of elements to show.
        """
        ideal_set = (10, 20, 50, 100)
        for index, value in enumerate(ideal_set):
            if size < value:
                return ideal_set[0:index]
        return ideal_set

    @abstractmethod
    def render_dbs(self):
        pass

    def create_gene_db(self, gene_ids, view_matrix, data_type):
        template = self.get_env().get_template('gene_db_template.html')
        log = voila_log()
        args = self.args
        metadata = self.view_metadata
        experiment_names = metadata['experiment_names']
        with ViewSpliceGraph(args) as sg, view_matrix(args) as m:
            for gene in sg.genes(gene_ids):
                log.debug('creating {}'.format(gene.id))

                with open(os.path.join(args.output, 'db', '{}.js'.format(gene.id)), 'w') as html:
                    html.write(
                        template.render(
                            gene=gene.get_experiment(experiment_names),
                            lsvs=tuple(getattr(m, data_type)(lsv_id) for lsv_id in m.view_gene_lsvs(gene.id))
                        )
                    )

    def create_metadata_db(self, view_matrix):
        args = self.args
        template = self.get_env().get_template('metadata_db_template.html')
        with view_matrix(args) as vm:
            with open(os.path.join(args.output, 'db', 'metadata.js'), 'w') as html:
                html.write(
                    template.render(metadata=vm.view_metadata)
                )

    def create_db_files(self, view_matrix, data_type):
        args = self.args

        with view_matrix(args) as m:
            gene_ids = tuple(m.view_gene_ids())

        try:
            os.makedirs(os.path.join(args.output, 'db'))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        self.create_metadata_db(view_matrix)

        multiple_results = []

        vp = VoilaPool()
        for genes in self.chunkify(gene_ids, vp.processes):
            multiple_results.append(vp.apply_async(self.create_gene_db, (genes, view_matrix, data_type)))

        [res.get() for res in multiple_results]

    @abstractmethod
    def create_summary(self, paged):
        pass

    def create_summaries(self, view_matrix):
        log = voila_log()
        log.info('Render Delta PSI HTML summaries')

        args = self.args

        with view_matrix(args) as m:
            log.debug('Get paginated genes')
            paged_genes = tuple(m.paginated_genes())

        log.debug('Set up pool')
        multiple_results = []
        with VoilaPool() as vp:
            for paged in self.chunkify(tuple(enumerate(paged_genes)), vp.processes):
                multiple_results.append(
                    vp.pool.apply_async(self.create_summary, (paged,)))

            log.debug('get pool results')
            for res in multiple_results:
                self.voila_links.update(res.get())
