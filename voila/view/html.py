import errno
import json
import os
import types
import uuid
from abc import abstractmethod, ABC
from distutils.dir_util import copy_tree

import jinja2
import numpy
from jinja2 import Environment, FileSystemLoader, StrictUndefined

import voila
from voila import constants
from voila.api.view_splice_graph_sqlite import ViewSpliceGraph
from voila.constants import EXEC_DIR
from voila.processes import VoilaPool
from voila.utils import utils_voila
from voila.utils.voila_log import voila_log


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (numpy.int64, numpy.bool_, numpy.uint8, numpy.uint32, numpy.float32)):
            return obj.item()
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        elif isinstance(obj, numpy.bytes_):
            return obj.decode('utf-8')
        elif isinstance(obj, types.GeneratorType):
            return tuple(obj)
        return json.JSONEncoder.default(self, obj)


class Html(ABC):
    def __init__(self, args, matrix):
        self.args = args
        self.voila_links = {}
        self._database_name = None
        self.db_id = uuid.uuid4().hex
        self.matrix = matrix

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

    def old_copy_static(self, index=True):
        html_dir = self.get_template_dir()
        if index:
            copy_tree(os.path.join(html_dir, 'static'), os.path.join(self.args.output, 'static'))
        copy_tree(os.path.join(html_dir, 'static'), os.path.join(self.args.output, 'summaries', 'static'))

    def copy_static(self):
        """
        Copy static files to output directory.
        :param index:
        :param args: command line arguments
        :return: None
        """
        args = self.args
        log = voila_log()
        log.info("Copy static files from Voila sources")
        html_dir = os.path.join(self.get_template_dir(), '../html')
        copy_tree(os.path.join(html_dir, 'js'), os.path.join(args.output, 'js'))
        copy_tree(os.path.join(html_dir, 'css'), os.path.join(args.output, 'css'))
        copy_tree(os.path.join(html_dir, 'img'), os.path.join(args.output, 'img'))

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

    def render_html(self, index_html, summary_html):
        exec_dir = os.path.dirname(os.path.abspath(voila.__file__))
        template_dir = os.path.join(exec_dir, 'html')
        env = jinja2.Environment(extensions=["jinja2.ext.do"], loader=jinja2.FileSystemLoader(template_dir),
                                 undefined=jinja2.StrictUndefined)

        index_template = env.get_template(index_html)
        summary_template = env.get_template(summary_html)

        with open(os.path.join(self.args.output, 'index.html'), 'w') as index:
            index.write(index_template.render({
                'db_id': self.db_id
            }))

        with open(os.path.join(self.args.output, 'summary.html'), 'w') as summary:
            summary.write(summary_template.render({
                'db_id': self.db_id
            }))

    def render_dbs(self):
        log = voila_log()

        log.debug('Create metadata file')
        with open(os.path.join(self.args.output, 'metadata.js'), 'w') as f:
            with self.matrix(self.args) as h:
                metadata = json.dumps(h.view_metadata, cls=NumpyEncoder)
                f.write('new PouchDB(\'voila_gene_{}\').bulkDocs(['.format(self.db_id))
                f.write(metadata)
                f.write(',')
                f.write(']);')

                f.write('new PouchDB(\'voila_lsv_{}\').bulkDocs(['.format(self.db_id))
                f.write(metadata)
                f.write(',')
                f.write(']);')

                f.write('\n')

                f.write('const lsvs_arr = [')

                for lsv_id in h.view_lsv_ids():
                    lsv = h.psi(lsv_id)
                    f.write(json.dumps({
                        '_id': lsv.lsv_id,
                        'target': lsv.target,
                        'binary': lsv.binary,
                        'exon_skipping': lsv.exon_skipping,
                        'A5SS': lsv.a5ss,
                        'A3SS': lsv.a3ss,
                        'gene_id': lsv.gene_id
                    }))
                    f.write(',')
                f.write('];')

        with VoilaPool() as pool:
            with self.matrix(self.args) as h:
                gene_ids = list(h.view_gene_ids())
                chunked_gene_ids = Html.chunkify(gene_ids, pool.processes)

            for p in [pool.apply_async(self.dbs, (gene_ids,)) for gene_ids in chunked_gene_ids]:
                p.get()

    def dbs(self, gene_ids):
        log = voila_log()
        for gene_id in gene_ids:
            with open(os.path.join(self.args.output, '{}.js'.format(gene_id)), 'w') as f:
                with self.matrix(self.args) as h, ViewSpliceGraph(self.args) as sg:
                    metadata = h.view_metadata
                    exp_name = metadata['experiment_names']
                    lsv_ids = h.view_gene_lsvs(gene_id)

                    f.write('new PouchDB(\'voila_gene_{}\').bulkDocs(['.format(self.db_id))

                    log.debug('Write DB Gene ID: {}'.format(gene_id))

                    gene = sg.gene(gene_id)
                    gene_exp = sg.gene_experiment(gene, exp_name)
                    text = json.dumps(gene_exp)

                    f.write(text)
                    f.write(',')

                    del gene
                    del gene_exp
                    del text

                    f.write(']);')
                    f.write('\n')

                    if lsv_ids:
                        f.write('new PouchDB(\'voila_lsv_{}\').bulkDocs(['.format(self.db_id))

                        for lsv_id in lsv_ids:
                            log.debug('Write DB LSV ID: {}'.format(lsv_id))

                            lsv = h.lsv(lsv_id).get_all()
                            lsv_dict = dict(lsv)
                            text = json.dumps(lsv_dict, cls=NumpyEncoder)

                            f.write(text)
                            f.write(',')

                            del lsv
                            del lsv_dict
                            del text

                        f.write(']);')
