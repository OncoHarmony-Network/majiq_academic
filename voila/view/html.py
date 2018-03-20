import json
import os
import types
import uuid
from distutils.dir_util import copy_tree

import numpy
from jinja2 import Environment, FileSystemLoader, StrictUndefined
from voila import constants
from voila.constants import EXEC_DIR
from voila.utils import utils_voila
from voila.utils.voila_log import voila_log
from voila.vlsv import Het, HetGroup


class Html:
    def __init__(self, args):
        self.args = args
        self.voila_links = {}
        self._database_name = None

    def add_to_voila_links(self, lsv_dict, page_name):
        for gene_id in lsv_dict.keys():
            self.voila_links[gene_id] = '{0}#{1}'.format(os.path.join(constants.SUMMARIES_SUBFOLDER, page_name),
                                                         gene_id)

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
    def voila_links(lsv_dict, page_name):
        for gene_id in lsv_dict.keys():
            yield gene_id, '{0}#{1}'.format(os.path.join(constants.SUMMARIES_SUBFOLDER, page_name), gene_id)

    @classmethod
    def get_page_name(cls, args, index):
        try:
            output_html = cls.get_output_html(args, args.voila_file)
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

        env = Environment(extensions=["jinja2.ext.do"],
                          loader=FileSystemLoader(
                              [os.path.join(cls.get_template_dir(), 'summaries'), cls.get_template_dir()]),
                          undefined=StrictUndefined)
        env.filters.update({
            'to_json': to_json,
            'to_dict': to_dict,
        })

        return env

    @staticmethod
    def get_template_dir():
        return os.path.join(EXEC_DIR, "templates/")

    @classmethod
    def copy_static(cls, args, index=True):
        """
        Copy static files to output directory.
        :param index:
        :param args: command line arguments
        :return: None
        """
        log = voila_log()
        log.info("Copy static files from Voila sources")
        static_dir = os.path.join(cls.get_template_dir(), 'static')
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
