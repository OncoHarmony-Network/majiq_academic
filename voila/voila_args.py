import argparse
import errno
import os

from majiq.src.io_utils import create_if_not_exists
from voila.io_voila import Voila
from voila.splice_graphics import SpliceGraph


class VoilaCantFindFile(argparse.ArgumentTypeError):
    def __init__(self, value):
        super(VoilaCantFindFile, self).__init__('cannot find file "{0}"'.format(value))


class VoilaArgs:
    SPLICE_GRAPH_HELP = 'Path to splice graph file.  File should be named "splicegraph.hdf5".  This is required unless ' \
                        'the --no-html flag is set.'

    def get_parents(self):
        pass

    @staticmethod
    def check_dir(value):
        """
        check if directory exists.  If not, create it.
        :param value: directory path
        :return: value
        """
        create_if_not_exists(value)
        return value

    @staticmethod
    def check_list_file(value):
        """
        Take file, which is a newline separated list of values, and convert it to a list of strings.  Raise error if file
        doesn't exist.
        :param value: file path
        :return: return list of strings
        """
        try:
            with open(value, 'r') as f:
                return [line.strip() for line in f]
        except IOError as e:
            if e.errno == errno.ENOENT:
                raise VoilaCantFindFile(value)
            else:
                raise

    @staticmethod
    def check_file(value):
        """
        Check if file exists.
        :param value: file path
        :return:
        """
        if not os.path.isfile(value):
            raise VoilaCantFindFile(value)
        return value

    @classmethod
    def check_splice_graph_file(cls, value):
        cls.check_file(value)
        with SpliceGraph(value, 'r') as sg:
            sg.check_version()
        return value

    @classmethod
    def check_voila_file(cls, value):
        cls.check_file(value)
        with Voila(value, 'r') as v:
            v.check_version()
        return value

    @staticmethod
    def required_argument(*args, **kwargs):
        parser = args[0]
        required = parser.add_argument_group('required arguments')
        kwargs['required'] = True
        required.add_argument(*args[1:], **kwargs)

    @staticmethod
    def splice_graphs_args(cls):
        """
        Splice graphs specific arguments.
        :return: parser
        """
        parser_splice_graphs = argparse.ArgumentParser(add_help=False)

        parser_splice_graphs.add_argument('splice_graph',
                                          type=cls.check_splice_graph_file,
                                          help=cls.SPLICE_GRAPH_HELP)

        parser_splice_graphs.add_argument('--limit',
                                          type=int,
                                          default=20,
                                          help='Limit the number of splice graphs shown.  Default is 20.')

        return parser_splice_graphs

    @staticmethod
    def multiproccess_args():
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument(
            '--max-processes',
            type=int,
            help='Max number of processes used.  Processes will never exceed the number of available processes.'
        )
        return parser

    @staticmethod
    def get_parser():
        return argparse.ArgumentParser(add_help=False)

    @classmethod
    def base_args(cls):
        """
        Arguments common to all analysis types.
        :return:  parser
        """

        parser = cls.get_parser()

        parser.add_argument(
            '-l', '--logger',
            default=None,
            help='Path for log files.'
        )

        parser.add_argument(
            '-s', '--silent',
            action='store_true',
            help='Do not write logs to standard out.'
        )
        return parser

    @classmethod
    def output_args(cls):
        parser = cls.get_parser()
        cls.required_argument(
            parser,
            '-o', '--output',
            dest='output',
            type=cls.check_dir,
            help='Path for output directory.'
        )
        return parser

    @staticmethod
    def html_args():
        html_parser = argparse.ArgumentParser(add_help=False)
        html_parser.add_argument('--no-html',
                                 action='store_true',
                                 help='Do not write html files.')

        html_parser.add_argument('--no-tsv',
                                 action='store_true',
                                 help='Do not generate tab-separated values output file.')
        return html_parser

    @classmethod
    def voila_file_args(cls):
        """
        Arguments needs for analysis types who use the majiq quantifier file.
        :return: parser
        """

        voila_file_parser = argparse.ArgumentParser(add_help=False)

        voila_file_parser.add_argument('voila_file',
                                       type=cls.check_voila_file,
                                       help='Location of majiq\'s voila file.  File should end with ".voila".')
        cls.splice_graph = True
        voila_file_parser.add_argument('--splice-graph',
                                       type=cls.check_splice_graph_file,
                                       help=cls.SPLICE_GRAPH_HELP)

        voila_file_parser.add_argument('--gtf',
                                       action='store_true',
                                       help='Generate GTF (GFF2) files for LSVs.')

        voila_file_parser.add_argument('--gff',
                                       action='store_true',
                                       help='Generate GFF3 files for LSVs.')

        return voila_file_parser

    @classmethod
    def gene_search_args(cls):
        """
        Arguments for analysis types who search for specific genes.
        :return: parser
        """

        gene_search_parser = argparse.ArgumentParser(add_help=False)
        gene_search_parser.add_argument('--gene-names-file',
                                        dest='gene_names',
                                        type=cls.check_list_file,
                                        default=[],
                                        help='Location of file that contains a list of common gene names which should '
                                             'remain in the results.  One name per line.')

        gene_search_parser.add_argument('--gene-names',
                                        nargs='*',
                                        default=[],
                                        help='Common gene names, separated by spaces, which should remain in the results. '
                                             'e.g. GENE1 GENE2 ...')
        return gene_search_parser

    @classmethod
    def lsv_type_search_args(cls):
        """
        Arguments for analysis types who search for specific LSV types.
        :return: parser
        """
        lsv_search_parser = argparse.ArgumentParser(add_help=False)
        lsv_search_parser.add_argument('--lsv-types-file',
                                       type=cls.check_list_file,
                                       dest='lsv_types',
                                       help='Location of file that contains a list of LSV types which should remain in the '
                                            'results.  One type per line')

        lsv_search_parser.add_argument('--lsv-types',
                                       nargs='*',
                                       default=[],
                                       help='LSV types which should remain in the results')

        return lsv_search_parser

    @classmethod
    def lsv_id_search_args(cls):
        """
        Arguments for analysis types who search for specific LSV IDs.
        :return: parser
        """
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('--lsv-ids-file',
                            type=cls.check_list_file,
                            dest='lsv_ids',
                            help='Location of file that contains a list of LSV IDs which should remain in '
                                 'the results.  One ID per line.')

        parser.add_argument('--lsv-ids',
                            nargs='*',
                            default=[],
                            help='LSV IDs, separated by spaces, which should remain in the results.  e.g. '
                                 'LSV_ID1 LSV_ID2 ...')

        return parser


