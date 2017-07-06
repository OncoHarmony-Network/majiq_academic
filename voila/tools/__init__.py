import importlib
import inspect
import os

from voila.tools.tool import Tool
from voila.voila_args import VoilaArgs


class ToolParserEmptyException(Exception):
    def __init__(self, tool_name):
        m = 'Parser is empty for tool {0}.  Check to see if its arguments method returns a parser.'.format(tool_name)
        super(ToolParserEmptyException, self).__init__(m)


class ToolClassNotFoundException(Exception):
    pass


class TooManyToolClasses(Exception):
    pass


class Tools(VoilaArgs):
    filter = ('__init__.py', 'tool.py')
    tool_subclass = Tool
    module = 'voila.tools'
    tool_dir = os.path.dirname(os.path.realpath(__file__))

    def __init__(self, args):
        self.get_tool(args.tool).run(args)

    @classmethod
    def add_arguments(cls, parser):
        subparser = parser.add_subparsers(dest='tool')
        subparser.required = True
        for tool_name in cls.tool_names():
            tool = cls.get_tool(tool_name)
            tool_parser = tool.arguments()
            base_parser = cls.base_args()
            if not tool_parser:
                raise ToolParserEmptyException(tool_name)
            subparser.add_parser(tool_name, parents=[base_parser, tool_parser], help=tool.help)

    @classmethod
    def validate_tool_file(cls, file_name):
        return file_name.endswith('.py') and file_name not in cls.filter

    @classmethod
    def tool_names(cls):
        for tool_file_name in filter(lambda x: cls.validate_tool_file(x), os.listdir(cls.tool_dir)):
            yield tool_file_name.split('.')[0]

    @classmethod
    def get_members(cls, tool_name):
        module = importlib.import_module('{0}.{1}'.format(cls.module, tool_name))
        for member in inspect.getmembers(module, inspect.isclass):
            if member[0] != 'Tool' and issubclass(member[1], cls.tool_subclass):
                yield member[1]

    @classmethod
    def get_tool(cls, tool_name):
        members = tuple(cls.get_members(tool_name))

        if len(members) > 1:
            m = 'Found more than one {0} class in {1}'.format(cls.tool_subclass.__name__, tool_name)
            raise TooManyToolClasses(m)

        if len(members) == 0:
            m = 'No subclass of {0} was found in {1}'.format(cls.tool_subclass.__name__, tool_name)
            raise ToolClassNotFoundException(m)

        return members[0]()

    @classmethod
    def arg_parents(cls):
        return ()
