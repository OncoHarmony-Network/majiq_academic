import argparse
from abc import ABC, abstractmethod


class ToolMethodNotImplemented(Exception):
    def __init__(self, cls, message):
        super(ToolMethodNotImplemented, self).__init__('{0}: {1}'.format(cls.__class__.__name__, message))


class MissingHelpMessageException(Exception):
    def __init__(self, cls):
        super(MissingHelpMessageException, self).__init__(cls.__class__.__name__)


class Tool(ABC):
    def __init__(self):
        if not hasattr(self, 'help'):
            raise MissingHelpMessageException(self)

    @abstractmethod
    def run(self, args):
        pass

    @abstractmethod
    def arguments(self):
        pass

    def get_parser(self, parent_args=()):
        return argparse.ArgumentParser(add_help=False, parents=parent_args)
