import os
import sys

__author__ = 'abarrera'


def we_are_frozen():
    # All of the modules are built-in to the interpreter, e.g., by cx_Freeze
    return hasattr(sys, "frozen")


def module_path():
    encoding = sys.getfilesystemencoding()
    return "See CJ, we need to take care of it"
    # if we_are_frozen():
    #     return os.path.dirname(str(sys.executable, encoding))
    # return os.path.dirname(str(__file__, encoding))
