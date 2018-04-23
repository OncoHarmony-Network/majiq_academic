import os
import unittest

from voila.tools import Tool


class Test(Tool):
    help = 'Run the available unit tests for all tools.'

    def arguments(self):
        return self.get_parser()

    def run(self, args):
        loader = unittest.TestLoader()
        start_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'unittests')
        suite = loader.discover(start_dir, '*')
        runner = unittest.TextTestRunner(verbosity=2)
        runner.run(suite)

