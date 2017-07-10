import unittest
from voila.tools.utils.find_files import get_voila_files
from voila import tools
import os


# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


class TestStringMethods(unittest.TestCase):

    def test_hello_world(self):
        tool_dir = os.path.dirname(tools.__file__)
        self.assertEqual(get_voila_files(directory=tool_dir, pattern="doesntexist"), [])
