import unittest

from voila.tools.find_voila_files import get_voila_files


class TestStringMethods(unittest.TestCase):

    def test_hello_world(self):
        self.assertEqual(get_voila_files(directory="./", pattern="doesntexist"), [])
