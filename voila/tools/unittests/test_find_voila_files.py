import unittest

from voila.tools.utils.find_files import get_voila_files


class TestStringMethods(unittest.TestCase):

    def test_hello_world(self):
        from os.path import expanduser
        home = expanduser("~")
        self.assertEqual(get_voila_files(directory=home, pattern="doesntexist"), [])
