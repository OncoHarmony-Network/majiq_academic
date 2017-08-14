import unittest

from voila.tools.utils import find_files


class TestStringMethods(unittest.TestCase):
    def test_hello_world(self):
        from os.path import expanduser
        home = expanduser("~")
        self.assertEqual(find_files.find_voila_files(directory=home,
                                                     pattern="THIS_PATTERN_DOESNT_EXIST",
                                                     file_type="tsv"), [])
