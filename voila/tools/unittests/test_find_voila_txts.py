import unittest

from voila.tools.find_voila_txts import find_voila_txts


class TestStringMethods(unittest.TestCase):

    def test_hello_world(self):
        self.assertEqual(find_voila_txts(directory="./", pattern="doesntexist"), [])
