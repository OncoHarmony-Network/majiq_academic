import unittest

from voila.tools.utils import print_hello


class TestPrintHello(unittest.TestCase):
    def test_func_inside_print_hello_file(self):
        self.assertEqual(print_hello.foo(), 'hello')