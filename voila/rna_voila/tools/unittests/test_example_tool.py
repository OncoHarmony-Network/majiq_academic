import unittest

from rna_voila.tools.example_tool import hello_world


class TestExampleStringMethods(unittest.TestCase):
    """
    This is the unit view for the example tool.  There should be a unit view for each method/function within your code.
    """

    def test_hello_world(self):
        self.assertEqual(hello_world(), 'Hello, world!')
        self.assertEqual(hello_world(upper=True), 'HELLO, WORLD!')


class TestStringMethods(unittest.TestCase):
    """
    From Python documentation
    """

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)
