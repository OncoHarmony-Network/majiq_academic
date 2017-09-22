import unittest
from voila.tools.hello_caleb import hello_world


# Caleb Matthew Radens
# radlinsky@gmail.com


__author__ = 'cradens'


class TestStringMethods(unittest.TestCase):

    def test_hello_world(self):
        self.assertEqual(hello_world(), 'Hello, world!')
        self.assertEqual(hello_world(capitalize=True), 'HELLO, WORLD!')