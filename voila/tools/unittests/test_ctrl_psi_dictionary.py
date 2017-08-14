import argparse
import unittest
from unittest import mock
from unittest.mock import patch

import sys

from voila.tools.ctrl_psi_dictionary import CtrlPsiDictionaryTool


class TestStringMethods(unittest.TestCase):
    def test_some_test(self):
        parser = CtrlPsiDictionaryTool().get_parser()
        sys.argv.append('--help')

        print(parser.parse_known_args())

