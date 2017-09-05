import os
import unittest
from unittest.mock import patch

from voila.tools.ctrl_psi_dictionary import getMaxPSI, makeDictionary, CtrlPsiDictionaryTool


class TestCtrlPsiDictionaryMethods(unittest.TestCase):
    def test_get_max_psi(self):
        self.assertEqual(getMaxPSI('100'), float(100))
        self.assertEqual(getMaxPSI('100;100'), float(100))
        self.assertEqual(getMaxPSI('100;100;100'), float(100))
        for val in (1, 2, 10, 11, 12, 50, 75, 99999):
            self.assertEqual(getMaxPSI(';'.join(str(x) for x in range(val))), float(val - 1))

    def test_make_dictionary(self):
        with patch('sys.argv', [
            '--input', 'voila/tools/unittests/data/ctrl_psi_dictionary/GTEX_Filelist',
            '--output', '/tmp/'
        ]) as argv:

            args = CtrlPsiDictionaryTool().arguments()
            makeDictionary(args.parse_args(argv))
            tmp_data = open(os.path.join('/tmp', 'Ctrl_Psi_Dictionary')).read()

            with open('voila/tools/unittests/data/ctrl_psi_dictionary/Ctrl_PSI_Dictionary_Unit_Test_Output',
                      'r') as output:
                self.assertEqual(tmp_data, output.read())
