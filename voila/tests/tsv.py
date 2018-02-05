import csv
import os
import unittest
from itertools import chain


class TestStringMethods(unittest.TestCase):
    def tsv_files(self):
        for f in (f for f in os.listdir('./') if f.endswith('.tsv')):
            tsvfile = open(f, 'r')
            yield csv.DictReader(tsvfile, delimiter='\t')
            tsvfile.close()

    def test_junction_count(self):
        for tsv in self.tsv_files():
            for row in tsv:
                e_psi = (k for k in row.keys() if k.endswith('E(PSI)'))
                per_junc = (k for k in row.keys() if k.endswith('per LSV junction'))

                junctions = int(row['Num. Junctions'])
                self.assertEqual(junctions, len(row['Junctions coords'].split(';')))
                for x in chain(e_psi, per_junc):
                    self.assertEqual(junctions, len(row[x].split(';')))


if __name__ == '__main__':
    unittest.main()
