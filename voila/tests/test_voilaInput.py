from unittest import TestCase

from voila.io_voila import voila_input_from_hdf5

LSV_ID = 'ENSG00000000419:49551404-49551773:source'
QUANTIFY_FILE = 'tests/data/Ra_TN_Ra_T17.deltapsi_quantify.pickle.h5'


def getVI():
    return voila_input_from_hdf5(QUANTIFY_FILE, None)


class TestVoilaInput(TestCase):
    vi = getVI()

    def test_get_lsvs(self):
        self.assertEqual(self.vi.view_lsv_ids(), self.vi.lsvs)

    def test_add_lsv(self):
        new_lsv = self.vi.lsvs[0]
        self.assertEqual(new_lsv, self.vi.lsvs[0])
        lsvs_length = len(self.vi.lsvs)
        self.vi.add_lsv(new_lsv)
        self.assertEqual(len(self.vi.lsvs), lsvs_length + 1)
        self.assertTrue(new_lsv in self.vi.lsvs)

    def test_get_metainfo(self):
        with open('tests/data/metainfo.txt', 'w') as m:
            m.write(self.vi.metainfo)

    def test_samples_metainfo(self):
        self.fail()

    def test_encode_metainfo(self):
        self.fail()

    def test_exclude(self):
        self.fail()

    def test_decode_metainfo(self):
        self.fail()

    def test_to_hdf5(self):
        self.fail()

    def test_from_hdf5(self):
        self.fail()
