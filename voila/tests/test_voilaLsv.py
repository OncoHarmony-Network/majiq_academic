from unittest import TestCase
from uuid import uuid4

import h5py

from voila.tests.test_voilaInput import LSV_ID, QUANTIFY_FILE
from voila.vlsv import VoilaLsv


def getLSV():
    with h5py.File(QUANTIFY_FILE, 'r') as h:
        return VoilaLsv((), None).from_hdf5(h['lsvs'][LSV_ID])


class TestVoilaLsv(TestCase):
    lsv = getLSV()

    def test_to_gff3(self):
        with open('tests/data/lsv.gff3', 'r') as g:
            self.assertEqual(VoilaLsv.to_gff3(getLSV()), g.read())

    def test_init_categories(self):
        self.assertEqual(VoilaLsv.init_categories(self.lsv.get_type()), self.lsv.get_categories())

    def test_is_delta_psi(self):
        self.assertTrue(self.lsv.is_delta_psi())

    def test_get_lsv_graphic(self):
        self.assertEqual(self.lsv.get_lsv_graphic(), self.lsv.lsv_graphic)

    def test_get_chrom(self):
        self.assertEqual(self.lsv.get_chrom(), self.lsv.lsv_graphic.chrom)

    def test_get_strand(self):
        self.assertEqual(self.lsv.get_strand(), self.lsv.lsv_graphic.strand)

    def test_set_chrom(self):
        uuid = uuid4().hex
        self.lsv.set_chrom(uuid)
        self.assertEqual(self.lsv.lsv_graphic.chrom, uuid)

    def test_set_strand(self):
        self.lsv.set_strand('-')
        self.assertEqual(self.lsv.lsv_graphic.strand, '-')
        self.lsv.set_strand('+')
        self.assertEqual(self.lsv.lsv_graphic.strand, '+')

    def test_set_id(self):
        uuid = uuid4().hex
        self.lsv.set_id(uuid)
        self.assertEqual(self.lsv.lsv_graphic.id, uuid)

    def test_get_id(self):
        self.assertEqual(self.lsv.get_id(), self.lsv.lsv_graphic.id)

    def test_get_gene_name(self):
        self.assertEqual(self.lsv.get_gene_name(), self.lsv.lsv_graphic.name)

    def test_set_type(self):
        uuid = uuid4().hex
        self.lsv.set_type(uuid)
        self.assertEqual(self.lsv.lsv_graphic.type, uuid)

    def test_get_type(self):
        self.assertEqual(self.lsv.get_type(), self.lsv.lsv_graphic.type)

    def test_get_bins(self):
        self.assertEqual(self.lsv.get_bins(), self.lsv.bins)

    def test_set_means(self):
        l = [uuid4().hex for _ in range(4)]
        self.lsv.set_means(l)
        self.assertTrue(self.lsv.means, l)

    def test_get_means(self):
        self.assertTrue(self.lsv.get_means(), self.lsv.means)

    def test_set_coords(self):
        l = [uuid4().hex for _ in range(2)]
        self.lsv.set_coords(l)
        self.assertTrue(self.lsv.lsv_graphic.coords, l)

    def test_get_coords(self):
        self.assertEqual(self.lsv.get_coordinates(), self.lsv.lsv_graphic.coords)

    def test_get_variances(self):
        self.assertEqual(self.lsv.get_variances(), self.lsv.variances)

    def test_set_excl_incl(self):
        s = [uuid4().hex for _ in range(4)]
        self.lsv.set_excl_incl(s)
        self.assertEqual(self.lsv.excl_incl, s)

    def test_get_excl_incl(self):
        self.assertEqual(self.lsv.get_excl_incl(), self.lsv.excl_incl)

    def test_get_extension(self):
        self.assertEqual(self.lsv.get_extension(), [49551404, 49552799])

    def test_get_categories(self):
        self.assertEqual(self.lsv.get_categories(), self.lsv.categories)

    def test_njuncs(self):
        self.assertEqual(self.lsv.njuncs(), self.lsv.categories['njuncs'])

    def test_nexons(self):
        self.assertEqual(self.lsv.nexons(), self.lsv.categories['nexons'])

    def test_categories2css(self):
        self.assertEqual(self.lsv.categories2css(), 'ir target')

    def test_set_bed12(self):
        uuid = uuid4().hex
        self.lsv.set_bed12(uuid)
        self.assertEqual(self.lsv.bed12_str, uuid)

    def test_get_bed12(self):
        self.assertEqual(self.lsv.get_bed12(), self.lsv.bed12_str)

    def test_get_gff3(self):
        self.assertEqual(self.lsv.get_gff3(), VoilaLsv.to_gff3(self.lsv))

    def test_to_JSON(self):
        lsv = getLSV()

        with open('tests/data/lsv.json', 'r') as json:
            self.assertEqual(lsv.to_json(), json.read())

    def test_is_lsv_changing(self):
        self.assertFalse(self.lsv.is_lsv_changing(0.5))

    def test_exclude(self):
        self.assertEqual(self.lsv.exclude(), ['lsv_graphic', 'categories', 'bins', 'psi1', 'psi2'])

    def test_to_hdf5(self):
        # see Voila Input
        self.assertTrue(True)

    def test_from_hdf5(self):
        # See Voila Input
        self.assertTrue(True)
