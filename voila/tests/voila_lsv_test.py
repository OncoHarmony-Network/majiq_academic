__author__ = 'abarrera'
from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises
from nose import with_setup
import os
try:
    import cPickle as pkl
except ImportError:
    import pickle as pkl
from voila.vlsv import VoilaLsv
import numpy as np

LSVGRAPHIC_FILENAME = os.path.join(os.path.dirname(__file__), 'data/lsv_graphic.pickle')


class TestVoilaLsv(object):
    @classmethod
    def setup_class(cls):
        """This method is run once for each class before any tests are run"""
        sg = pkl.load(open(LSVGRAPHIC_FILENAME))
        binsl = []
        num_juncs = 3
        for ii in xrange(num_juncs):
            binsl.append(np.random.rand(40))
            binsl[-1] /= binsl[-1].sum()

        cls.psi1 = []
        for ii in xrange(num_juncs):
            cls.psi1.append(np.random.rand(20))
            cls.psi1[-1] /= cls.psi1[-1].sum()

        cls.psi2 = []
        for ii in xrange(num_juncs):
            cls.psi2.append(np.random.rand(20))
            cls.psi2[-1] /= cls.psi2[-1].sum()

        cls.sg = sg
        cls.binsl = binsl

    # @classmethod
    # def teardown_class(cls):
    #     """This method is run once for each class _after_ all tests are run"""

    def setUp(self):
        """This method is run once before _each_ test method is executed"""

    def teardown(self):
        """This method is run once after _each_ test method is executed"""

    def test_init_psi(self):
        a = VoilaLsv(TestVoilaLsv.binsl, self.sg)
        assert_equal(a.is_delta_psi(), False)
        assert_not_equal(a.is_delta_psi(), True)

    def test_init_dpsi(self):
        a = VoilaLsv(TestVoilaLsv.binsl, self.sg, psi1=TestVoilaLsv.psi1, psi2=TestVoilaLsv.psi2)
        assert_equal(a.is_delta_psi(), True)
        assert_not_equal(a.is_delta_psi(), False)

    # def test_return_true(self):
    #     a = VoilaLsv(TestVoilaLsv.binsl, self.sg)
    #     assert_equal(a.return_true(), True)
    #     assert_not_equal(a.return_true(), False)
    #
    # def test_raise_exc(self):
    #     a = A()
    #     assert_raises(KeyError, a.raise_exc, "A value")
    #
    # @raises(KeyError)
    # def test_raise_exc_with_decorator(self):
    #     a = A()
    #     a.raise_exc("A message")