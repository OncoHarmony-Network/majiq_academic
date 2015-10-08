__author__ = 'abarrera'
from nose.tools import assert_equal
from nose.tools import assert_not_equal
from nose.tools import assert_raises
from nose.tools import raises
import os
try:
    import cPickle as pkl
except ImportError:
    import pickle as pkl
from voila.vlsv import VoilaLsv
import numpy as np

LSVGRAPHIC_FILENAME = os.path.join(os.path.dirname(__file__), 'data/lsv_splicegraph.pickle')


class VoilaLsvTest(object):
    @classmethod
    def setup_class(cls):
        """This method is run once for each class before any tests are run"""
        sg = pkl.load(open(LSVGRAPHIC_FILENAME))
        binsl = []
        for ii in xrange(3):
            binsl.append(np.random.rand(40))
            binsl[-1] /= binsl[-1].sum()

        test_instance = cls()
        test_instance.sg = sg
        test_instance.binsl = binsl

        return test_instance

    @classmethod
    def teardown_class(cls):
        """This method is run once for each class _after_ all tests are run"""

    def setUp(self):
        """This method is run once before _each_ test method is executed"""

    def teardown(self):
        """This method is run once after _each_ test method is executed"""

    def initTest(self):
        a = VoilaLsv(self.binsl, self.sg)
        assert_equal(a.psi1, None)
        assert_not_equal(a.bins, None)

    # def test_return_true(self):
    #     a = A()
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