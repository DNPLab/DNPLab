import unittest
from numpy.testing import assert_array_equal
import dnplab as dnp
import numpy as np

#FIX -> Should only use Fit functions, not import or processing functions
class dnpFit_tester(unittest.TestCase):
    def setUp(self):
        self.x = np.r_[0:10:1024j]
        self.y = np.exp(-1*self.x/10)
        self.data = dnp.DNPData(self.y, ['t2'], [self.x])


    def test_fit_functions(self):
        self.out = efit.exponential_fit(self.data, type="T1")
        self.assertAlmostEqual(self.out.attrs["T1"], 2.140702947551208, places=4)

        self.out = efit.exponential_fit(self.data, type="T2")
        self.assertAlmostEqual(self.out.attrs["T2"], 1.0682212598985381, places=4)

        self.out = efit.exponential_fit(self.data, type="T2", stretched=True)
        self.assertAlmostEqual(self.out.attrs["T2"], 0.8938213879865939, places=4)

        self.out = efit.exponential_fit(self.data, type="mono")
        self.assertAlmostEqual(self.out.attrs["tau"], 2.140702798915825, places=4)
