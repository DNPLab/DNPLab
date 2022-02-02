import unittest
from numpy.testing import assert_allclose, assert_array_equal
import dnplab as dnp
import numpy as np

# FIX -> Should only use Fit functions, not import or processing functions
class dnpFit_tester(unittest.TestCase):
    def setUp(self):
        self.dim = "x"
        self.x = np.r_[-1:1:100j]
        self.y = dnp.lineshape.gaussian(self.x, 0, 0.1)
        self.p0 = (0, 1)
        self.data = dnp.DNPData(self.y, [self.dim], [self.x])

    def test_fit_function(self):
        out = dnp.fit(dnp.lineshape.gaussian, self.data, dim=self.dim, p0=self.p0)
        fit = out["fit"]
        print(fit)
        print(fit.coords["x"])

        assert_allclose(self.data.values, fit.values)
        assert_array_equal(self.data.coords["x"], fit.coords["x"])
