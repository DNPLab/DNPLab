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

    def test_fit_t2(self):

        t = np.linspace(0,10,100)
        y = 2*np.exp(-t/2)

        yError1 = np.zeros( len(y) )
        yError2 = np.zeros( len(y) )

        yError1[:30] = 0.05 * y[:30]
        yError2[30:] = 0.03 * y[30:]
        y = y + yError1 + yError2

        data = dnp.DNPData(y, ['tTest'], [t])
        out = dnp.fit_t2(data, 'tTest', [2,1] )
        out2 = dnp.fit_t2(data, 'tTest', [2,1,3] )

        with self.assertRaises(TypeError):
            dnp.fit_t2(data,[2,1,3] )

        with self.assertRaises(ValueError):
            dnp.fit_t2(data, "doesNotExist", [2,1,3] )

        #additional kwarg for curve_fit
        out3 = dnp.fit_t2(data, "tTest", [10,10,5], bounds=(4, np.inf) )
