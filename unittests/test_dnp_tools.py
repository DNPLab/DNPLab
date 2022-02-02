import unittest
import dnplab as dnp
from numpy.testing import assert_array_equal
import numpy as np


class dnpTools_tester(unittest.TestCase):
    def setUp(self):
        x = np.r_[0:10]
        y = x ** 2.0
        self.data = dnp.DNPData(y, ["t2"], [x])

    def test_signal_to_noise(self):

        dnp.signal_to_noise(
            self.data,
            dim="t2",
        )

    def test_integrate(self):
        dnp.integrate(self.data, dim="t2")

    def test_mr_properties(self):

        info_1H = dnp.mr_properties("1H")
        self.assertEqual(info_1H, 26.7522128)

        info_1H = dnp.mr_properties("1H", 0.35)
        self.assertEqual(info_1H, 14902114.17018196)

        info_2H = dnp.mr_properties("2H", "qmom")
        self.assertEqual(info_2H, 0.286)

        info_6Li = dnp.mr_properties("6Li", "natAbundance")
        self.assertEqual(info_6Li, 7.59)

        info_6Li = dnp.mr_properties("6Li", "relSensitivity")
        self.assertEqual(info_6Li, 0.000645)
