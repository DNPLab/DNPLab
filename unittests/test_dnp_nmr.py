import unittest
from numpy.testing import assert_array_equal
import dnplab as dnp
import numpy as np


class dnpNMR_tester(unittest.TestCase):
    def setUp(self):
        pts = 1024
        omega = 50*np.pi
        tau = 0.1
        t2 = np.r_[0:1:1j*pts]
        y = np.exp(1j*t2*omega) * np.exp(-1 * t2 / tau)
        self.data = dnp.DNPData(y, ['t2'], [t2])
        self.data.attrs['nmr_frequency'] = 400e6

    def test_basic_nmr_processing(self):

        data = dnp.remove_background(self.data)

        data = dnp.left_shift(data)

        data = dnp.apodize(data, lw = 1)

        data = dnp.fourier_transform(data)

        data = dnp.autophase(data)

    def test_calc_enhancement(self):
        pass


class dnpNMR_tester_sim(unittest.TestCase):
    def setUp(self):
        p1 = np.array([0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0])
        p2 = np.array([0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0, 0, 0])
        p3 = np.array([0, 0, 0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0])
        self.data = dnp.DNPData(
            np.array([p1, p2, p3]).T,
            ["x", "t2"],
            [np.arange(0, len(p1)), np.arange(0, 3)],
        )

    def test_align(self):
        self.aligned_data = dnp.align(self.data, dim="x")
        assert_array_equal(
            self.aligned_data.values,
            np.array(
                [
                    [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0],
                    [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0],
                    [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0],
                ]
            ).T,
        )
