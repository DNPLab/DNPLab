import unittest
from numpy.testing import assert_array_equal
import dnplab as dnp
import numpy as np
import logging

# logging.basicConfig(filename='phase_debug.log', encoding='utf-8', level=logging.INFO)
logger = logging.getLogger(__name__)


class dnpNMR_tester(unittest.TestCase):
    def setUp(self):
        pts = 1024
        omega = 50 * dnp.pi
        tau = 0.1
        t2 = np.r_[0 : 1 : 1j * pts]
        y = np.exp(1j * t2 * omega) * np.exp(-1 * t2 / tau)
        self.data = dnp.DNPData(y, ["t2"], [t2])
        self.data.dnplab_attrs["frequency"] = 400e6

    def test_fourier_transform_functions(self):
        self.data = dnp.fourier_transform(self.data, zero_fill_factor=4)
        self.data = dnp.inverse_fourier_transform(self.data)


if __name__ == "__main__":
    unittest.main()
    pass
