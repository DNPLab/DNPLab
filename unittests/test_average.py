import unittest
import dnplab as dnp
from numpy.testing import assert_array_equal
import numpy as np

import logging

logger = logging.getLogger(__name__)


class dnpTools_tester(unittest.TestCase):
    def setUp(self):
        avg_dim = "Average"
        avg_x = np.r_[0:4]
        dim = "x"
        x = np.r_[0:10]
        y = np.array([x**2 * (i + 1) for i in avg_x])
        self.data = dnp.DNPData(y, [avg_dim, dim], [avg_x, x])
        self.avg_values = np.mean(y, axis=0)

    def test_average_function(self):
        avg_data = dnp.average(self.data)
        assert_array_equal(self.data.coords["x"], avg_data.coords["x"])
        assert_array_equal(avg_data.values, self.avg_values)


if __name__ == "__main__":
    unittest.main()
    pass
