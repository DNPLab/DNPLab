import unittest
import numpy as np
from numpy.testing import assert_array_equal as assertArrayEqual
from .testing import get_gauss_3d
from dnplab.core.data import DNPData


class dnpDataTester(unittest.TestCase):
    def setUp(self):
        self.x, self.y, self.z, self.gauss_3d = get_gauss_3d(0.1)
        self.data = DNPData(self.gauss_3d, ["x", "y", "z"], [self.x, self.y, self.z])

    def test_DNPData(self):
        assertArrayEqual(self.data.coords["x"], self.x)
        assertArrayEqual(self.data.dims, ["x", "y", "z"])

    def test_coords_get_set(self):
        self.data.new_dim("r", np.r_[0:10])
        assertArrayEqual(self.data.dims, ["x", "y", "z", "r"])
        self.data.rename("r", "s")
        assertArrayEqual(self.data.dims, ["x", "y", "z", "s"])

    def test_coords_sort_reorder(self):
        self.data.reorder(["y", "z", "x"])
        assertArrayEqual(self.data.dims, ["y", "z", "x"])
        assertArrayEqual(self.data.coords["z"], self.z)
        self.data.sort_dims()
        assertArrayEqual(self.data.dims, ["x", "y", "z"])
        assertArrayEqual(self.data.coords["z"], self.z)


if __name__ == "__main__":
    pass
