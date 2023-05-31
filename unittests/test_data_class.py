import unittest
import numpy as np
from numpy.testing import assert_array_equal as assertArrayEqual
from .testing import get_gauss_3d
import dnplab as dnp


class dnpDataTester(unittest.TestCase):
    def setUp(self):
        self.x, self.y, self.z, self.gauss_3d = get_gauss_3d(0.1)
        self.data = dnp.DNPData(
            self.gauss_3d, ["x", "y", "z"], [self.x, self.y, self.z]
        )

        self.size = 10
        self.d = np.random.random((self.size, self.size))

        self.ax = np.arange(self.size)
        self.Ddata = dnp.DNPData(self.d, dims=["1", "2"], coords=[self.ax, self.ax])

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

    def test_000_checkProcDim(self):
        f = np.sin
        f2 = np.max

        tdata = f2(self.Ddata, axis="1")

        if not ("numpy.amax" in tdata.proc_attrs[0]):
            self.fail(
                "Assertion that np.max is in proc_attrs is False, proc_attrs:{0}".format(
                    tdata.proc_attrs
                )
            )
        if not (tdata.proc_attrs[-1][1]["axis"] == ["1"]):
            self.fail(
                "Assertion that axis is {0} is False, tdata.proc_attrs[-1]:{1}".format(
                    "1", tdata.proc_attrs[-1][1]
                )
            )

        tdata2 = f2(self.Ddata)

        try:
            float(tdata2)
        except ValueError:
            self.fail("tdata2 is not something like a float: {0}".format(tdata2))

    def test_001_checkkwargsConfig(self):
        data1 = dnp.DNPData(self.gauss_3d, ["x", "y", "z"], [self.x, self.y, self.z])
        self.assertEqual(5, data1.max_print_attrs)
        self.assertEqual(False, data1.print_values)
        data2 = dnp.DNPData(
            self.gauss_3d,
            ["x", "y", "z"],
            [self.x, self.y, self.z],
            max_print_attrs=10,
            print_values=True,
        )
        self.assertEqual(10, data2.max_print_attrs)
        self.assertEqual(True, data2.print_values)


if __name__ == "__main__":
    pass
