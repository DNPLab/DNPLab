import unittest
from numpy.testing import assert_array_equal as assertArrayEqual
from .testing import get_gauss_3d
from dnpLab.dnpData import dnpData


class OdnpDataTester(unittest.TestCase):
    def setUp(self):
        self.x, self.y, self.z, self.gauss_3d = get_gauss_3d(0.1)
        self.dnpdata = dnpData(self.gauss_3d, [self.x, self.y, self.z], ['x', 'y', 'z'])

    def test_dnpdata(self):
        assertArrayEqual(self.dnpdata.getAxes('x'), self.x)
        assertArrayEqual(self.dnpdata.dims, ['x', 'y', 'z'])

    def test_coords_get_set(self):
        self.dnpdata.addAxes('r', range(10))
        assertArrayEqual(self.dnpdata.dims, ['x', 'y', 'z', 'r'])
        self.dnpdata.rename('r', 's')
        assertArrayEqual(self.dnpdata.dims, ['x', 'y', 'z', 's'])

    def test_coords_sort_reorder(self):
        self.dnpdata.reorder(['y', 'z', 'x'])
        assertArrayEqual(self.dnpdata.dims, ['y', 'z', 'x'])
        assertArrayEqual(self.dnpdata.getAxes('z'), self.z)
        self.dnpdata.sort()
        assertArrayEqual(self.dnpdata.dims, ['x', 'y', 'z'])
        assertArrayEqual(self.dnpdata.getAxes('z'), self.z)


if __name__ == '__main__':
    pass

