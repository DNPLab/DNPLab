import unittest
from numpy.testing import assert_array_equal as assertArrayEqual
from .testing import get_gauss_3d
from odnpLab.odnpData import odnpData


class OdnpDataTester(unittest.TestCase):
    def setUp(self):
        self.x, self.y, self.z, self.gauss_3d = get_gauss_3d(0.1)
        self.odnpdata = odnpData(self.gauss_3d, [self.x, self.y, self.z], ['x', 'y', 'z'])

    def test_odnpdata(self):
        assertArrayEqual(self.odnpdata.get_axes('x'), self.x)
        assertArrayEqual(self.odnpdata.axesLabels, ['x', 'y', 'z'])

    def test_axes_get_set(self):
        self.odnpdata.add_axes('r', range(10))
        assertArrayEqual(self.odnpdata.axesLabels, ['x', 'y', 'z', 'r'])
        self.odnpdata.rename('r', 's')
        assertArrayEqual(self.odnpdata.axesLabels, ['x', 'y', 'z', 's'])

    def test_axes_sort_reorder(self):
        self.odnpdata.reorder(['y', 'z', 'x'])
        assertArrayEqual(self.odnpdata.axesLabels, ['y', 'z', 'x'])
        assertArrayEqual(self.odnpdata.get_axes('z'), self.z)
        self.odnpdata.sort()
        assertArrayEqual(self.odnpdata.axesLabels, ['x', 'y', 'z'])
        assertArrayEqual(self.odnpdata.get_axes('z'), self.z)


if __name__ == '__main__':
    pass

