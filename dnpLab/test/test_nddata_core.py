import unittest
from numpy.testing import assert_array_equal as assertArrayEqual
from dnpLab.core import nddata
import numpy as np
import random

test_dims = ['x', 'y', 'z', 'p', 'q', 'r']

class dnpLab_nddata_core_tester(unittest.TestCase):
    def setUp(self):
        self.dims = test_dims
        random.sample(test_dims, random.randint(1,len(test_dims)))

    def test_generating(self):
        for ix in range(1000):
            random_dims = random.sample(test_dims, random.randint(1,len(test_dims)))

            random_coords = [np.r_[0:random.randint(1,6)] for dim in random_dims]
            shape = [coord.size for coord in random_coords]

            random_values = np.random.randn(*shape)
            data = nddata.nddata_core(random_values, random_dims, random_coords)
            self.assertTrue(data._self_consistent())

if __name__ == '__main__':
    pass
