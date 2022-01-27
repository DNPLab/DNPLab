import unittest
import dnplab as dnp
import os
from numpy.testing import assert_array_equal
import numpy as np


class save_h5_tester(unittest.TestCase):
    def setUp(self):
        self.x = np.r_[0:10]
        self.y = self.x**2
        self.data = dnp.DNPData(self.y, ['x'], [self.x])

        self.ws = {'data':self.data}


    def test_h5_save(self):

        dnp.save(
            self.data,
            os.path.join(".", "unittests", "test_save_DNPData.h5"),
            overwrite=True,
        )

        os.remove(os.path.join(".", "unittests", "test_save_DNPData.h5"))

    def test_h5_save_ws(self):

        dnp.save(
            self.ws,
            os.path.join(".", "unittests", "test_save_DNPData_dict.h5"),
            overwrite=True,
        )

        os.remove(os.path.join(".", "unittests", "test_save_DNPData_dict.h5"))


if __name__ == "__main__":
    pass
