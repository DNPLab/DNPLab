import unittest
import os
import dnplab as dnp
from numpy.testing import assert_array_equal
import numpy as np

import logging
import sys
import pathlib

logger = logging.getLogger(__name__)


class dnpTools_tester(unittest.TestCase):
    def setUp(self):
        x = np.r_[0:10]
        y = x**2.0
        self.data = dnp.DNPData(y, ["t2"], [x])
        self.testdata = os.path.join(".", "data", "csv")
        p = pathlib.Path(self.testdata)
        self.data = dnp.io.load_csv.load_csv(
            p.joinpath("csv_example.csv"),
            skiprows=1,
            maxrows=1000,
            tcol=0,
            real=1,
            imag=3,
            convert_time=lambda x: float(x.replace(",", ".")) / 1e6,
        )
        self.data.attrs["nmr_frequency"] = 14.86e6

    def test_integrate(self):
        dnp.integrate(self.data, dim="t2")

    def test_cumulative_integrate(self):
        dnp.cumulative_integrate(self.data, dim="t2")

    def test_mr_properties(self):
        info_1H = dnp.mr_properties("1H")
        self.assertEqual(info_1H, 42577469.05766274)

        info_1H = dnp.mr_properties("1H", 0.35)
        self.assertEqual(info_1H, 14902114.170181958)

        info_2H = dnp.mr_properties("2H", "qmom")
        self.assertEqual(info_2H, 0.286)

        info_6Li = dnp.mr_properties("6Li", "natAbundance")
        self.assertEqual(info_6Li, 7.59)

        info_6Li = dnp.mr_properties("6Li", "relSensitivity")
        self.assertEqual(info_6Li, 0.000645)
