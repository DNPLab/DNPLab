import unittest
import os
import dnplab as dnp
from numpy.testing import assert_array_equal
import numpy as np

import logging
import sys
import pathlib


# define logger for unittest output
#logger = logging.getLogger(__name__)
#logger.setLevel(logging.DEBUG)
#handler = logging.FileHandler('./log_test_dnp_tool.log')
#handler.setLevel(logging.DEBUG)
#logger.addHandler(handler)


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
        )
        self.data.attrs["nmr_frequency"] = 14.86e6

    def test_000_funcionality_signal_to_noise(self):
        """
        check only whether the function raises no error with DNPData input, not whether rsults are useful
        alot of simple tests lumped together

        Missing: check whether signal and noise are scalar values

        note that these tests are not really unittests but more integration tests
        """
        f = dnp.processing.signal_to_noise

        self.assertRaises(ValueError, f, self.data, (-0.1,0.015), (0.015, 0.022))

        data = dnp.fourier_transform(self.data)

        try:
            snr = f(data, (-0.1,0.015), (0.015, 0.022))
        except ValueError as e:
            self.fail("signal_to_noise reported ValueError {0}".format(e))
        self.assertTrue(not np.isnan(snr))

        snr = f(
            data,
            (-0.001,0.0015),
            (0.0015, 0.0022),
        )

    def test_001_using_different_dimensions(self):
        f = dnp.processing.signal_to_noise
        data = dnp.fourier_transform(self.data)

        # some input checks, just to check that no errors are thrown:
        snr = f(data, [(-0.001,0.001)], [(0.001, 0.0016)])
        snr = f(
            data, [(-0.001,0.0015)], [(0.0014, 0.0022)], remove_background=(-0.001640, -0.000200), deg=3
        )  # works with degree
        snr = f(
            data, [(-0.001,0.0015)], [(0.0014, 0.0022)], remove_background=(-0.001640, -0.000200)
        )  # works without degree
        snr = f(
            data, [(-0.001,0.0015)], [(0.0014, 0.0022)], remove_background=[(-0.001640, -0.000200)]
        )  # works with list as intended
        snr = f(
            data, [(-0.001215, 0.001041)], [(0.0014, 0.0022)], remove_background=[(-0.001640, -0.000200)]
        )
        snr = f(
            data,
            [(-0.001215, 0.001041)],
            [(0.0012, 0.00162)],
            remove_background=[(-0.0013001, -0.0005000)],
        )
        # with defaults
        snr = f(data)
        # with slices
        snr = f(data, slice(0, None), remove_background=[(-0.001640, -0.000200)])
        self.assertRaises(ValueError, f, data, [slice(0, None), (0.0015, 0.0022)])
        snr2 = f(data, (0, 0.001000), remove_background=[(-0.001640, -0.000200)])

        snr = f(
            data,
            slice(0, None),
            [slice(0, 100), slice(500, 600)],
            remove_background=[(-0.001640, -0.000200)],
        )

    def test_integrate(self):
        dnp.integrate(self.data, dim="t2")

    def test_mr_properties(self):
        info_1H = dnp.mr_properties("1H")
        self.assertEqual(info_1H, 26.7522128)

        info_1H = dnp.mr_properties("1H", 0.35)
        self.assertEqual(info_1H, 14902114.17018196)

        info_2H = dnp.mr_properties("2H", "qmom")
        self.assertEqual(info_2H, 0.286)

        info_6Li = dnp.mr_properties("6Li", "natAbundance")
        self.assertEqual(info_6Li, 7.59)

        info_6Li = dnp.mr_properties("6Li", "relSensitivity")
        self.assertEqual(info_6Li, 0.000645)
