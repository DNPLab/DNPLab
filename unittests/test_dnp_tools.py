import unittest
import os
import dnplab as dnp
from numpy.testing import assert_array_equal
import numpy as np

import logging
import sys
import pathlib


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

    def test_000_funcionality_signal_to_noise(self):
        """
        check only whether the function raises no error with DNPData input, not whether rsults are useful
        alot of simple tests lumped together

        Missing: check whether signal and noise are scalar values
        """
        f = dnp.processing.signal_to_noise

        self.assertRaises(ValueError, f, self.data, (300, 400), (500, 600))

        data = dnp.fourier_transform(self.data)

        try:
            data, snr = f(data, (300, 400), (500, 600))
        except ValueError as e:
            self.fail("signal_to_noise reported ValueError {0}".format(e))
        self.assertTrue(not np.isnan(snr))

        dat, snr, signal, noise = f(
            data, (300, 400), (500, 600), fullOutput=True, detrend=False
        )

        self.assertTrue(not np.isnan(snr))
        self.assertTrue(not np.isnan(signal))
        self.assertTrue(not np.isnan(noise))

        # some input checks:
        dat, snr = f(data, [(300, 400)], [(500, 600)])
        self.assertTrue(data._values.size>0) #check for existing data
        dat, snr = f(
            data, [(300, 400)], [(500, 600)], remove_background=(100, 200), deg=3
        )  # works with degree
        dat, snr = f(
            data, [(300, 400)], [(500, 600)], remove_background=(100, 200)
        )  # works without degree
        dat, snr = f(
            data, [(300, 400)], [(500, 600)], remove_background=[(100, 200)]
        )  # works with list as intended

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
