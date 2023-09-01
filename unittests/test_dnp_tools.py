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
            convert_time=lambda x: float(x.replace(",", ".")) / 1e6,
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

        self.assertRaises(ValueError, f, self.data, (300, 400), (500, 600))

        data = dnp.fourier_transform(self.data)

        try:
            snr = f(data, (300, 400), (500, 600))
        except ValueError as e:
            self.fail("signal_to_noise reported ValueError {0}".format(e))
        self.assertTrue(not np.isnan(snr))

        snr = f(
            data,
            (300, 400),
            (500, 600),
        )

    def test_001_using_different_dimensions(self):
        f = dnp.processing.signal_to_noise
        data = dnp.fourier_transform(self.data)

        # some input checks, just to check that no errors are thrown:
        snr = f(data, [(300, 400)], [(500, 600)])
        snr = f(
            data, [(300, 400)], [(500, 600)], remove_background=(100, 200), deg=3
        )  # works with degree
        snr = f(
            data, [(300, 400)], [(500, 600)], remove_background=(100, 200)
        )  # works without degree
        snr = f(
            data, [(300, 400)], [(500, 600)], remove_background=[(100, 200)]
        )  # works with list as intended
        snr = f(
            data, [(-121.5, 104.1)], [(632.5, 1264.2)], remove_background=[(100, 200)]
        )
        snr = f(
            data,
            [(-121.5, 104.1)],
            [(632.5, 1264.2)],
            remove_background=[(-1300.1, -500.0)],
        )
        # with defaults
        snr = f(data)
        # with slices
        snr = f(data, slice(0, None), remove_background=[(100, 200)])
        #self.assertRaises(ValueError, f, data, [slice(0, None), (100, 300)])
        # with more than one signal region:
        snr = f( data, [slice(0, None), (100, 300)])
        self.assertEqual(len(snr),2)
        self.assertEqual(len(snr[0]),1)
        self.assertEqual(len(snr[1]),1)

        snr2 = f(data, (0, 1000), remove_background=[(100, 200)])

        snr = f(
            data,
            slice(0, None),
            [slice(0, 100), slice(500, 600)],
            remove_background=[(100, 200)],
        )

    def test_002_SNR_on_higherDimensionalData(self):
        import dnplab as dnp

        coords3 = [np.arange(0, 100), np.arange(0, 20), np.arange(0, 40)]
        data3 = np.random.random((100, 20, 40))
        DNPObj3 = dnp.DNPData(data3, ["t2", "t3", "t4"], coords3)
        f = dnp.processing.signal_to_noise
        snr = f(DNPObj3, [(10, 20),(30,40),(50,60)], [(80, 90)],dim="t2")
        self.assertEqual(len(snr),3)
        self.assertEqual(len(snr[0]),800)
        self.assertEqual(len(snr[1]),800)
        self.assertEqual(len(snr[2]),800)

    def test_integrate(self):
        dnp.integrate(self.data, dim="t2")

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
