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
        snr2 = f(data, (0, 1000), remove_background=[(100, 200)])

        snr = f(
            data,
            slice(0, None),
            [slice(0, 100), slice(500, 600)],
            remove_background=[(100, 200)],
        )

    def test_002_check_if_data_is_unmodified(self):
        f = dnp.processing.signal_to_noise
        data = dnp.fourier_transform(self.data)
        s1 = data.shape

        # some input checks, just to check that no errors are thrown:
        snr = f(data, [(300, 400)], [(500, 600)])
        snr = f(
            data, [(300, 400)], [(500, 600)], remove_background=(100, 200), deg=3
        )  # works with degree
        snr = f(
            data, [(300, 400)], [(500, 600)], remove_background=(100, 200)
        )  # works without degree
        snr = f(data, [(300, 400)], [(500, 600)], remove_background=[(100, 200)])
        self.assertFalse(np.any(np.isnan(snr)))
        s2 = data.shape
        self.assertEqual(s1, s2)

    def test_003_multiple_regions_behaviour(self):
        from dnplab.core.base import ABCData
        ABCData.add_proc_attrs=lambda *args,**kwargs: '1'
        d1 = 512
        d2 = 10

        d = np.random.random((d1, d2))
        ax1 = np.arange(d1)
        ax2=np.arange(d2)
        Ddata = ABCData(d, dims=["f2", "2"], coords=[ax1, ax2])

        f = dnp.processing.signal_to_noise
        snr = f(Ddata,[(10,20),(50,60),(100,200)],[(400,500)])

        self.assertEqual(len(snr),3)
        self.assertEqual(len(snr[0]),10)
        self.assertEqual(len(snr[1]),10)
        self.assertFalse(np.any(np.isnan(snr)))

    def test_004_multiple_dimensions_beaviour(self):
        from dnplab.core.base import ABCData
        ABCData.add_proc_attrs=lambda *args,**kwargs: '1'
        d1 = 512
        d2 = 10
        d = np.random.random((d1, d2))
        d[:20,...]=np.linspace(0,10,20)[:,np.newaxis] # noise = 3.03488489333442
        d[150,...]=200 #signal
        snr_r1 = 200/ 3.03488489333442
        ax1 = np.arange(d1)
        ax2=np.arange(d2)
        Ddata = ABCData(d, dims=["f2", "2"], coords=[ax1, ax2])
        f = dnp.processing.signal_to_noise
        s1 = Ddata.shape

        # some input checks, just to check that no errors are thrown:
        snr = f(Ddata, [(100, 300)], [(0, 20)])
        self.assertEqual(len(snr),10)
        self.assertEqual(type(snr[0]),np.float64)
        self.assertTrue(np.isclose(snr_r1,snr[0]))
        snr = f(
            Ddata, [(300, 400)], [(500, 600)], remove_background=(100, 200), deg=3
        )  # works with degree
        self.assertEqual(len(snr),10)
        self.assertEqual(type(snr[0]),np.float64)
        snr = f(
            Ddata, [(300, 400)], [(500, 600)], remove_background=(100, 200)
        )  # works without degree
        self.assertEqual(len(snr),10)
        self.assertEqual(type(snr[0]),np.float64)
        snr = f(Ddata, [(300, 400)], [(500, 600)], remove_background=[(100, 200)])
        self.assertEqual(len(snr),10)
        self.assertEqual(type(snr[0]),np.float64)
        s2 = Ddata.shape
        self.assertEqual(s1, s2)
        # no nan values are returned
        self.assertFalse(np.any(np.isnan(snr)))

    def test_005_multiple_dimensions_and_multiple_regions(self):
        from dnplab.core.base import ABCData
        ABCData.add_proc_attrs=lambda *args,**kwargs: '1'
        d1 = 512
        d2 = 5
        d3=10
        d = np.random.random((d1, d2,d3))
        d[:20,...]=np.linspace(0,10,20)[:,np.newaxis,np.newaxis] # noise = 3.03488489333442
        d[150,...]=200 #signal
        snr_r1 = 200/ 3.03488489333442
        ax1 = np.arange(d1)
        ax2=np.arange(d2)
        ax3=np.arange(d3)
        Ddata = ABCData(d, dims=["f2", "2","3"], coords=[ax1, ax2,ax3])
        f = dnp.processing.signal_to_noise

         # works with degree & 2 regions & noise region specified
        snr = f(
            Ddata, [(100,200),(300, 400)], [(0, 20)], remove_background=(0, 95), deg=1
        )
        self.assertEqual(len(snr),2)
        self.assertEqual(len(snr[0]),d2*d3)
        self.assertEqual(len(snr[1]),d3*d2)
        self.assertFalse(np.any(np.isnan(snr)))
        self.assertTrue( np.all(np.isclose(k,snr_r1) for k in snr[0]) )

         # works with degree & 2 regions & noise region not specified
        snr = f(
            Ddata, [(100,200),(300, 400)], remove_background=(0, 95), deg=1
        )
        self.assertEqual(len(snr),2)
        self.assertEqual(len(snr[0]),d2*d3)
        self.assertEqual(len(snr[1]),d3*d2)
        self.assertFalse(np.any(np.isnan(snr)))


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
