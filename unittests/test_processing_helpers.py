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

        self.rnd_data=self.data.copy()

    def test_999_quicktests(self):

        # create_complex
        data_r=np.ones(100)
        data_c=np.ones(100)*2
        #
        x = np.r_[0:100]
        y = np.array([x**2.0,x**3.0]).T

        #leads to warning
        rnd_data=dnp.DNPData(y, ["t2",'bla'], [x])
        dnp.processing.create_complex(rnd_data,data_r,data_c)

        #needs integrals:
        x = np.r_[0:100]
        y = np.array([x**2.0,x**3.0]).T
        data = dnp.DNPData(y, ["t2","Power"], [x,np.array([0,1])])

        ft_data=dnp.fourier_transform(data)
        integrals = dnp.integrate(ft_data)
        dnp.processing.calculate_enhancement(integrals)

        #makes data not consistent!
        # smooth
        dnp.processing.smooth(data,window_length=3,polyorder=2)

        # left_shift
        dnp.processing.left_shift(self.data,shift_points=5)

        # normalize
        dnp.processing.normalize(self.data)

        #reference
        dnp.processing.reference(self.data,dim='t2')

        #pseudo_modulation
        dnp.processing.pseudo_modulation(self.data,0.1,dim='t2')


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
        self.assertTrue(not np.isnan(snr._values))

        snr = f(
            data,
            (300, 400),
            (500, 600),
        )
        # dnpdata object as output
        self.assertTrue(type(snr), type(self.data))

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
        self.assertEqual(snr.shape, (1,))

        # with defaults
        snr = f(data)
        # with slices
        snr = f(data, slice(0, None), remove_background=[(100, 200)])

        # with more than one signal region:
        snr = f(data, [slice(0, None), (100, 300)])
        self.assertEqual(snr.shape[0], 2)

        self.assertEqual(snr.shape, (2,))

        # multiple noise regions
        snr2 = f(data, (0, 1000), remove_background=[(100, 200)])
        snr = f(
            data,
            slice(0, None),
            [slice(0, 100), slice(500, 600)],
            remove_background=[(100, 200)],
        )

    def test_002_SNR_on_higherDimensionalData(self):
        coords3 = [np.arange(0, 100), np.arange(0, 20), np.arange(0, 40)]
        data3 = np.random.random((100, 20, 40))
        DNPObj3 = dnp.DNPData(data3, ["t2", "t3", "t4"], coords3)
        f = dnp.processing.signal_to_noise

        # single snr region
        snr0 = f(DNPObj3, (10, 20), [(80, 90)], dim="t2")
        logger.info("snr0 (single regions) value shape is {0}".format(snr0.shape))
        self.assertEqual(len(snr0.shape), 3)
        self.assertEqual(snr0.shape, (1, 20, 40))

        snr = f(DNPObj3, [(10, 20), (30, 40), (50, 60)], [(80, 90)], dim="t2")
        self.assertEqual(snr.shape[0], 3)
        self.assertEqual(len(snr.shape), 3)
        self.assertEqual(snr.shape[1], 20)
        self.assertEqual(snr.shape[2], 40)

    def test_003_correct_snr_attribution(self):
        # create artificial testdata
        data = np.empty((100, 5, 8))
        for u in range(100):
            for k in range(5):
                for l in range(8):
                    # idea: [0,1,2,3,4] + [l*10+k+u*100 if x==50 else 0 for x in range(95)] along u
                    if u < 5:
                        data[u, k, l] = u
                    elif u == 50:
                        data[u, k, l] = l * 10 + k + u * 100
                    else:
                        data[u, k, l] = 0
        dims = ["f2", "a1", "a2"]
        coords = [np.arange(100), np.arange(5), np.arange(8)]
        DNPObj = dnp.DNPData(data, dims, coords)

        snr = dnp.processing.signal_to_noise(DNPObj, (45, 55), (0, 5), dim="f2")

        noise = np.std(np.arange(5))
        signal_10_2_5 = 5 * 10 + 2 + 100 * 10

        self.assertTrue(
            snr["signal_region", 0, "a1", 2, "a2", 5], signal_10_2_5 / noise
        )
