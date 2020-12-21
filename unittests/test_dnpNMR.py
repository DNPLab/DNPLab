import unittest
from dnplab import dnpdata
from numpy.testing import assert_array_equal
import dnplab.dnpImport as wrapper
import dnplab.dnpNMR as nmr
import dnplab as dnp
import numpy as np
import os


class dnpNMR_tester(unittest.TestCase):
    def setUp(self):
        path = os.path.join(".", "data", "topspin", "304")
        self.data = wrapper.load(path, data_type="topspin")
        self.ws = dnp.create_workspace()
        self.ws["raw"] = self.data
        self.ws.copy("raw", "proc")

    def test_basic_nmr_processing(self):
        shape_data = np.shape(self.ws)
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 1].real), 5.447918664876271, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 1].real), -4.128415225672084, places=4
        )

        nmr.remove_offset(self.ws)
        self.assertAlmostEqual(shape_data, np.shape(self.ws))
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 1].real), 5.445267326190019, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 1].real), -4.131066564358336, places=4
        )

        self.ws.copy("proc", "temp")
        nmr.window(self.ws, type="hamming")
        self.assertEqual(max(self.ws["proc"].attrs["window"]), 1.0)
        self.assertAlmostEqual(
            min(self.ws["proc"].attrs["window"]), 0.07671999999999995, places=4
        )

        self.ws.copy("temp", "proc")
        self.ws = nmr.window(self.ws, type="hann")
        self.assertEqual(max(self.ws["proc"].attrs["window"]), 1.0)
        self.assertEqual(min(self.ws["proc"].attrs["window"]), 0.0)

        self.ws.copy("temp", "proc")
        nmr.window(self.ws, type="lorentz_gauss", linewidth=5, gauss_linewidth=10)
        self.assertAlmostEqual(
            max(self.ws["proc"].attrs["window"]), 1.1895922020471337, places=4
        )

        self.ws.copy("temp", "proc")
        self.ws = nmr.window(
            self.ws, type="hamming", linewidth=5, gauss_linewidth=10, inverse=True
        )
        self.assertAlmostEqual(
            max(self.ws["proc"].attrs["window"]), 13.03441084462983, places=4
        )
        self.ws.copy("temp", "proc")
        self.ws = nmr.window(self.ws, type="sin2")
        self.assertEqual(max(self.ws["proc"].attrs["window"]), 1.0)

        self.ws.copy("temp", "proc")
        nmr.window(self.ws, type="exponential", linewidth=5)
        self.assertAlmostEqual(
            min(self.ws["proc"].attrs["window"]), 0.00035733315645396175, places=4
        )
        self.ws.pop("temp")
        self.assertAlmostEqual(shape_data, np.shape(self.ws))
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 1].real), 5.390978190372195, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 1].real), -3.9647928998767163, places=4
        )

        nmr.fourier_transform(self.ws, zero_fill_factor=2)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.assertEqual(len(self.ws["proc"].values[:, 4].real), 15844)
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 4].real), 61.3994271072369, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 4].real), -76.53300374124703, places=4
        )

        self.ws.copy("proc", "temp")
        nmr.autophase(
            self.ws,
            method="search",
            order="first",
            pivot=len(self.ws["proc"].values[:, 7]) / 2,
            delta=np.pi / 2,
            reference_slice=None,
            force_positive=False,
        )
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 3].real), 242.36299886442168, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 3].real), -65.59778997862813, places=4
        )
        phs0 = self.ws["proc"].attrs["phase_0"]
        phs1 = self.ws["proc"].attrs["phase_1"]
        self.assertEqual(len(phs1), len(self.ws["proc"].values))
        self.ws.copy("proc", "keep")
        self.ws.copy("temp", "proc")

        nmr.autophase(
            self.ws,
            method="arctan",
            order="zero",
            reference_slice=1,
            force_positive=True,
        )
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 7].real), 398.66277628048533, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 7].real), -47.32339420356964, places=4
        )
        self.ws.copy("temp", "proc")

        nmr.autophase(
            self.ws,
            method="manual",
            order="first",
            pivot=len(self.ws["proc"].values[:, 7]) / 2,
            delta=np.pi / 2,
            phase=phs1 * (45 * np.pi / 180),
        )
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 3].real), 239.3422920589204, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 3].real), -64.37971134351001, places=4
        )
        self.ws.copy("temp", "proc")

        nmr.autophase(
            self.ws,
            method="manual",
            order="zero",
            phase=phs0 * (45 * np.pi / 180),
        )
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 3].real), 238.27990936049952, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 3].real), -66.28951430959444, places=4
        )
        self.ws.pop("temp")
        self.ws.copy("keep", "proc")
        self.ws.pop("keep")

        self.ws = nmr.baseline(self.ws, type="poly", order=1, reference_slice=None)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.assertAlmostEqual(
            max(self.ws["proc"].values[:, 2].real), 351.9904580259595, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values[:, 2].real), -73.68315636898302, places=4
        )
        self.assertAlmostEqual(
            len(self.ws["proc"].attrs["baseline"]), len(self.ws["proc"].values)
        )
        nmr.integrate(self.ws, dim="f2", integrate_center=0, integrate_width=50)
        self.assertEqual((8,), np.shape(self.ws["proc"].values))
        self.assertAlmostEqual(
            max(self.ws["proc"].values.real), 2031.4439405548092, places=4
        )
        self.assertAlmostEqual(
            min(self.ws["proc"].values.real), -1655.6511804282302, places=4
        )


class dnpNMR_tester_sim(unittest.TestCase):
    def setUp(self):
        p1 = np.array([0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0])
        p2 = np.array([0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0, 0, 0])
        p3 = np.array([0, 0, 0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0])
        self.data = dnpdata(
            np.array([p1, p2, p3]).T,
            [np.arange(0, len(p1)), np.arange(0, 3)],
            ["x", "t2"],
        )
        self.ws = dnp.create_workspace()
        self.ws["raw"] = self.data
        self.ws.copy("raw", "proc")

    def test_align(self):
        nmr.align(self.ws, dim="x")
        assert_array_equal(
            self.ws["proc"].values,
            np.array(
                [
                    [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0],
                    [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0],
                    [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0],
                ]
            ).T,
        )
