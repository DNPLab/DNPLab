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
        self.assertEqual(max(self.ws["proc"].values[:, 1].real), 5.447918664876271)
        self.assertEqual(min(self.ws["proc"].values[:, 1].real), -4.128415225672084)

        nmr.remove_offset(self.ws)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.assertEqual(max(self.ws["proc"].values[:, 1].real), 5.445267326190019)
        self.assertEqual(min(self.ws["proc"].values[:, 1].real), -4.131066564358336)

        self.ws.copy("proc", "temp")
        nmr.window(self.ws, type="hamming")
        self.assertEqual(max(self.ws["proc"].attrs["window"]), 1.0)
        self.assertEqual(min(self.ws["proc"].attrs["window"]), 0.07671999999999995)

        self.ws.copy("temp", "proc")
        self.ws = nmr.window(self.ws, type="hann")
        self.assertEqual(max(self.ws["proc"].attrs["window"]), 1.0)
        self.assertEqual(min(self.ws["proc"].attrs["window"]), 0.0)

        self.ws.copy("temp", "proc")
        nmr.window(self.ws, type="lorentz_gauss", linewidth=5, gauss_linewidth=10)
        self.assertEqual(max(self.ws["proc"].attrs["window"]), 1.1895922020471337)

        self.ws.copy("temp", "proc")
        self.ws = nmr.window(
            self.ws, type="hamming", linewidth=5, gauss_linewidth=10, inverse=True
        )
        self.assertEqual(max(self.ws["proc"].attrs["window"]), 13.03441084462983)
        self.ws.copy("temp", "proc")
        self.ws = nmr.window(self.ws, type="sin2")
        self.assertEqual(max(self.ws["proc"].attrs["window"]), 1.0)
        self.assertEqual(min(self.ws["proc"].attrs["window"]), 3.749399456654644e-33)

        self.ws.copy("temp", "proc")
        nmr.window(self.ws, type="exponential", linewidth=5)
        self.assertEqual(min(self.ws["proc"].attrs["window"]), 0.00035733315645396175)
        self.ws.pop("temp")
        self.assertEqual(shape_data, np.shape(self.ws))
        self.assertEqual(max(self.ws["proc"].values[:, 1].real), 5.390978190372195)
        self.assertEqual(min(self.ws["proc"].values[:, 1].real), -3.9647928998767163)

        nmr.fourier_transform(self.ws, zero_fill_factor=2)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.assertEqual(len(self.ws["proc"].values[:, 4].real), 15844)
        self.assertEqual(max(self.ws["proc"].values[:, 4].real), 61.3994271072369)
        self.assertEqual(min(self.ws["proc"].values[:, 4].real), -76.53300374124703)
        nmr.autophase(
            self.ws,
            method="arctan",
            order="first",
            pivot=len(self.ws["proc"].values[:, 7]) / 2,
            delta=np.pi / 2,
            reference_slice=None,
            force_positive=False,
        )
        self.assertEqual(shape_data, np.shape(self.ws))
        self.assertEqual(max(self.ws["proc"].values[:, 3].real), 235.73708636474333)
        self.assertEqual(min(self.ws["proc"].values[:, 3].real), -63.73460092653687)
        phs0 = self.ws["proc"].attrs["phase_0"]
        self.assertEqual(
            len(self.ws["proc"].attrs["phase_1"]), len(self.ws["proc"].values)
        )
        self.ws = nmr.baseline(self.ws, type="poly", order=1, reference_slice=None)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.assertEqual(max(self.ws["proc"].values[:, 2].real), 371.9768176047327)
        self.assertEqual(min(self.ws["proc"].values[:, 2].real), -71.66041320221538)
        self.assertEqual(
            len(self.ws["proc"].attrs["baseline"]), len(self.ws["proc"].values)
        )
        nmr.integrate(self.ws, dim="t2", integrate_center=0, integrate_width=50)
        self.assertEqual((8,), np.shape(self.ws["proc"].values))
        self.assertEqual(max(self.ws["proc"].values.real), 2135.920732152984)
        self.assertEqual(min(self.ws["proc"].values.real), -1726.9213059614328)


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
