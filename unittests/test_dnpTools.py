import unittest
import dnplab as dnp
from dnplab import dnpdata
from numpy.testing import assert_array_equal
import dnplab.dnpImport as wrapper
import dnplab.dnpNMR as nmr
import numpy as np
import os


class dnpTools_tester(unittest.TestCase):
    def setUp(self):
        path = os.path.join(".", "data", "topspin", "32")
        data = wrapper.load(path, data_type="topspin")
        self.ws = dnp.create_workspace()
        self.ws.add("raw", data)
        self.ws.copy("raw", "proc")

        path = os.path.join(".", "data", "topspin", "5")
        data = wrapper.load(path, data_type="topspin")
        self.ws_off = dnp.create_workspace()
        self.ws_off.add("raw", data)
        self.ws_off.copy("raw", "proc")

    def test_calc_enhancement(self):

        nmr.remove_offset(self.ws)
        nmr.window(self.ws, linewidth=15)
        nmr.fourier_transform(self.ws, zero_fill_factor=2)
        nmr.autophase(self.ws, method="arctan")

        dnp.dnpTools.calculate_enhancement(
            self.ws,
            off_spectrum=1,
            on_spectra="all",
            integrate_center=0,
            integrate_width="full",
            method="integrate",
            dim="f2",
        )

        self.assertAlmostEqual(self.ws["enhancement"].values[0], 0.93468899, places=6)
        self.assertAlmostEqual(self.ws["enhancement"].values[-1], -1.10402673, places=6)

        nmr.remove_offset(self.ws_off)
        nmr.window(self.ws_off, linewidth=15)
        nmr.fourier_transform(self.ws_off, zero_fill_factor=2)
        nmr.autophase(self.ws_off, method="arctan")

        dnp.dnpTools.calculate_enhancement(
            self.ws_off,
            off_spectrum=self.ws_off["proc"],
            on_spectra=self.ws["proc"],
            integrate_center="max",
            integrate_width="full",
            method="amplitude",
            dim="f2",
        )

        self.assertAlmostEqual(
            self.ws_off["enhancement"].values[0], 0.02584709, places=6
        )
        self.assertAlmostEqual(
            self.ws_off["enhancement"].values[-1], -0.02576615, places=6
        )

    def test_mr_properties(self):

        info_1H = dnp.dnpTools.mr_properties("1H")
        self.assertEqual(info_1H, 26.7522128)

        info_1H = dnp.dnpTools.mr_properties("1H", 0.35)
        self.assertEqual(info_1H, 14902114.17018196)

        info_2H = dnp.dnpTools.mr_properties("2H", "qmom")
        self.assertEqual(info_2H, 0.286)

        info_6Li = dnp.dnpTools.mr_properties("6Li", "natAbundance")
        self.assertEqual(info_6Li, 7.59)

        info_6Li = dnp.dnpTools.mr_properties("6Li", "relSensitivity")
        self.assertEqual(info_6Li, 0.000645)
