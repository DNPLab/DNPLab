import unittest
import dnplab as dnp
from dnplab import dnpdata
from numpy.testing import assert_array_equal
import dnplab.dnpImport as wrapper
import dnplab.dnpNMR as nmr
import numpy as np
import os
import copy


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

    def test_signal_to_noise(self):

        ws = copy.deepcopy(self.ws)

        nmr.remove_offset(ws)
        nmr.window(ws, linewidth=15)
        nmr.fourier_transform(ws, zero_fill_factor=2)
        nmr.autophase(ws, method="arctan")

        dnp.dnpTools.signal_to_noise(
            ws,
            signal_center=0,
            signal_width="full",
            noise_center="default",
            noise_width="default",
        )

        self.assertEqual(len(ws["proc"].attrs["s_n"]), 8)
        self.assertAlmostEqual(ws["proc"].attrs["s_n"][4], 10.91728842, places=6)

        nmr.remove_offset(self.ws_off)
        nmr.window(self.ws_off, linewidth=15)
        nmr.fourier_transform(self.ws_off, zero_fill_factor=2)
        nmr.autophase(self.ws_off, method="search")

        dnp.dnpTools.signal_to_noise(
            self.ws_off,
            signal_center=-10,
            signal_width=100,
            noise_center=200,
            noise_width=50,
        )

        self.assertAlmostEqual(
            self.ws_off["proc"].attrs["s_n"], 3.3986939859030096, places=6
        )

    def test_baseline(self):

        ws = copy.deepcopy(self.ws)
        shape_data = np.shape(ws)

        nmr.remove_offset(ws)
        nmr.window(ws, linewidth=15)
        nmr.fourier_transform(ws, zero_fill_factor=2)
        nmr.autophase(ws, method="search")

        ws = dnp.dnpTools.baseline(ws, type="polynomial", order=1, reference_slice=None)
        self.assertEqual(shape_data, np.shape(ws))
        self.assertAlmostEqual(
            max(ws["proc"].values[:, 2].real), 29.894686521140628, places=4
        )
        self.assertAlmostEqual(
            min(ws["proc"].values[:, 2].real), -921.7417731961074, places=4
        )
        self.assertAlmostEqual(
            len(ws["proc"].attrs["baseline"]), len(ws["proc"].values)
        )

    def test_integrate(self):

        ws = copy.deepcopy(self.ws)

        nmr.remove_offset(ws)
        nmr.window(ws, linewidth=15)
        nmr.fourier_transform(ws, zero_fill_factor=2)
        nmr.autophase(ws, method="arctan")

        ws2 = copy.deepcopy(ws)
        ws3 = copy.deepcopy(ws)
        dnp.dnpTools.integrate(ws, dim="f2", integrate_center=0, integrate_width=50)
        self.assertEqual((8,), np.shape(ws["integrate"].values))
        self.assertAlmostEqual(
            max(ws["integrate"].values.real), 6170.447249940133, places=4
        )
        self.assertAlmostEqual(
            min(ws["integrate"].values.real), -7188.897892203664, places=4
        )

        dnp.dnpTools.integrate(
            ws2, dim="f2", integrate_center=[-10, 0, 10], integrate_width=50
        )
        self.assertEqual((3, 8), np.shape(ws2["integrate"].values))
        self.assertAlmostEqual(
            max(ws2["integrate"].values.real[1]), 6170.447249940133, places=4
        )
        self.assertAlmostEqual(
            min(ws2["integrate"].values.real[1]), -7188.897892203664, places=4
        )

        dnp.dnpTools.integrate(
            ws3, dim="f2", integrate_center=[-10, 0, 10], integrate_width=[10, 50, 10]
        )
        self.assertEqual((3, 8), np.shape(ws2["integrate"].values))
        self.assertAlmostEqual(
            max(ws3["integrate"].values.real[1]), 6170.447249940133, places=4
        )
        self.assertAlmostEqual(
            min(ws3["integrate"].values.real[1]), -7188.897892203664, places=4
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
