import os
import unittest

import numpy as np
from numpy.testing import assert_array_equal

import dnplab as dnp
import dnplab.dnpImport as wrapper
import dnplab.dnpNMR as nmr
from dnplab import dnpdata


class dnpNMR_tester(unittest.TestCase):
    def setUp(self):
        path = os.path.join(".", "data", "vnmrj", "10mM_tempol_in_water_array.fid")
        self.data = wrapper.load(path, data_type="vnmrj")
        self.ws = dnp.create_workspace()
        self.ws["raw"] = self.data
        self.ws.copy("raw", "proc")

    def test_basic_nmr_processing(self):
        shape_data = np.shape(self.ws)
        self.ws = nmr.remove_offset(self.ws)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.ws = nmr.window(self.ws)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.ws = nmr.fourier_transform(self.ws)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.ws = nmr.autophase(self.ws, order="first")
        self.assertEqual(shape_data, np.shape(self.ws))
        phs0 = self.ws["proc"].attrs["phase_0"]
        self.assertEqual(
            len(self.ws["proc"].attrs["phase_1"]), len(self.ws["proc"].values)
        )
        self.ws = nmr.baseline(self.ws)
        self.assertEqual(shape_data, np.shape(self.ws))
        self.assertEqual(
            len(self.ws["proc"].attrs["baseline"]), len(self.ws["proc"].values)
        )

    def test_integrate(self):
        values = self.ws["proc"].values
        # Doing so must not change the workspace at all
        data = self.ws["proc"]
        nmr.integrate(data)
        np.testing.assert_array_equal(self.ws["proc"].values, values)


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
