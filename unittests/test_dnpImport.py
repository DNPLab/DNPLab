import unittest
import dnplab.dnpImport as wrapper
import dnplab.dnpIO.topspin as topspin
import dnplab.dnpIO.prospa as prospa
import dnplab.dnpIO.vnmrj as vnmrj
import os
from numpy.testing import assert_array_equal


class import_topspin_tester(unittest.TestCase):
    def setUp(self):
        self.testdata = os.path.join(".", "data", "topspin")

    def test_dir_data_type(self):
        self.assertEqual(
            "fid", topspin.dir_data_type(os.path.join(self.testdata, str(1)))
        )
        self.assertEqual(
            "serPhaseCycle", topspin.dir_data_type(os.path.join(self.testdata, str(5)))
        )
        self.assertEqual(
            "ser", topspin.dir_data_type(os.path.join(self.testdata, str(28)))
        )
        self.assertEqual(
            "", topspin.dir_data_type(os.path.join(self.testdata, str(33)))
        )

    def test_import_topspin_exp1_is_fid(self):
        data = wrapper.load(os.path.join(self.testdata, str(1)), data_type="topspin")
        self.assertEqual(data.dims[0], "t2")
        self.assertEqual(data.values.size, 8147)
        self.assertAlmostEqual(data.values.min(), -5 - 4.168734491315137j)
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)

    def test_import_topspin_exp5_is_2d_phcyc(self):
        data = wrapper.load(os.path.join(self.testdata, str(5)), data_type="topspin")
        self.assertEqual(data.values.shape[0], 11912)
        self.assertEqual(data.dims, ["t2"])
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)
        self.assertAlmostEqual(data.values[365], -0.182861328125 - 0.71875j)

    def test_import_topspin_exp28_is_2d(self):
        data = wrapper.load(os.path.join(self.testdata, str(28)), data_type="topspin")
        self.assertEqual(data.values.shape, (7922, 8))
        self.assertEqual(data.dims, ["t2", "t1"])
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)
        self.assertAlmostEqual(data.values[365, 6], -0.110595703125 + 0.47705078125j)

    def test_import_topspin_jcamp_dx(self):
        attrs = topspin.topspin_jcamp_dx(os.path.join(self.testdata, "1", "acqus"))
        self.assertEqual(attrs["DIGTYP"], 9)
        self.assertAlmostEqual(attrs["O1"], 1413.27)
        assert_array_equal(attrs["XGAIN"], [0, 0, 0, 0])
        assert_array_equal(attrs["TPOAL"], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])


class prospa_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data = os.path.join(".", "data", "prospa", "toluene_10mM_Tempone")

    def test_import_prospa_exp_is_1d(self):
        datas = [
            wrapper.load(
                os.path.join(self.test_data, "%i" % expNum, "data.csv"),
                data_type="prospa",
            )
            for expNum in [1, 21, 42]
        ]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (16384,))
            self.assertEqual(data.dims, ["t2"])
            self.assertAlmostEqual(data.attrs["nmr_frequency"], 14244500.0)
        self.assertAlmostEqual(datas[0].values[365], -0.217937 + 0.24907j)
        self.assertAlmostEqual(datas[1].values[365], 0.0400292 - 0.0756107j)
        self.assertAlmostEqual(datas[2].values[365], 1.09858 - 2.57966j)


class vnmrj_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data2Ds = [
            os.path.join(".", "data", "vnmrj", s)
            for s in ["10mM_tempol_in_water_array.fid"]
        ]
        self.test_data1Ds = [
            os.path.join(".", "data", "vnmrj", s)
            for s in [
                "10mM_tempol_in_water_mw_40dBm.fid",
                "10mM_tempol_in_water_mw_off.fid",
            ]
        ]

    def test_import_vnmrj_1d(self):
        datas = [
            wrapper.load(path=path, data_type="vnmrj") for path in self.test_data1Ds
        ]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (131072,))
            self.assertEqual(
                data.dims,
                [
                    "t2",
                ],
            )
            self.assertAlmostEqual(data.attrs["nmr_frequency"], 14244283.4231)
        self.assertAlmostEqual(datas[0].values[365], (-20378767 - 2734659j))
        self.assertAlmostEqual(datas[1].values[365], (-950662 + 138458j))

    def test_import_vnmrj_2d(self):
        datas = [
            wrapper.load(path=path, data_type="vnmrj") for path in self.test_data2Ds
        ]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (131072, 5))
            self.assertEqual(data.dims, ["t2", "t1"])
            self.assertAlmostEqual(data.attrs["nmr_frequency"], 14244283.4231)
        self.assertAlmostEqual(datas[0].values[365, 3], (-1263136 + 1063328.5j))


class specman_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data_2D = os.path.join(".", "data", "specman", "test_specman2D.exp")
        self.test_data_4D = os.path.join(".", "data", "specman", "test_specman4D.d01")

    def test_import_specman_2D(self):
        data = wrapper.load(self.test_data_2D, data_type="specman")
        self.assertEqual(data.dims, ["t", "x"])
        self.assertEqual(data.values.shape, (1500, 80))

    def test_import_specman_4D(self):
        data = wrapper.load(self.test_data_4D, data_type="specman")
        self.assertEqual(data.dims, ["t", "x", "y", "z"])
        self.assertEqual(data.values.shape, (1500, 40, 5, 3))


class delta_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data_1D = os.path.join(".", "data", "delta", "50percCHCL3.jdf")
        self.test_data_2D = os.path.join(".", "data", "delta", "lineshape_drift.jdf")

    def test_import_delta_1D(self):
        data = wrapper.load(self.test_data_1D, data_type="delta")
        self.assertEqual(data.dims, ["t2"])
        self.assertEqual(data.values.shape, (16384,))

    def test_import_delta_2D(self):
        data = wrapper.load(self.test_data_2D, data_type="delta")
        self.assertEqual(data.dims, ["t2", "t1"])
        self.assertEqual(data.values.shape, (1024, 256))


if __name__ == "__main__":
    unittest.main()
    pass
