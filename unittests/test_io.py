import unittest
import dnplab as dnp
import os
from numpy.testing import assert_array_equal
import logging

logger = logging.getLogger(__name__)


class import_topspin_tester(unittest.TestCase):
    def setUp(self):
        self.testdata = os.path.join(".", "data", "topspin")

    def test_import_topspin_exp1_is_fid(self):
        data = dnp.load(os.path.join(self.testdata, str(1)), data_format="topspin")
        self.assertEqual(data.dims[0], "t2")
        self.assertEqual(data.values.size, 8192)
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)

    def test_import_topspin_exp5_is_2d_phcyc(self):
        data = dnp.load(os.path.join(self.testdata, str(5)), data_format="topspin")
        self.assertEqual(data.values.shape[0], 11973)
        self.assertEqual(data.dims, ["t2", "t1"])
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)

    def test_import_topspin_exp28_is_2d(self):
        data = dnp.load(os.path.join(self.testdata, str(28)), data_format="topspin")
        self.assertEqual(data.values.shape, (7983, 8))
        self.assertEqual(data.dims, ["t2", "t1"])
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)

    def test_import_topspin_jcamp_dx(self):
        attrs = dnp.io.topspin.load_topspin_jcamp_dx(
            os.path.join(self.testdata, "1", "acqus")
        )
        self.assertEqual(attrs["DIGTYP"], 9)
        self.assertAlmostEqual(attrs["O1"], 1413.27)
        assert_array_equal(attrs["XGAIN"], [0, 0, 0, 0])
        assert_array_equal(attrs["TPOAL"], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])


class prospa_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data = os.path.join(".", "data", "prospa", "toluene_10mM_Tempone")

    def test_import_prospa_exp_is_1d(self):
        datas = [
            dnp.load(
                os.path.join(self.test_data, "%i" % expNum, "data.csv"),
                data_format="prospa",
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
        datas = [dnp.load(path=path, data_format="vnmrj") for path in self.test_data1Ds]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (131072,))
            self.assertEqual(
                data.dims,
                [
                    "t2",
                ],
            )
            self.assertAlmostEqual(data.attrs["nmr_frequency"], 14244283.4231)
        self.assertAlmostEqual(datas[0].values[365], (-20378767 + 2734659j))
        self.assertAlmostEqual(datas[1].values[365], (-950662 - 138458j))

    def test_import_vnmrj_2d(self):
        datas = [dnp.load(path=path, data_format="vnmrj") for path in self.test_data2Ds]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (131072, 5))
            self.assertEqual(data.dims, ["t2", "t1"])
            self.assertAlmostEqual(data.attrs["nmr_frequency"], 14244283.4231)
        self.assertAlmostEqual(datas[0].values[365, 3], (-1263136 - 1063328.5j))


class specman_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data_2D = os.path.join(".", "data", "specman", "test_specman2D.exp")
        self.test_data_4D = os.path.join(".", "data", "specman", "test_specman4D.d01")
        self.test_data_field_monitor = os.path.join(".", "data", "specman", "test_specman_field_monitor.exp")

    def test_import_specman_2D(self):
        data = dnp.load(self.test_data_2D, data_format="specman")
        self.assertEqual(data.dims, ["x0", "x1", "x2"])
        self.assertEqual(data.values.shape, (4500, 252, 2))

    def test_import_specman_4D(self):
        data = dnp.load(self.test_data_4D, data_format="specman")
        self.assertEqual(data.dims, ["x0", "x1", "x2", "x3", "x4"])
        self.assertEqual(data.values.shape, (1500, 40, 5, 3, 2))

    def test_import_specman_2D_with_autodetect(self):
        data = dnp.load(self.test_data_2D, data_format="specman", autodetect_dims = True, autodetect_coords = True)
        self.assertEqual(data.dims, ["t2", "t", "x"])
        self.assertEqual(data.values.shape, (4500, 252, 2))
        
    def test_import_specman_4D_with_autodetect(self):
        data = dnp.load(self.test_data_4D, data_format="specman", autodetect_dims = True, autodetect_coords = True)
        self.assertEqual(data.dims, ["t2", "Fr_pump", "offset1", "tsquare", "x"])
        self.assertEqual(data.values.shape, (1500, 40, 5, 3, 2))

    def test_import_specman_field_monitor(self):
        data = dnp.load(self.test_data_field_monitor, data_format="specman", autodetect_dims = True, autodetect_coords = True)
        self.assertEqual(data.dims, ["tau", "Field", "x"])
        self.assertEqual(data.values.shape, (101, 101, 2))
        self.assertEqual(data.coords['tau'][0], 3.0000000000000004e-07)

class bes3t_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data_HYSCORE = os.path.join(".", "data", "bes3t", "HYSCORE.DSC")
        self.test_data_DEER = os.path.join(".", "data", "bes3t", "DEER.DSC")
        self.test_data_ESE = os.path.join(".", "data", "bes3t", "2D_ESE.DTA")
        self.test_data_1D = os.path.join(".", "data", "bes3t", "1D_CW.DTA")
        self.test_data_2D = os.path.join(".", "data", "bes3t", "2D_CW.YGF")

    def test_import_bes3t_HYSCORE(self):
        data = dnp.load(self.test_data_HYSCORE, data_format="xepr")
        self.assertEqual(data.dims, ["t2", "t1"])
        self.assertEqual(data.values.shape, (175, 175))
        self.assertEqual(max(data.coords["t2"]), 3520.0)
        self.assertEqual(max(data.coords["t1"]), 3520.0)

    def test_import_bes3t_DEER(self):
        data = dnp.load(self.test_data_DEER, data_format="xepr")
        self.assertEqual(data.dims, ["t2"])
        self.assertEqual(data.values.shape, (504,))
        self.assertEqual(data.attrs["frequency"], 33.85)

    def test_import_bes3t_ESE(self):
        data = dnp.load(self.test_data_ESE, data_format="xepr")
        self.assertEqual(data.dims, ["t2", "t1"])
        self.assertEqual(data.values.shape, (512, 50))
        self.assertEqual(data.attrs["frequency"], 9.296)

    def test_import_bes3t_1D(self):
        data = dnp.load(self.test_data_1D, data_format="xenon")
        self.assertEqual(data.dims, ["B0"])
        self.assertEqual(data.values.shape, (2250,))
        self.assertEqual(data.attrs["frequency"], 9.804448)

    def test_import_bes3t_2D(self):
        data = dnp.load(self.test_data_2D, data_format="xenon")
        self.assertEqual(data.dims, ["B0", "t1"])
        self.assertEqual(data.values.shape, (1600, 100))
        self.assertEqual(data.attrs["frequency"], 9.627213)


class winepr_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data_ESP = os.path.join(".", "data", "parspc", "ExampleESP.par")
        self.test_data_1D = os.path.join(".", "data", "parspc", "Example1D.spc")
        self.test_data_2D = os.path.join(".", "data", "parspc", "Example2D.spc")

    def test_import_winepr_ESP(self):
        data = dnp.load(self.test_data_ESP, data_format="esp")
        self.assertEqual(data.dims, ["t2"])
        self.assertEqual(data.values.shape, (1024,))
        self.assertEqual(data.attrs["conversion_time"], 81.92)

    def test_import_winepr_1D(self):
        data = dnp.load(self.test_data_1D, data_format="winepr")
        self.assertEqual(data.dims, ["B0"])
        self.assertEqual(data.values.shape, (512,))
        self.assertEqual(data.attrs["temperature"], 294.2)

    def test_import_winepr_2D(self):
        data = dnp.load(self.test_data_2D, data_format="winepr")
        self.assertEqual(data.dims, ["B0", "t1"])
        self.assertEqual(data.values.shape, (1024, 15))
        self.assertEqual(data.attrs["frequency"], 9.79)


class delta_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data_1D = os.path.join(".", "data", "delta", "50percCHCL3.jdf")
        self.test_data_2D = os.path.join(".", "data", "delta", "lineshape_drift.jdf")

    def test_import_delta_1D(self):
        data = dnp.load(self.test_data_1D, data_format="delta")
        self.assertEqual(data.dims, ["t2"])
        self.assertEqual(data.values.shape, (16384,))
        self.assertEqual(max(data.coords["t2"]), 0.262128)

    def test_import_delta_2D(self):
        data = dnp.load(self.test_data_2D, data_format="delta")
        self.assertEqual(data.dims, ["t2", "t1"])
        self.assertEqual(data.values.shape, (8192, 256))
        self.assertEqual(max(data.coords["t2"]), 0.5451929600000001)
        self.assertEqual(max(data.coords["t1"]), 11.953125)


class csv_import_tester(unittest.TestCase):
    def setUp(self):
        self.testdata = os.path.join(".", "data", "csv")

    def test_import_csv_arrLNA_fid(self):
        import pathlib

        p = pathlib.Path(self.testdata)
        data = dnp.io.load_csv.load_csv(
            p.joinpath("csv_example.csv"),
            skiprows=1,
            maxrows=115,
            tcol=0,
            real=1,
            imag=3,
        )
        self.assertEqual(data.dims[0], "t2")
        self.assertEqual(data.values[1], 5e3 + 1j * 25000)
        self.assertEqual(data.coords[0][1], 20)
        self.assertEqual(data.values.size, 115)

    def test_remove_data_csv_arrLNA_fid(self):
        import pathlib

        p = pathlib.Path(self.testdata)
        data = dnp.io.load_csv.load_csv(
            p.joinpath("csv_example.csv"),
            skiprows=1,
            maxrows=115,
            tcol=None,
            real=1,
            imag=3,
        )
        self.assertEqual(data.dims[0], "t2")
        self.assertEqual(data.values[1], 5e3 + 1j * 25000)
        self.assertEqual(data.coords[0][100], 100)
        self.assertEqual(data.values.size, 115)

    def test_set_imag_to_zero(self):
        import pathlib

        p = pathlib.Path(self.testdata)
        data = dnp.io.load_csv.load_csv(
            p.joinpath("csv_example.csv"),
            skiprows=1,
            maxrows=115,
            tcol=None,
            real=1,
            imag=None,
        )
        self.assertEqual(data.dims[0], "t2")
        self.assertEqual(data.values[1], 5e3)
        self.assertEqual(data.coords[0][100], 100)
        self.assertEqual(data.values.size, 115)


class dnplab_configparse_tester(unittest.TestCase):
    def test_000_escape_split(self):
        # config
        import sys
        import configparser
        from pathlib import Path

        p = Path(__file__).parent.joinpath("dnplab")
        sys.path.insert(0, p)
        p = str(p)
        from dnplab import config as dnpconfig

        cfg = configparser.ConfigParser(
            converters={
                "list": lambda x: list(x.strip("[").strip("]").split(",")),
                "args_kwargs": dnpconfig.config._kwarg_converter,
            }
        )

        string1 = "Contact Time t$_c$ [s]"
        string2 = r"abc=1,def\=2,ghi=3"

        cfg_file = str(Path(__file__).parent.joinpath("data_testconfig.cfg"))
        cfg.read(cfg_file)

        args2, kwargs2 = cfg.getargs_kwargs("UNITTEST_EXAMPLE", "test1")
        logger.info("{0}\n{1}".format(args2, kwargs2))
        self.assertEqual(len(args2), 1)
        self.assertEqual(len(kwargs2), 2)
        self.assertEqual(kwargs2["ghi"], "3")
        self.assertEqual(args2[0], "def=2")

        args1, kwargs1 = cfg.getargs_kwargs("UNITTEST_EXAMPLE", "test0")
        logger.info("{0}\n{1}".format(args1, kwargs1))
        self.assertEqual(len(args1), 1)
        self.assertEqual(len(kwargs1), 0)
        self.assertEqual(args1[0], string1)


if __name__ == "__main__":
    unittest.main()
    pass
