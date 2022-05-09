import unittest
import dnplab as dnp
import os


class load_wrapper_tester(unittest.TestCase):
    def setUp(self):
        self.topspin_dir = os.path.join(".", "data", "topspin")
        self.prospa_dir = os.path.join(".", "data", "prospa", "toluene_10mM_Tempone")
        self.vnmrj_dir = os.path.join(
            ".", "data", "vnmrj", "10mM_tempol_in_water_array.fid"
        )

        self.list_prospa_dir = [
            "./data/prospa/toluene_10mM_Tempone/1/data.1d",
            "./data/prospa/toluene_10mM_Tempone/2/data.1d",
            "./data/prospa/toluene_10mM_Tempone/3/data.1d",
            "./data/prospa/toluene_10mM_Tempone/4/data.1d",
            "./data/prospa/toluene_10mM_Tempone/5/data.1d",
            "./data/prospa/toluene_10mM_Tempone/6/data.1d",
            "./data/prospa/toluene_10mM_Tempone/7/data.1d",
            "./data/prospa/toluene_10mM_Tempone/8/data.1d",
            "./data/prospa/toluene_10mM_Tempone/9/data.1d",
            "./data/prospa/toluene_10mM_Tempone/10/data.1d",
        ]
        self.list_index = range(10)
        self.dim_name = "test_dim"

    def test_topspin(self):
        data = dnp.load(os.path.join(self.topspin_dir, str(1)), data_type="topspin")
        self.assertEqual(data.dims[0], "t2")
        self.assertEqual(data.values.size, 8148)
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)

    def test_prospa(self):
        dnp.load(os.path.join(self.prospa_dir, str(1), "data.csv"), data_type="prospa")

    def test_vnmrj(self):
        dnp.load(self.vnmrj_dir, data_type="vnmrj")

    def test_import_list(self):
        dnp.load(
            self.list_prospa_dir,
            data_type="prospa",
            dim=self.dim_name,
            coord=self.list_index,
        )


if __name__ == "__main__":
    unittest.main()
    pass
