import unittest
from dnplab import import_data


class import_data_tester(unittest.TestCase):
    def setUp(self):
        self.topspin_dir = "./data/topspin/"
        self.prospa_dir = "./data/prospa/toluene_10mM_Tempone"
        self.vnmrj_dir = "./data/vnmrj/10mM_tempol_in_water_array.fid"

    def test_topspin(self):
        data = import_data(self.topspin_dir, "topspin", 1)  # last arg is expNum
        self.assertEqual(data.dims[0], "t2")
        self.assertEqual(data.values.size, 8147)
        self.assertAlmostEqual(data.values.min(), -5 - 4.168734491315137j)
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)

    def test_prospa(self):
        import_data(self.prospa_dir + "/1/data.csv", "prospa")

    def test_vnmrj(self):
        import_data(self.vnmrj_dir, "vnmrj")


if __name__ == "__main__":
    unittest.main()
    pass
