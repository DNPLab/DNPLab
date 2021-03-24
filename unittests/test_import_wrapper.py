import unittest
from dnplab import load
import os


class load_wrapper_tester(unittest.TestCase):
    def setUp(self):
        self.topspin_dir = os.path.join(".", "data", "topspin")
        self.prospa_dir = os.path.join(".", "data", "prospa", "toluene_10mM_Tempone")
        self.vnmrj_dir = os.path.join(
            ".", "data", "vnmrj", "10mM_tempol_in_water_array.fid"
        )

    def test_topspin(self):
        data = load(os.path.join(self.topspin_dir, str(1)), data_type="topspin")
        self.assertEqual(data.dims[0], "t2")
        self.assertEqual(data.values.size, 8147)
        self.assertAlmostEqual(data.values.min(), -5 - 4.168734491315137j)
        self.assertAlmostEqual(data.attrs["nmr_frequency"], 14831413.270000001)

    def test_prospa(self):
        load(os.path.join(self.prospa_dir, str(1), "data.csv"), data_type="prospa")

    def test_vnmrj(self):
        load(self.vnmrj_dir, data_type="vnmrj")


if __name__ == "__main__":
    unittest.main()
    pass
