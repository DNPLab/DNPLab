import unittest
import dnplab.dnpImport as wrapper
import dnplab.dnpSave as saver
import os
import time
from numpy.testing import assert_array_equal


class save_h5_tester(unittest.TestCase):
    def setUp(self):
        self.testdata = wrapper.load(
            os.path.join(
                ".",
                "data",
                "topspin",
                "hydrationGUI_Results",
                "hydration_parameters.h5",
            )
        )

    def test_h5_save(self):

        saver.save(
            self.testdata, os.path.join(".", "unittests", "test.h5"), overwrite=True
        )

        saved_data = wrapper.load(os.path.join(".", "unittests", "test.h5"))

        self.assertEqual(
            self.testdata["hydration_inputs"].keys(),
            saved_data["hydration_inputs"].keys(),
        )

        self.assertEqual(
            self.testdata["hydration_results"]["ksigma"],
            saved_data["hydration_results"]["ksigma"],
        )

        os.remove(os.path.join(".", "unittests", "test.h5"))


if __name__ == "__main__":
    unittest.main()
    pass
