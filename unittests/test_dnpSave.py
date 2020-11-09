import unittest
import dnplab.dnpImport as wrapper
import dnplab.dnpSave as saver
from dnplab import create_workspace
import os
import time
from numpy.testing import assert_array_equal


class save_h5_tester(unittest.TestCase):
    def setUp(self):
        self.testdata_hyd = wrapper.load(
            os.path.join(
                ".",
                "data",
                "topspin",
                "hydrationGUI_Results",
                "hydration_parameters.h5",
            )
        )
        self.testdata_wsob = wrapper.load(
            os.path.join(
                ".",
                "data",
                "topspin",
                "304",
            )
        )

    def test_h5_save_hyd(self):

        saver.save(
            self.testdata_hyd,
            os.path.join(".", "unittests", "test_hyd.h5"),
            overwrite=True,
        )

        saved_data = wrapper.load(os.path.join(".", "unittests", "test_hyd.h5"))

        self.assertEqual(
            self.testdata_hyd["hydration_inputs"].keys(),
            saved_data["hydration_inputs"].keys(),
        )

        self.assertEqual(
            self.testdata_hyd["hydration_results"]["ksigma"],
            saved_data["hydration_results"]["ksigma"],
        )

        os.remove(os.path.join(".", "unittests", "test_hyd.h5"))

    def test_h5_save_ws(self):

        saver.save(
            self.testdata_wsob,
            os.path.join(".", "unittests", "test_ob.h5"),
            overwrite=True,
        )

        saved_data_ob = wrapper.load(os.path.join(".", "unittests", "test_ob.h5"))

        self.assertEqual(
            max(self.testdata_wsob.values[0]),
            max(saved_data_ob["data"].values[0]),
        )

        self.assertEqual(self.testdata_wsob.values.shape, saved_data_ob["data"].shape)

        self.assertEqual(self.testdata_wsob.dims, saved_data_ob["data"].dims)

        os.remove(os.path.join(".", "unittests", "test_ob.h5"))

        ws_ws = create_workspace()
        ws_ws.add("raw", self.testdata_wsob)
        ws_ws.copy("raw", "proc")

        saver.save(ws_ws, os.path.join(".", "unittests", "test_ws.h5"), overwrite=True)

        saved_data_ws = wrapper.load(os.path.join(".", "unittests", "test_ws.h5"))

        self.assertEqual(
            saved_data_ws["raw"].attrs.keys(),
            saved_data_ws["proc"].attrs.keys(),
        )

        self.assertEqual(
            ws_ws["proc"].attrs.keys(),
            saved_data_ws["proc"].attrs.keys(),
        )

        self.assertEqual(ws_ws["proc"].shape, saved_data_ws["proc"].shape)

        self.assertEqual(ws_ws["proc"].dims, saved_data_ws["proc"].dims)

        os.remove(os.path.join(".", "unittests", "test_ws.h5"))


if __name__ == "__main__":
    unittest.main()
    pass
