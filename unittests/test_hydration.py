import unittest
import numpy as np
import dnplab as dnp


TESTSET = {
    "T1": np.array(
        [
            2.020153734009,
            2.276836030132750,
            2.3708172489377400,
            2.4428968088189100,
            2.5709096032675700,
        ]
    ),
    "T1_power": np.array(
        [
            0.000589495934876689,
            0.024242327290569100,
            0.054429505156431400,
            0.0862844940360515,
            0.11617812912435900,
        ]
    ),
    "E": np.array(
        [
            0.57794113752189,
            -0.4688718613022250,
            -0.5464528159680670,
            -1.0725090541762200,
            -1.4141203961920700,
            -1.695789643686440,
            -1.771840068080760,
            -1.8420812985152700,
            -1.97571340381877,
            -2.091405209753480,
            -2.1860546327712800,
            -2.280712535872610,
            -2.4709892163826400,
            -2.5184316153191200,
            -2.556110148443770,
            -2.576413132701720,
            -2.675593912859120,
            -2.8153300703866400,
            -2.897475156648710,
            -3.0042154567120800,
            -3.087886507216510,
        ]
    ),
    "E_power": np.array(
        [
            0.0006454923080882520,
            0.004277023425898170,
            0.004719543572446050,
            0.00909714298712173,
            0.01344187403986090,
            0.01896059941058610,
            0.02101937603827090,
            0.022335737104727900,
            0.026029715703921800,
            0.02917012237740640,
            0.0338523245243911,
            0.03820738749745440,
            0.04733370907740660,
            0.05269608016472140,
            0.053790874615060400,
            0.05697639350179900,
            0.06435487925718170,
            0.07909179437004270,
            0.08958910066880800,
            0.1051813598911370,
            0.11617812912435900,
        ]
    ),
}


class TestHydration(unittest.TestCase):
    def setUp(self):
        self.data = {
            "E_array": TESTSET["E"],
            "E_powers": TESTSET["E_power"],
            "T1_array": TESTSET["T1"],
            "T1_powers": TESTSET["T1_power"],
            "field": 348.5,
            "spin_C": 125,
            "T10": 2.0,
            "T100": 2.5,
            "smax_model": "tethered",
            "interpolate_method": "second_order",
        }

        self.constants = {
            "ksigma_bulk": 95.4,
            "krho_bulk": 353.4,
            "klow_bulk": 366,
            "tcorr_bulk": 54,
            "D_H2O": 2.3e-9,
            "D_SL": 4.1e-10,
            "delta_T1_water": False,
            "T1_water": False,
            "macro_C": False,
        }

        self.ws = dnp.create_workspace("hydration_inputs", self.data)
        self.ws.add("hydration_constants", self.constants)

    def test_hydration_return_dict(self):

        result = dnp.dnpHydration.hydration(self.ws)

        self.assertEqual(len(self.data["E_powers"]), 21)
        self.assertEqual(len(self.data["E_array"]), 21)
        self.assertEqual(len(self.data["T1_powers"]), 5)
        self.assertEqual(len(self.data["T1_array"]), 5)
        self.assertEqual(len(result["interpolated_T1"]), 21)

        self.assertAlmostEqual(min(result["interpolated_T1"]), 2.051727288206873, places = 6)
        self.assertAlmostEqual(max(result["interpolated_T1"]), 2.5383409868931706)
        self.assertEqual(len(result["ksigma_array"]), 21)
        self.assertAlmostEqual(min(result["ksigma_array"]), 2.5002212753387965, places = 6)
        self.assertAlmostEqual(max(result["ksigma_array"]), 19.573744921860015, places = 6)
        self.assertEqual(len(result["ksigma_fit"]), 21)
        self.assertAlmostEqual(min(result["ksigma_fit"]), 1.9391281665269349, places = 6)
        self.assertAlmostEqual(max(result["ksigma_fit"]), 19.17592802884648, places = 6)
        self.assertEqual(len(result["uncorrected_Ep"]), 21)
        self.assertAlmostEqual(min(result["uncorrected_Ep"]), -2.90847384506599, places = 6)
        self.assertAlmostEqual(max(result["uncorrected_Ep"]), 0.7434984782673344, places = 6)

        self.assertAlmostEqual(result["ksigma"], 20.17803815171555, places = 6)
        self.assertAlmostEqual(result["klow"], 1286.2512443126634, places = 6)
        self.assertAlmostEqual(result["tcorr"], 485.8952027395348, places = 6)
        self.assertAlmostEqual(result["Dlocal"], 3.01176054373284e-10, places = 6)

        self.assertTrue(
            [x in self.ws["hydration_results"].keys() for x in result.keys()]
        )

    def test_hydration_inplace(self):

        self.ws["hydration_inputs"]["spin_C"] = 100
        self.ws["hydration_inputs"]["smax_model"] = "free"
        self.ws["hydration_inputs"]["interpolate_method"] = "linear"
        self.ws["hydration_constants"]["D_SL"] = 4.1e-9

        dnp.dnpHydration.hydration(self.ws)

        self.assertEqual(len(self.ws["hydration_results"]["interpolated_T1"]), 21)
        self.assertAlmostEqual(
            min(self.ws["hydration_results"]["interpolated_T1"]), 2.0978389591800344, places = 6
        )
        self.assertAlmostEqual(
            max(self.ws["hydration_results"]["interpolated_T1"]), 2.5855460713039324, places = 6
        )
        self.assertEqual(len(self.ws["hydration_results"]["ksigma_array"]), 21)
        self.assertAlmostEqual(
            min(self.ws["hydration_results"]["ksigma_array"]), 3.0565812706454305, places = 6
        )
        self.assertAlmostEqual(
            max(self.ws["hydration_results"]["ksigma_array"]), 24.020476541485717, places = 6
        )
        self.assertEqual(len(self.ws["hydration_results"]["ksigma_fit"]), 21)
        self.assertAlmostEqual(
            min(self.ws["hydration_results"]["ksigma_fit"]), 2.3390612006821225, places = 6
        )
        self.assertAlmostEqual(
            max(self.ws["hydration_results"]["ksigma_fit"]), 24.280175919600392, places = 6
        )
        self.assertEqual(len(self.ws["hydration_results"]["uncorrected_Ep"]), 21)
        self.assertAlmostEqual(
            min(self.ws["hydration_results"]["uncorrected_Ep"]), -2.9084738857665413, places = 6
        )
        self.assertAlmostEqual(
            max(self.ws["hydration_results"]["uncorrected_Ep"]), 0.7434984946695101, places = 6
        )

        self.assertAlmostEqual(
            self.ws["hydration_results"]["ksigma"], 73.98621213840886, places = 6
        )
        self.assertAlmostEqual(self.ws["hydration_results"]["klow"], 1494.0321716770457, places = 6)
        self.assertAlmostEqual(self.ws["hydration_results"]["tcorr"], 230.5905102559904, places = 6)

        self.assertAlmostEqual(
            self.ws["hydration_results"]["Dlocal"], 1.4987607235715452e-09, places = 6
        )
