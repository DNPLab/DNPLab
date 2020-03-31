import unittest
import numpy as np
from odnpLab.hydration import HydrationCalculator, HydrationParameter


class TestHydration(unittest.TestCase):
    def setUp(self):
        self.T1p = np.array([1.995006087,2.059663367,2.149840696,2.245130396,2.302170516])
        self.T1_powers = np.array([0.00062242,0.025845202,0.058096435,0.092231066,0.123976276])
        self.Ep = np.array([0.296341264,-1.648327271,-1.754501023,-2.69317187,-3.3022277,-3.810733497,-4.010433129,-4.108039458,-4.353963952,-4.541108629,-4.770356192,-4.938621066,-5.25776408,-5.372534317,-5.445489014,-5.529224875,-5.695695725,-5.989497644,-6.143749818,-6.353057052,-6.541485408])
        self.E_powers = np.array([0.000643822,0.004483313,0.004731359,0.0091599,0.013577705,0.019220847,0.02205228,0.022732423,0.02664021,0.029912122,0.034827366,0.039436018,0.04900578,0.055400375,0.055803875,0.059311184,0.067142499,0.082665971,0.093849966,0.110199413,0.123976276])

        hp = HydrationParameter()
        hp.field = 348.5
        hp.slC = 200e-6
        hp.T100 = 2.50
        hp.T10 = 1.95
        hp.ksig_bulk = 95.4
        self.hp = hp
        self.hc = HydrationCalculator(T1=self.T1p, T1_power=self.T1_powers,
                                 E=self.Ep, E_power=self.E_powers,
                                 hp=hp)

    def test_interpT1_linear_almost_1p99(self):
        self.hc.hp.t1InterpMethod = 'linear'
        self.hc.run()
        self.assertAlmostEqual(self.hc.results.T1fit[0], 1.995006087, places=2)

    def test_interpT1_2ord_almost_1p99(self):
        self.hc.hp.t1InterpMethod = '2ord'
        self.hc.run()
        self.assertAlmostEqual(self.hc.results.T1fit[0], 1.99, places=2)

    def test_interpT1_2ord_smaller_than_linear(self):
        """Second order interpolation should gives smaller T1[0] than linear?"""
        self.hc.hp.t1InterpMethod = '2ord'
        self.hc.run()
        second = self.hc.results.T1fit[0]
        self.hc.hp.t1InterpMethod = 'linear'
        self.hc.run()
        linear = self.hc.results.T1fit[0]
        self.assertLess(second, linear)

    def _run_tethered_2ord(self):
        self.hc.hp.smaxMod = 'tethered'
        self.hc.hp.t1InterpMethod = '2ord'
        self.hc.hp.disableJRot()
        self.hc.run()

    def test_no_trot_tethered_ksigma_is_25p99(self):
        self._run_tethered_2ord()
        self.assertAlmostEqual(self.hc.results.k_sigma, 25.99, places=2)

    def test_no_trot_tethere_ksigma_bulk_ratio_is_3p67(self):
        self._run_tethered_2ord()
        self.assertAlmostEqual(self.hc.results.ksigma_kbulk_invratio, 3.67, places=2)

    def test_no_trot_tethere_klow_is_879p5(self):
        self._run_tethered_2ord()
        self.assertAlmostEqual(self.hc.results.k_low, 879.5, places=1)

    def test_no_trot_tethere_tcorr_is_326p0(self):
        self._run_tethered_2ord()
        self.assertAlmostEqual(self.hc.results.tcorr, 326.0, places=1)

    def test_no_trot_tethere_Dlocal_is_4p49Em10(self):
        self._run_tethered_2ord()
        self.assertAlmostEqual(self.hc.results.dLocal, 4.49e-10, places=12)

    def test_kpho_is_564p1(self):
        self._run_tethered_2ord()
        self.assertAlmostEqual(self.hc.results.k_rho, 564.1, places=1)

    def _run_expert_mode(self):
        """Expert Mode"""
        self.hc.hp.smaxMod = 'tethered'
        self.hc.hp.t1InterpMethod = '2ord'
        self.hc.hp.enableJRot(tauRot=0.020, percentBound=50)
        self.hc.run()

    def test_expert_ksigma_is_25p53(self):
        self._run_expert_mode()
        self.assertAlmostEqual(self.hc.results.k_sigma, 25.53)


if __name__ == '__main__':
    unittest.main()
