import unittest
import numpy as np
from odnpLab.hydration import HydrationCalculator, HydrationParameter
import os

simple_set_path = os.path.join('data', 'test_hydration', 'simple_set_1')


class TestHydration(unittest.TestCase):
    def setUp(self):
        self.T1p = np.array([1.995006087,2.059663367,2.149840696,2.245130396,2.302170516])
        self.T1_powers = np.array([0.00062242,0.025845202,0.058096435,0.092231066,0.123976276])
        self.Ep = np.array([1.992144472,2.00329461,2.004013328,2.016809974,2.029505447,2.045611122,2.053641678,2.055565442,2.066577408,2.07574152,2.089406727,2.102102039,2.128069065,2.145097958,2.146163264,2.155375272,2.17562144,2.214312695,2.240869791,2.277460538,2.306022])
        self.E_powers = np.array([0.000643822,0.004483313,0.004731359,0.0091599,0.013577705,0.019220847,0.02205228,0.022732423,0.02664021,0.029912122,0.034827366,0.039436018,0.04900578,0.055400375,0.055803875,0.059311184,0.067142499,0.082665971,0.093849966,0.110199413,0.123976276])

        # self.E_powers = np.genfromtxt(os.path.join(simple_set_path, 'E_powers.txt'))
        # self.Ep = np.genfromtxt(os.path.join(simple_set_path, 'Ep.txt'))
        # self.T1_powers = np.genfromtxt(os.path.join(simple_set_path, 'T1_powers.txt'))
        # self.T1p = np.genfromtxt(os.path.join(simple_set_path, 'T1p.txt'))

        hp = HydrationParameter()
        hp.slC = 200
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

    def test_no_trot_tethered_ksigma_is_25p53(self):
        self.hc.hp.smaxMod = 'tethered'
        self.hc.hp.t1InterpMethod = 'linear'
        self.hc.run()
        self.assertAlmostEqual(self.hc.results.ksigma, 25.53)

    # def test_T10_is_1p33(self):
    #     # TODO: implement more assertions
    #     T10 = self.hc.T1fit[0]
    #     self.assertAlmostEqual(T10, 1.33)
    #
    # def test_ksi_is_0p0326(self):
    #     self.assertAlmostEqual(self.hc.results.ksi, 0.0326)


if __name__ == '__main__':
    unittest.main()
