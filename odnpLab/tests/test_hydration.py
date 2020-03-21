import unittest
import numpy as np
from odnpLab.hydration import HydrationCalculator, HydrationParameter
import os

simple_set_path = os.path.join('data', 'test_hydration', 'simple_set_1')


class TestHydration(unittest.TestCase):
    def setUp(self):
        self.E_powers = np.genfromtxt(os.path.join(simple_set_path, 'E_powers.txt'))
        self.Ep = np.genfromtxt(os.path.join(simple_set_path, 'Ep.txt'))
        self.T1_powers = np.genfromtxt(os.path.join(simple_set_path, 'T1_powers.txt'))
        self.T1p = np.genfromtxt(os.path.join(simple_set_path, 'T1p.txt'))

        hp = HydrationParameter()
        hp.slC = 100
        hp.ksig_bulk = 95.4
        hp.T100 = 2.5
        hp.smaxMod = 'tethered'
        hp.fitopt = 'linear'
        hc = HydrationCalculator(T1=self.T1p, T1_power=self.T1_powers,
                                 E=self.Ep, E_power=self.E_powers,
                                 hp=hp)

        self.hc, self.hp = hc, hp

    def test_T10_is_1p33(self):
        # TODO: implement more assertions
        T10 = self.hc.T1fit[0]
        self.assertAlmostEqual(T10, 1.33)

    def test_ksi_is_0p0326(self):
        self.assertAlmostEqual(self.hc.results.ksi, 0.0326)


if __name__ == '__main__':
    unittest.main()
