import unittest
import numpy as np
from ..hydration.scratch import *
import os

simple_set_path = os.path.join('data', 'test_hydration', 'simple_set_1')


class TestHydration(unittest.TestCase):
    def setUp(self):
        self.E_powers = np.genfromtxt(os.path.join(simple_set_path, 'E_powers.txt'))
        self.Ep = np.genfromtxt(os.path.join(simple_set_path, 'Ep.txt'))
        self.T1_powers = np.genfromtxt(os.path.join(simple_set_path, 'T1_powers.txt'))
        self.T1p = np.genfromtxt(os.path.join(simple_set_path, 'T1p.txt'))
        self.conc = 100
        self.k_sig_bulk = 95.4
        self.T100 = 2.5

    def test_getT1p(self):
        # TODO: implement more assertions
        t1p_fun = getT1p(self.T1p, self.T1_powers)
        T10 = t1p_fun(0)
        self.assertAlmostEqual(T10, 1.33)

        mis_matched_T1p = np.array([range(len(self.T1_powers)-1)])
        self.assertRaises(AssertionError, getT1p, mis_matched_T1p, self.T1_powers)

    def test_getT1CorrectedEnhancement(self):
        # TODO: implement this
        t1p_fun = getT1p(self.T1p, self.T1_powers)
        e_corr = getT1CorrectedEnhancement(self.Ep, self.E_powers, t1p_fun)
        pass

    def test_getKsigma(self):
        t1p_fun = getT1p(self.T1p, self.T1_powers)
        e_corr = getT1CorrectedEnhancement(self.Ep, self.E_powers, t1p_fun)
        ksigma = getKsigma(e_corr, self.E_powers)

        T10 = t1p_fun(0)
        k_rho = ((1 / T10) - (1 / self.T100)) / self.conc
        ksi = ksigma / k_rho

        self.assertAlmostEqual(ksi, 0.0326)

        JH = wip_getJH(ksi)
        self.assertAlmostEqual(JH, 0.94)



if __name__ == '__main__':
    unittest.main()
