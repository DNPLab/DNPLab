import sys
sys.path.append('.../odnplab')

import unittest
import numpy as np
import odnpLab
from odnpLab.hydration import HydrationCalculator, HydrationParameter


class TestHydration(unittest.TestCase):
    def setUp(self):
        self.T1p = np.array([2.0850606305680700
        ,2.1588239679329700
        ,2.246270535391770
        ,2.3131083281827900
        ,2.361693123732750])
        self.T1_powers = np.array([0.0005813884968087980
        ,0.023929259151532000
        ,0.05364586645098500
        ,0.08510236988564890
        ,0.11448695555368400])
        self.Ep = np.array([0.5244962257098100
        ,-1.0460060878468600
        ,-1.1923729834996700
        ,-2.0283543307318300
        ,-2.5623585026035300
        ,-3.1083332432507900
        ,-3.2258438385993300
        ,-3.347961336229490
        ,-3.67138811609969
        ,-3.7817316135687500
        ,-4.032510453619460
        ,-4.2417758835183300
        ,-4.552046251816710
        ,-4.733158273060310
        ,-4.777088808443740
        ,-4.8489491678359600
        ,-5.024100136084710
        ,-5.324997995063750
        ,-5.51186544738343
        ,-5.807764442039210
        ,-5.9835905676571500])
        self.E_powers = np.array([0.0006370105281515720
        ,0.004214559147178920
        ,0.004652040210109040
        ,0.008966203931561900
        ,0.013256981804212100
        ,0.018721615363004000
        ,0.02068371710482890
        ,0.022057273192711100
        ,0.025735671896052500
        ,0.028814716528459300
        ,0.03343506130920710
        ,0.03775305788018570
        ,0.04682777569441930
        ,0.052001376894013300
        ,0.05321444208836340
        ,0.05628331436209730
        ,0.06358027888841540
        ,0.07825239654597910
        ,0.08867115988141430
        ,0.1042006843923840
        ,0.11448695555368400])

        hp = HydrationParameter()
        hp.field = 348.5
        hp.slC = 125e-6
        hp.T100 = 2.0
        hp.T10 = 1.5
        self.hp = hp
        self.hc = HydrationCalculator(T1=self.T1p, T1_power=self.T1_powers,
                                 E=self.Ep, E_power=self.E_powers,
                                 hp=hp)

    def test_interpT1_linear_almost_1p99(self):
        self.hc.hp.t1InterpMethod = 'linear'
        self.hc.run()
        self.assertAlmostEqual(self.hc.results.T1interp[0], 2.097, places=3)

    def test_interpT1_2ord_almost_1p99(self):
        self.hc.hp.t1InterpMethod = '2ord'
        self.hc.run()
        self.assertAlmostEqual(self.hc.results.T1interp[0], 2.085, places=3)

    def _run_2ord(self):
        self.hc.hp.smaxMod = 'tethered'
        self.hc.hp.t1InterpMethod = '2ord'
        self.hc.run()
        
    def _run_linear(self):
        self.hc.hp.smaxMod = 'tethered'
        self.hc.hp.t1InterpMethod = 'linear'
        self.hc.run()
    
    def test_2ord_krho(self):
        self._run_2ord()
        self.assertAlmostEqual(self.hc.results.k_rho, 1333.33, places=2)
    
    def test_2ord_ksigma(self):
        self._run_2ord()
        self.assertAlmostEqual(self.hc.results.k_sigma, 38.04, places=2)
    
    def test_2ord_ksi(self):
        self._run_2ord()
        self.assertAlmostEqual(self.hc.results.ksi, 0.0285, places=4)

    def test_2ord_klow(self):
        self._run_2ord()
        self.assertAlmostEqual(self.hc.results.k_low, 2133.46, places=2)

    def test_2ord_tcorr(self):
        self._run_2ord()
        self.assertAlmostEqual(self.hc.results.tcorr, 448.97, places=2)

    def test_2ord_Dlocal(self):
        self._run_2ord()
        self.assertAlmostEqual(self.hc.results.dLocal, 3.26e-10, places=12)
        
    def test_linear_krho(self):
        self._run_linear()
        self.assertAlmostEqual(self.hc.results.k_rho, 1333.33, places=2)
    
    def test_linear_ksigma(self):
        self._run_linear()
        self.assertAlmostEqual(self.hc.results.k_sigma, 38.20, places=2)
    
    def test_linear_ksi(self):
        self._run_linear()
        self.assertAlmostEqual(self.hc.results.ksi, 0.0286, places=4)

    def test_linear_klow(self):
        self._run_linear()
        self.assertAlmostEqual(self.hc.results.k_low, 2133.09, places=2)

    def test_linear_tcorr(self):
        self._run_linear()
        self.assertAlmostEqual(self.hc.results.tcorr, 447.77, places=2)

    def test_linear_Dlocal(self):
        self._run_linear()
        self.assertAlmostEqual(self.hc.results.dLocal, 3.27e-10, places=12)

if __name__ == '__main__':
    unittest.main()
