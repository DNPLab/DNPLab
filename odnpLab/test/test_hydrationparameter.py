import sys
sys.path.append('.../odnplab')

import unittest
import odnpLab
from odnpLab.hydration import HydrationParameter


class TestHydrationParameter(unittest.TestCase):

    def test_default_values(self):
        hp = HydrationParameter()
        self.assertEqual(hp.field, 348.5)
        self.assertEqual(hp.slC, 100)
        self.assertEqual(hp.smaxMod, 'tethered')
        self.assertEqual(hp.ksig_bulk, 95.4)
        self.assertEqual(hp.T10, 1.5)
        self.assertEqual(hp.T100, 2.5)
        self.assertEqual(hp.tcorr_bulk, 54)
        self.assertEqual(hp.dH2O, 2.3e-9)
        self.assertEqual(hp.dSL, 4.1e-10)
        self.assertEqual(hp.k_low_bulk, 366)
        self.assertEqual(hp.t1InterpMethod, '2ord')
        
    def test_struct_like(self):
        hp = HydrationParameter()
        hp.field = 400
        hp.smaxMod = 'free'
        self.assertEqual(hp['field'], hp.field)
        self.assertEqual(hp['smaxMod'], hp.smaxMod)
        
        with self.assertRaises(ValueError):
            hp.smaxMod = 'incorrect'
        
    def test_dict_like(self):
        hp = HydrationParameter()
        hp['field'] = 300
        hp['smaxMod'] = 'tethered'
        self.assertEqual(hp.field, hp['field'])
        self.assertEqual(hp.smaxMod, hp['smaxMod'])

        with self.assertRaises(ValueError):
            hp['smaxMod'] = 'incorrect'
            
if __name__ == '__main__':
    unittest.main()
