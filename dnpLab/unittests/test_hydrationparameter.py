import sys
sys.path.append('.../odnplab')
import unittest
from dnpLab.dnpHydration import HydrationParameter


class TestHydrationParameter(unittest.TestCase):

    def test_default_values(self):
    
        hp = HydrationParameter()
        
        self.assertEqual(hp.field, 348.5)
        self.assertEqual(hp.spin_C, 100)
        self.assertEqual(hp.smax_model, 'tethered')
        self.assertEqual(hp.ksigma_bulk, 95.4)
        self.assertEqual(hp.T10, 1.5)
        self.assertEqual(hp.T100, 2.5)
        self.assertEqual(hp.tcorr_bulk, 54)
        self.assertEqual(hp.D_H2O, 2.3e-9)
        self.assertEqual(hp.D_SL, 4.1e-10)
        self.assertEqual(hp.klow_bulk, 366)
        self.assertEqual(hp.t1_interp_method, 'second_order')
        
    def test_struct_like(self):
    
        hp = HydrationParameter()
        
        hp.field = 300
        hp.spin_C = 200
        hp.ksigma_bulk = 95
        hp.T10 = 1
        hp.T100 = 2
        hp.tcorr_bulk = 50
        hp.D_H2O = 3e-10
        hp.D_SL = 5e-10
        hp.klow_bulk = 350
        hp.smax_model = 'free'
        hp.t1_interp_method = 'second_order'
        
        self.assertEqual(hp['field'], hp.field)
        self.assertEqual(hp['spin_C'], hp.spin_C)
        self.assertEqual(hp['ksigma_bulk'], hp.ksigma_bulk)
        self.assertEqual(hp['T10'], hp.T10)
        self.assertEqual(hp['T100'], hp.T100)
        self.assertEqual(hp['tcorr_bulk'], hp.tcorr_bulk)
        self.assertEqual(hp['D_H2O'], hp.D_H2O)
        self.assertEqual(hp['D_SL'], hp.D_SL)
        self.assertEqual(hp['klow_bulk'], hp.klow_bulk)
        self.assertEqual(hp['smax_model'], hp.smax_model)
        self.assertEqual(hp['t1_interp_method'], hp.t1_interp_method)
        
        with self.assertRaises(ValueError):
            hp.smax_model = 'notinlist'
        with self.assertRaises(ValueError):
            hp.t1_interp_method = 'notinlist'
        
    def test_dict_like(self):
    
        hp = HydrationParameter()
        
        hp['field'] = 300
        hp['spin_C'] = 200
        hp['ksigma_bulk'] = 95
        hp['T10'] = 1
        hp['T100'] = 2
        hp['tcorr_bulk'] = 50
        hp['D_H2O'] = 3e-10
        hp['D_SL'] = 5e-10
        hp['klow_bulk'] = 350
        hp['smax_model'] = 'free'
        hp['t1_interp_method'] = 'second_order'
        
        self.assertEqual(hp['field'], hp.field)
        self.assertEqual(hp['spin_C'], hp.spin_C)
        self.assertEqual(hp['ksigma_bulk'], hp.ksigma_bulk)
        self.assertEqual(hp['T10'], hp.T10)
        self.assertEqual(hp['T100'], hp.T100)
        self.assertEqual(hp['tcorr_bulk'], hp.tcorr_bulk)
        self.assertEqual(hp['D_H2O'], hp.D_H2O)
        self.assertEqual(hp['D_SL'], hp.D_SL)
        self.assertEqual(hp['klow_bulk'], hp.klow_bulk)
        self.assertEqual(hp['smax_model'], hp.smax_model)
        self.assertEqual(hp['t1_interp_method'], hp.t1_interp_method)

        with self.assertRaises(ValueError):
            hp['smax_model'] = 'notinlist'
        with self.assertRaises(ValueError):
            hp['t1_interp_method'] = 'notinlist'


if __name__ == '__main__':
    unittest.main()
