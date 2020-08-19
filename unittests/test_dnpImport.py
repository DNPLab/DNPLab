import unittest
import dnpLab.dnpImport.topspin as topspin
import dnpLab.dnpImport.prospa as prospa
import dnpLab.dnpImport.vnmrj as vnmrj
import os
from numpy.testing import assert_array_equal

class import_topspin_tester(unittest.TestCase):
    def setUp(self):
        self.testdata = './data/topspin/'

    def test_dir_data_type(self):
        self.assertEqual('fid', topspin.dir_data_type(self.testdata, expNum=1))
        self.assertEqual('fid', topspin.dir_data_type(self.testdata, expNum=2))
        self.assertEqual('fid', topspin.dir_data_type(self.testdata, expNum=3))
        self.assertEqual('',    topspin.dir_data_type(self.testdata, expNum=4))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=5))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=6))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=7))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=8))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=9))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=10))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=11))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=12))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=13))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=14))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=15))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=16))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=17))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=18))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=19))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=20))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=21))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=22))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=23))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=24))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=25))
        self.assertEqual('serPhaseCycle', topspin.dir_data_type(self.testdata, expNum=26))
        self.assertEqual('',    topspin.dir_data_type(self.testdata, expNum=27))
        self.assertEqual('ser', topspin.dir_data_type(self.testdata, expNum=28))
        self.assertEqual('ser', topspin.dir_data_type(self.testdata, expNum=29))
        self.assertEqual('ser', topspin.dir_data_type(self.testdata, expNum=30))
        self.assertEqual('ser', topspin.dir_data_type(self.testdata, expNum=31))
        self.assertEqual('ser', topspin.dir_data_type(self.testdata, expNum=32))
        self.assertEqual('',    topspin.dir_data_type(self.testdata, expNum=33))

    def test_import_topspin_exp1_is_fid(self):
        data = topspin.import_topspin(self.testdata, expNum=1)
        self.assertEqual(data.dims[0], 't2')
        self.assertEqual(data.values.size, 8147)
        self.assertAlmostEqual(data.values.min(), -5-4.168734491315137j)
        self.assertAlmostEqual(data.attrs['nmr_frequency'], 14831413.270000001)

    def test_import_topspin_exp5_is_2d_phcyc(self):
        data = topspin.import_topspin(self.testdata, expNum=5)
        self.assertEqual(data.values.shape[0], 11912)
        self.assertEqual(data.dims, ['t2'])
        self.assertAlmostEqual(data.attrs['nmr_frequency'], 14831413.270000001)
        self.assertAlmostEqual(data.values[365], -0.182861328125-0.71875j)

    def test_import_topspin_exp28_is_2d(self):
        data = topspin.import_topspin(self.testdata, expNum=28)
        self.assertEqual(data.values.shape, (7922, 8))
        self.assertEqual(data.dims, ['t2', 't1'])
        self.assertAlmostEqual(data.attrs['nmr_frequency'], 14831413.270000001)
        self.assertAlmostEqual(data.values[365, 6], -0.110595703125+0.47705078125j)
    
    def test_import_topspin_jcamp_dx(self):
        attrs = topspin.topspin_jcamp_dx(os.path.join(self.testdata,'1','acqus'))
        self.assertEqual(attrs['DIGTYP'], 9)
        self.assertAlmostEqual(attrs['O1'], 1413.27)
        assert_array_equal(attrs['XGAIN'], [0, 0, 0, 0])
        assert_array_equal(attrs['TPOAL'], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])


class prospa_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data = './data/prospa/toluene_10mM_Tempone'

    def test_import_prospa_exp_is_1d(self):
        datas = [prospa.import_prospa(self.test_data + '/%i/data.csv'%expNum) \
                 for expNum in[1, 21, 42]]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (16384, ))
            self.assertEqual(data.dims, ['t2'])
            self.assertAlmostEqual(data.attrs['nmr_frequency'], 14244500.0)
        self.assertAlmostEqual(datas[0].values[365], -0.217937+0.24907j)
        self.assertAlmostEqual(datas[1].values[365], 0.0400292-0.0756107j)
        self.assertAlmostEqual(datas[2].values[365], 1.09858-2.57966j)


class vnmrj_import_tester(unittest.TestCase):
    def setUp(self):
        self.test_data2Ds = ['./data/vnmrj/'+s for s in [
            '10mM_tempol_in_water_array.fid'
        ]]
        self.test_data1Ds = ['./data/vnmrj/' + s for s in [
            '10mM_tempol_in_water_mw_40dBm.fid',
            '10mM_tempol_in_water_mw_off.fid'
        ]]

    def test_import_vnmrj_1d(self):
        datas = [vnmrj.import_vnmrj(path=path) for path in self.test_data1Ds]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (131072, ))
            self.assertEqual(data.dims, ['t2', ])
            self.assertAlmostEqual(data.attrs['nmr_frequency'], 14244283.4231)
        self.assertAlmostEqual(datas[0].values[365], (-20378767-2734659j))
        self.assertAlmostEqual(datas[1].values[365], (-950662+138458j))

    def test_import_vnmrj_2d(self):
        datas = [vnmrj.import_vnmrj(path=path) for path in self.test_data2Ds]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (131072, 5))
            self.assertEqual(data.dims, ['t2', 'array'])
            self.assertAlmostEqual(data.attrs['nmr_frequency'], 14244283.4231)
        self.assertAlmostEqual(datas[0].values[365, 3], (-1263136+1063328.5j))


if __name__ == '__main__':
    unittest.main()
    pass
