import unittest
import dnpLab.dnpImport.bruker as bruker
import dnpLab.dnpImport.kea as kea
import dnpLab.dnpImport.varian as varian


class ImportBrukerTester(unittest.TestCase):
    def setUp(self):
        self.testdata = './data/test_set/'

    def test_dirDataType(self):
        self.assertEqual('fid', bruker.dirDataType(self.testdata, expNum=1))
        self.assertEqual('fid', bruker.dirDataType(self.testdata, expNum=2))
        self.assertEqual('fid', bruker.dirDataType(self.testdata, expNum=3))
        self.assertEqual('',    bruker.dirDataType(self.testdata, expNum=4))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=5))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=6))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=7))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=8))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=9))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=10))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=11))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=12))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=13))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=14))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=15))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=16))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=17))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=18))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=19))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=20))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=21))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=22))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=23))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=24))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=25))
        self.assertEqual('serPhaseCycle', bruker.dirDataType(self.testdata, expNum=26))
        self.assertEqual('',    bruker.dirDataType(self.testdata, expNum=27))
        self.assertEqual('ser', bruker.dirDataType(self.testdata, expNum=28))
        self.assertEqual('ser', bruker.dirDataType(self.testdata, expNum=29))
        self.assertEqual('ser', bruker.dirDataType(self.testdata, expNum=30))
        self.assertEqual('ser', bruker.dirDataType(self.testdata, expNum=31))
        self.assertEqual('ser', bruker.dirDataType(self.testdata, expNum=32))
        self.assertEqual('',    bruker.dirDataType(self.testdata, expNum=33))

    def test_importBruker_exp1_is_fid(self):
        data = bruker.importBruker(self.testdata, expNum=1)
        self.assertEqual(data.dims[0], 't')
        self.assertEqual(data.values.size, 8147)
        self.assertAlmostEqual(data.values.min(), -5-4.168734491315137j)
        self.assertAlmostEqual(data.attrs['nmrFreq'], 14831413.270000001)

    def test_importBruker_exp5_is_2d_phcyc(self):
        data = bruker.importBruker(self.testdata, expNum=5)
        self.assertEqual(data.values.shape[0], 11912)
        self.assertEqual(data.dims, ['t'])
        self.assertAlmostEqual(data.attrs['nmrFreq'], 14831413.270000001)
        self.assertAlmostEqual(data.values[365], -0.182861328125-0.71875j)

    def test_importBruker_exp28_is_2d(self):
        data = bruker.importBruker(self.testdata, expNum=28)
        self.assertEqual(data.values.shape, (7922, 8))
        self.assertEqual(data.dims, ['t', 't1'])
        self.assertAlmostEqual(data.attrs['nmrFreq'], 14831413.270000001)
        self.assertAlmostEqual(data.values[365, 6], -0.110595703125+0.47705078125j)


class KEAImportTester(unittest.TestCase):
    def setUp(self):
        self.test_data = './data/kea/toluene_10mM_Tempone'

    def test_importKEA_exp_is_1d(self):
        datas = [kea.importKea(self.test_data + '/%i/data.csv'%expNum) \
                 for expNum in[1, 21, 42]]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (16384, ))
            self.assertEqual(data.dims, ['t'])
            self.assertAlmostEqual(data.attrs['nmrFreq'], 14244500.0)
        self.assertAlmostEqual(datas[0].values[365], -0.217937+0.24907j)
        self.assertAlmostEqual(datas[1].values[365], 0.0400292-0.0756107j)
        self.assertAlmostEqual(datas[2].values[365], 1.09858-2.57966j)


class VNMRJImportTester(unittest.TestCase):
    def setUp(self):
        self.test_data2Ds = ['./data/vnmrj/'+s for s in [
            '10mM_tempol_in_water_array.fid'
        ]]
        self.test_data1Ds = ['./data/vnmrj/' + s for s in [
            '10mM_tempol_in_water_mw_40dBm.fid',
            '10mM_tempol_in_water_mw_off.fid'
        ]]

    def test_importVNMRJ1D(self):
        datas = [varian.importVarian(path=path) for path in self.test_data1Ds]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (131072, ))
            self.assertEqual(data.dims, ['t', ])
            self.assertAlmostEqual(data.attrs['nmrFreq'], 14244283.4231)
        self.assertAlmostEqual(datas[0].values[365], (-20378767-2734659j))
        self.assertAlmostEqual(datas[1].values[365], (-950662+138458j))

    def test_importVNMRJ2D(self):
        datas = [varian.importVarian(path=path) for path in self.test_data2Ds]
        for i, data in enumerate(datas):
            self.assertEqual(data.values.shape, (131072, 5))
            self.assertEqual(data.dims, ['t', 'x'])
            self.assertAlmostEqual(data.attrs['nmrFreq'], 14244283.4231)
        self.assertAlmostEqual(datas[0].values[365, 3], (-1263136+1063328.5j))


if __name__ == '__main__':
    unittest.main()
    pass
