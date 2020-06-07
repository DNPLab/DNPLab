import unittest
import odnpLab.odnpImport.bruker as bruker
import odnpLab.odnpImport.kea as kea


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
        self.assertEqual(data.axesLabels[0], 't')
        self.assertEqual(data.data.size, 8147)
        self.assertAlmostEqual(data.data.min(), -5-4.168734491315137j)
        self.assertAlmostEqual(data.params['nmrFreq'], 14831413.270000001)

    def test_importBruker_exp5_is_2d_phcyc(self):
        data = bruker.importBruker(self.testdata, expNum=5)
        self.assertEqual(data.data.shape[0], 11912)
        self.assertEqual(data.axesLabels, ['t'])
        self.assertAlmostEqual(data.params['nmrFreq'], 14831413.270000001)
        self.assertAlmostEqual(data.data[365], -0.182861328125-0.71875j)

    def test_importBruker_exp28_is_2d(self):
        data = bruker.importBruker(self.testdata, expNum=28)
        self.assertEqual(data.data.shape, (7922, 8))
        self.assertEqual(data.axesLabels, ['t', 't1'])
        self.assertAlmostEqual(data.params['nmrFreq'], 14831413.270000001)
        self.assertAlmostEqual(data.data[365, 6], -0.110595703125+0.47705078125j)


class KEAImportTester(unittest.TestCase):
    def setup(self):
        self.test_data = './data/kea/toluene_10mM_Tempone'
    def test_importKEA(self):
        data = kea.importKea(self.test_data)

if __name__ == '__main__':
    unittest.main()
    pass
