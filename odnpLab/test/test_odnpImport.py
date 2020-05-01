import unittest
import odnpLab.odnpImport.bruker as bruker


class ODNPImportBrukerTester(unittest.TestCase):
    def setUp(self):
        self.testdata = './data/20190821_TW_4OH-TEMPO_500uM/'

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

    def test_bruker_import(self):
        data = bruker.importBruker(self.testdata, expNum=1)
        print(data)
