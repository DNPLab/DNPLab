
import unittest
from dnplab import dnpdata
from numpy.testing import assert_array_equal
import dnplab.dnpImport.vnmrj as vnmrj
import dnplab.dnpNMR as nmr
import dnplab as dnp
import numpy as np


class dnpNMR_tester(unittest.TestCase):
    def setUp(self):
        path = './data/vnmrj/10mM_tempol_in_water_array.fid'
        self.data = vnmrj.import_vnmrj(path)
        self.ws = dnp.create_workspace()
        self.ws['raw'] = self.data
        self.ws.copy('raw','proc')

    def test_basic_nmr_processing(self):

        self.ws = nmr.remove_offset(self.ws, {})
        self.ws = nmr.window(self.ws, {})
        self.ws = nmr.fourier_transform(self.ws, {})
        self.ws = nmr.autophase(self.ws, {})

    def test_integrate(self):
        values = self.ws['proc'].values
        # Doing so must not change the workspace at all
        data = self.ws['proc']
        nmr.integrate(data, {})
        np.testing.assert_array_equal(self.ws['proc'].values, values)


class dnpNMR_tester_sim(unittest.TestCase):
    def setUp(self):
        p1 = np.array([0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0])
        p2 = np.array([0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0, 0, 0])
        p3 = np.array([0, 0, 0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0])
        self.data = dnpdata(np.array([p1,p2,p3]).T, [np.arange(0, len(p1)), np.arange(1,3)], ['x', 't2'])
        self.ws = dnp.create_workspace()
        self.ws['raw'] = self.data
        self.ws.copy('raw','proc')

    def test_align(self):
        nmr.align(self.ws, {'dim': 'x'})
        assert_array_equal(self.ws['proc'].values, np.array([
            [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0],
            [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0],
            [0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 0, 0, 0]
        ]).T)
