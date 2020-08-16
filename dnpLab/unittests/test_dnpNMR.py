
import unittest
from numpy.testing import assert_array_equal
from dnpLab.core import nddata
import dnpLab.dnpImport.vnmrj as vnmrj
import dnpLab.dnpNMR as nmr
import dnpLab as dnp
import numpy as np
import random

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
