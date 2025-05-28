import unittest
import os
import dnplab as dnp
from numpy.testing import assert_array_equal
import numpy as np

import logging
import sys
import pathlib
import warnings

logger = logging.getLogger(__name__)


class dnp_processing_tester(unittest.TestCase):

    def setUp(self):
        x = np.linspace(0,2,500)
        y = 2 * x + 1
        # x^2 + x + C
        self.data = dnp.DNPData(y, dims = ['x0'], coords = [x] )
        self.fullIntegral = 6
        self.halfIntegral = 2

    def test_000_integration_regions(self):

        # test value and that regions accepts a single tuple
        integralFull = dnp.processing.integrate(self.data, 'x0')
        integralHalf = dnp.processing.integrate(self.data, 'x0', regions = (0.,1.) )
        integralHalf2 = dnp.processing.integrate(self.data, 'x0', regions = [(0.,1.)] )
        self.assertTrue (np.isclose( integralFull._values, self.fullIntegral, rtol=1e-3, atol=1e-3) )
        self.assertTrue (np.isclose( integralHalf._values, self.halfIntegral, rtol=1e-2, atol=7e-3) )
        self.assertTrue (np.isclose( integralHalf2._values, self.halfIntegral, rtol=1e-2, atol=7e-3) )

        # test multiple regions
        integral2 = dnp.processing.integrate(self.data, 'x0', regions = [(0.,0.5),(1,2.0)] )
        v1 = integral2._values[0]
        v2 = integral2._values[1]

        self.assertTrue( np.isclose(v1, 0.75, rtol=0.01, atol = 0.001) )
        self.assertTrue( np.isclose(v2, 4, rtol=0.01, atol = 0.001) )

