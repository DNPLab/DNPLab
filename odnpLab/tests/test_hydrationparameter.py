#!/usr/bin/env python
# coding: utf-8

import unittest
from odnpLab.hydration.parameter import HydrationParameter


class TestHydrationParameter(unittest.TestCase):

    def test_default_values(self):
        hp = HydrationParameter()
        self.assertEqual(hp.smaxMod, 'tethered')
        self.assertEqual(hp.dH2O, 2.3e-9)

    def test_should_value_error_smaxmod(self):
        hp = HydrationParameter()
        with self.assertRaises(ValueError):
            hp.smaxMod = 'gear2'

    def test_dict_like_get(self):
        hp = HydrationParameter()
        hp.field = 380
        self.assertEqual(hp['field'], 380)
        self.assertEqual(hp['smaxMod'], hp.smaxMod)

    def test_dict_like_set(self):
        hp = HydrationParameter()
        hp['field'] = 400
        self.assertEqual(hp.field, 400)

        hp['smaxMod'] = 'free'
        self.assertEqual(hp.smaxMod, 'free')

        with self.assertRaises(ValueError):
            hp['smaxMod'] = 'gear2'
