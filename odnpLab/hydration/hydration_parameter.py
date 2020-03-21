#!/usr/bin/env python
# coding: utf-8
from odnpLab.parameter import Parameter


class HydrationParameter(Parameter):
    """Hydration Parameters Getting and Setting"""

    def __init__(self):
        """"""
        super().__init__()

        self.field = 380
        # static magnetic field in mT, needed to find omega_e and _H

        self.slC = 100
        # (Eq. 1-2) unit is M, spin label concentration for scaling relaxations to
        # get "relaxivities"

        self.__smaxMod = 'tethered'  # either 'tethered' or 'free'

        self.ksig_bulk = 95.4  # unit is s^-1 M^-1
        # The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
        # Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
        # 2015, 137, 12013−12023. Figure 3 caption

        self.T10 = 1.33  # this is the T1 with spin label but at 0 mw power, unit is sec
        self.T100 = 2.5  # this is the T1 without spin label and without mw power, unit is sec

        self.tcorr_bulk = 54  # (section 2.5), "corrected" bulk tcorr, unit is ps

        self.dH2O = 2.3e-9
        # (Eq. 19-20) bulk water diffusivity, unit is d^2/s where d is distance in
        # meters. *didnt use m to avoid confusion with mass

        self.dSL = 4.1e-10
        # (Eq. 19-20) spin label diffusivity, unit is d^2/s where d is distance in
        # meters. *didnt use m to avoid confusion with mass

        self.k_low_bulk = 366  # unit is s^-1 M^-1
        # The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
        # Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
        # 2015, 137, 12013−12023. Figure 3 caption

        # fit option, either linear or 2nd order  TODO: insert ref
        self.__fitOpt = 'linear'

    @property
    def fitopt(self):
        return self.__fitOpt

    @fitopt.setter
    def fitopt(self, value: str):
        if value in ['2ord', 'linear']:
            self.__fitOpt = value
        else:
            raise ValueError('fitopt should be either `linear` or `2ord`')

    @property
    def smaxMod(self):
        return self.__smaxMod

    @smaxMod.setter
    def smaxMod(self, value: str):
        if value == 'tethered':
            self.__smaxMod = value
        elif value == 'free':
            self.__smaxMod = value
        else:
            raise ValueError('smaxMod should be either `tethered` or `free`')

    def __getitem__(self, key):
        if key in ['smaxMod']:
            return self.smaxMod
        return self.__dict__[key]

    def __setitem__(self, key, value):
        if key in ['smaxMod']:
            self.smaxMod = value
        else:
            self.__dict__[key] = value

    def print(self):
        # TODO: implement this
        raise NotImplementedError

