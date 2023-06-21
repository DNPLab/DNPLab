
from __future__ import division
import numpy as _np


base_units = [
    's', # second
    'm', # meter
    'kg', # kilogram
    'A', # Amps
    'K', # Kelvin
    'mol', # mol
    'cd', # candela
    ]

derived_units = [
    'rad', # radians
    'sr', # steradian
    'Hz', # frequency
    'N', # newton
    'Pa', # Pascal
    'J', # Joule
    'W', # Watt
    'C', # Coulomb
    'V', # Volt
    'F', # Farad
    'ohm', # Ohm
    '\u03A9', # Ohm
    'S', # Siemens
    'Wb', # Weber
    'T', # Tesla
    'H', # Henry
    'C', # Degrees C
    'lm', # Lumen
    'lx', # Lux
    'Bq', # becquerel
    'Gy', # Gray
    'Sv', # Sievert
    'kat', # katal
    ]

prefix = {
    'q': 1e-30, # quecto
    'r': 1e-27, # ronto
    'y': 1e-24, # yocto
    'z': 1e-21, # zepto
    'a': 1e-18, # atto
    'f': 1e-15, # femto
    'p': 1e-12, # pico
    'n': 1e-9, # nano
    'u': 1e-6, # micro
    '\u03BC':1e-6, # unicode for mu, micro
    'm': 1e-3, # mili
    '': 1e0, # None
    'k': 1e3, # kilo
    'M': 1e6, # Mega
    'G': 1e9, # Giga
    'T': 1e12, # Tera
    'P': 1e15, # Peta
    'E': 1e18, # Exa
    'Z': 1e21, # Zetta
    'Y': 1e24, # Yotta
    'R': 1e27, # Ronna
    'Q': 1e30, # Quetta
          }

unicode_maps = {
    'ohm' : '\u03A9',
    'u' : '\u03BC'
}

special_units = ['ppm', 'kg', 'gauss']

class Unit(object):
    def __init__(self, symbol = None, prefix = None, autodetect_prefix = True, convert_to_base = True):
        """Units class for storing units of Values and Coords
        """

        self.dimensionless = False

        self.unit
        if symbol == None:
            self.symbol = ''
            self.dimensionless = True
    

    @property
    def symbol(self):
        return self._symbol
    
    @symbol.setter
    def symbol(self, symbol):
        self._symbol = symbol
    
    @symbol.getter
    def symbol(self, symbol):
        return self._symbol
