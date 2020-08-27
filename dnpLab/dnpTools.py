'''dnpTools
Collection of tools and functions useful to process DNP-NMR data
'''

import inspect
import numpy as np


from .mrProperties import gmrProperties
from .mrProperties import radicalProperties

def mrProperties(nucleus, *args):
    '''Return magnetic resonance property of specified isotope.
    This function is model after the Matlab function gmr written by Mirko Hrovat
    https://www.mathworks.com/matlabcentral/fileexchange/12078-gmr-m-nmr-mri-properties

    Reference: R.K.Harris et. al., Pure and Applied Chemistry, 2001, 73:1795-1818.
    Electron value comes from 1998 CODATA values, http://physics.nist.gov/cuu/Constants .
       or  http://physics.nist.gov/PhysRefData/codata86/codata86.html
       or  http://www.isis.rl.ac.uk/neutronSites/constants.htm
    Xenon gyromagnetic ratio was calculated from 27.661 MHz value from Bruker's web site.

    Args:

        nucleus: String defining the nucleus e.g. 1H, 13C, etc.

        args: If numerical value is given, it is interpreted as the B0 value in Tesla and Larmor frequency is returned. As string the following values are valid:

        gamma: Return Gyromagnetic Ration [radians/T/s]
        spin: Spin number of selected nucleus [1]
        qmom: Quadrupole moment [fm^2} (100 barns)
        natAbundance: Natural abundance [%]
        relSensitivity: Relative sensitiviy with respect to 1H at constant B0
        moment: Magnetic dipole moment in terms of the nuclear magneton, uN, |u|/uN = |gamma|*hbar[I(I + 1)]^1/2/uN , hbar=h/2pi.
        qlw: quadrupolar line-width factor as defined by: Qlw = Q^2(2I + 3)/[I^2(2I � 1)]

    Returns:

        .. code-block:: python

            dnp.dnpTools.mrProperties('1H')
            26.7522128                          # 1H Gyromagnetic Ratio

            dnp.dnpTools.mrProperties('1H', 0.35)
            14902114.17018196                   # 1H Larmor Frequency at 0.35 T

            dnp.dnpTools.mrProperties('2H', 'qmom')
            0.286                               # Nuclear Quadrupole Moment

            value = dnp.dnpTools.mrProperties('6Li', 'natAbundance')
            7.59                                # Natural Abundance

            value = dnp.dnpTools.mrProperties('6Li', 'relSensitivity')
            0.000645                            # Relative sensitivity


    '''

    if isinstance(nucleus, str):
        if nucleus in gmrProperties:
            gmr = gmrProperties.get(nucleus)[1]
        else:
            print("Isotope doesn't exist in list")
            return
    else:            
        print("ERROR: String expected")            

    if len(args) == 0:
        return gmr

    elif len(args) == 1:

        if isinstance(args[0], str):
            if args[0] == 'gamma':
                return gmrProperties.get(nucleus)[1]

            if args[0] == 'spin':
                return gmrProperties.get(nucleus)[0]

            elif args[0] == 'qmom':
                return gmrProperties.get(nucleus)[2]

            elif args[0] == 'natAbundance':
                return gmrProperties.get(nucleus)[3]

            elif args[0] == 'relSensitivity':
                return gmrProperties.get(nucleus)[4]

            elif args[0] == 'moment':
                return gmrProperties.get(nucleus)[5]

            elif args[0] == 'qlw':
                return gmrProperties.get(nucleus)[6]

            else:
                print("Keyword not recognize")

        else:
            vLarmor = args[0] * gmr * 1e7 / 2 / np.pi
            return vLarmor

    elif len(args) > 1:
        print("Too many input arguments")


def getRadicalProperties(name):
    '''Return properties of different radicals

    Args:

    Returns:

    '''

    if isinstance(name, str):
        if name in radicalProperties:
            giso = radicalProperties.get(name)[0]
        else:
            print("Radical doesn't exist in list")
            return
    else:            
        print("ERROR: String expected")            
    
    return giso
