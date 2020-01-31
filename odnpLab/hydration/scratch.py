# hydration function(s)
# goal: take odnpImport and output tcorr and many observables

import numpy as np
from scipy import interpolate
from scipy import optimize


def getT1p(T1: np.array, power: np.array):
    """
    returns a function to calculate T1 at arbitrary power

    :param T1: np.array of T1
    :param power: np.array of power in Watt unit, same length as T1
    :return
        a function that takes power and returns T1.
        points inside the data range will be linearly interpolated
        points outside the data range will be extrapolated
    """
    assert len(T1) == len(power)
    return interpolate.interp1d(power, T1,
                                kind='linear',
                                fill_value='extrapolate')


def calcODNP(Ep: np.array, T1p: callable):
    """
    returns all calculated values

    :param E: np.array of enhancements
    :param power: np.array of power in Watt unit, same length as E
    :param t1p_fun: callable function that returns T1 at any power
    :return:
        float k_sigma
    """
    # Following: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    # equations are labeled (#) for where they appear in the paper, sections are specified in some cases

    field = 348.5  # static magnetic field in mT, needed to find omega_e and _H

    slC = 200e-6
    # (Eq. 1-2) unit is M, spin label concentration for scaling relaxations to
    # get "relaxivities"

    s_max = 1  # (section 2.2) maximal saturation factor
    # for s_max eventually we can specify if bulk or not and use 1 or,
    # s_max = 1-(2/(3+(3*(slC*1e-6*198.7)))), from:
    # M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. &
    # J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

    omega_e = 1.76085963023e11 * (field / 1000)
    # gamma_e is from NIST, field is converted to T. The field cancels jin the
    # following ratio but you need these individually for the spectral density
    # function later. unit is Hz

    omega_H = 2.6752218744e8 * (field / 1000)
    # gamma_H is from NIST, field is converted to T. The field cancels in the
    # following ratio but you need these individually for the spectral density
    # function later. unit is Hz

    wRatio = omega_e / omega_H  # (Eq. 4-6) ratio of omega_e and omega_H, ~= S_0/I_0

    # Ep will be the series of enhancements
    # T1p will be the series of T1s after interpolated to match the enhancement
    # series. Power not needed anymore, it is only used to line up and
    # interpolate the T1s to the Enhancements

    ksig_smax = (1 - Ep) / (slC * wRatio * T1p)
    # (Eq. 42) this calculates the series of k_sigma*s(p) which approximates to
    # k_sigma*s_max

    k_sigma = max(ksig_smax) / s_max
    # (Eq. 43) this takes the maximum of the k_sigma*s(p) series and divides by
    # the s_max to isolate k_sigma, unit is s^-1 M^-1

    ksig_bulk = 95.4  # unit is s^-1 M^-1
    # The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
    # Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
    # 2015, 137, 12013−12023. Figure 3 caption

    T10 = 1.33  # this is the T1 with spin label but at 0 mw power, unit is sec
    T100 = 2.5  # this is the T1 without spin label and without mw power, unit is sec
    k_rho = ((1/T10) - (1/T100)) / slC  # (Eq. 36) "self" relaxivity, unit is s^-1 M^-1

    ksi = k_sigma / k_rho  # (3) this is the coupling factor, unitless

    tcorr = getTcorr(ksi, omega_e, omega_H)
    # (Eq. 21-23) this calls the fit to the spectral density functions. The fit
    # optimizes the value of tcorr in the calculation of ksi, the correct tcorr
    # is the one for which the calculation of ksi from the spectral density
    # functions matches the ksi found experimentally. tcorr unit is ps

    tcorr_bulk = 54  # (section 2.5), "corrected" bulk tcorr, unit is ps

    dH2O = 2.3e-9
    # (Eq. 19-20) bulk water diffusivity, unit is d^2/s where d is distance in
    # meters. *didnt use m to avoid confusion with mass

    dSL = 4.1e-10
    # (Eq. 19-20) spin label diffusivity, unit is d^2/s where d is distance in
    # meters. *didnt use m to avoid confusion with mass

    dLocal = (tcorr_bulk / tcorr) * (dH2O + dSL)
    # (Eq. 19-20) local diffusivity, i.e. diffusivity of the water near the spin label

    ############################################################################
    # This is defined in its most compact form in:
    # Frank, JM and Han, SI;  Chapter Five - Overhauser Dynamic Nuclear Polarization
    # for the Study of Hydration Dynamics, Explained. Methods in Enzymology, Volume 615, 2019
    #
    # But also explained well in:
    # Franck, JM, et. al.; "Anomalously Rapid Hydration Water Diffusion Dynamics
    # Near DNA Surfaces" J. Am. Chem. Soc. 2015, 137, 12013−12023.

    k_low = ((5 * k_rho) - (7 * k_sigma)) / 3
    # section 6, (Eq. 13). this describes the relatively slowly diffusing water
    # near the spin label, sometimes called "bound" water
    ############################################################################

    klow_bulk = 366  # unit is s^-1 M^-1
    # The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
    # Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
    # 2015, 137, 12013−12023. Figure 3 caption

    # these quantities are often used
    print 'k_sigma = ', k_sigma
    print '[k_sig/k_sig,bulk]^-1 = ', 1/(k_sigma/ksig_bulk)
    print 'k_low = ', k_low
    print 'k_low/klow_bulk = ', k_low/klow_bulk
    print 'tcorr = ', tcorr
    print 'tcorr retardation = ', tcorr/tcorr_bulk
    print 'coupling factor = ', ksi
    print 'Diffusivity = ', dLocal
    
    
    return k_sigma k_rho ksi tcorr dLocal k_low


def getTcorr(ksi: float):
    """
    returns correlation time tcorr

    :param ksi: float of epsilon
    :return:
        float tcorr
    """

    def get_ksi(tcorr: float):
        """
        returns ksi for any given tcorr
        :param tcorr: float
        :return:
            float ksi
        """

        # Again using: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

        # (Eq. 22), difference, sum and H terms for "z"
        zdiff = np.sqrt(1j * (omega_e - omega_H) * tcorr)
        zsum  = np.sqrt(1j * (omega_e + omega_H) * tcorr)
        zH    = np.sqrt(1j * omega_H * tcorr)

        # (Eq. 21) the three forms of the FFHS spectral density function corresponding to the three "z" terms above
        Jdiff = np.real((1 + (zdiff / 4)) / (1 + zdiff + ((4 * (zdiff ** 2)) / 9) + ((zdiff ** 3) / 9)))
        Jsum  = np.real((1 + (zsum  / 4)) / (1 + zsum  + ((4 * (zsum ** 2)) / 9)  + ((zsum ** 3) / 9)))
        JH    = np.real((1 + (zH    / 4)) / (1 + zH    + ((4 * (zH ** 2)) / 9)    + ((zH ** 3) / 9)))
        
        # (Eq. 23) calculation of ksi from the spectral density functions
        ksi_tcorr = ((6 * Jdiff) - Jsum) / ((6 * Jdiff) + (3 * JH) + Jsum)

        return ksi_tcorr

    # root finding
    # see https://docs.scipy.org/doc/scipy/reference/optimize.html
    results = optimize.root_scalar(
        lambda tcorr: ((get_ksi(tcorr) - ksi) ** 2),
        method='newton',
        x0=500)
    assert results.converged
    return results.root
