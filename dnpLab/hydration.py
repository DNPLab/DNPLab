""" Hydration module

This module calculates hydration related quantities using processed ODNP data.

"""

import numpy as np
from scipy import interpolate
from scipy import optimize
from dnpLab.utils import AttrDict, Parameter


class FitError(Exception):
    """Exception of Failed Fitting"""
    pass


class HydrationParameter(Parameter):
    """Hydration Parameters

        Franck, JM, et. al.; "Anomalously Rapid Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc. 2015, 137, 12013−12023.
    """

    field = None
    """float: Static magnetic field in mT, needed to find omega_e and _H"""

    spin_C = None
    """float: (Eq. 1-2) unit is microM, spin label concentration for scaling
    relaxations to get "relaxivities" """

    __smax_model = 'tethered'  # either 'tethered' or 'free'
    """str: Method used to determine smax"""

    T10 = None
    """float: T1 with spin label but at 0 mw E_power, unit is sec"""

    T100 = None
    """float: T1 without spin label and without mw E_power, unit is sec"""

    ksigma_bulk = 95.4
    """float: unit is s^-1 M^-1 (Figure 3 caption)"""

    klow_bulk = 366
    """float: unit is s^-1 M^-1 (Figure 3 caption)"""

    tcorr_bulk = 54
    """float: Corrected bulk tcorr, unit is ps, (section 2.5)"""

    D_H2O = 2.3e-9
    """float: (Eq. 19-20) bulk water diffusivity, unit is d^2/s where d is 
    distance in meters."""

    D_SL = 4.1e-10
    """float: (Eq. 19-20) spin label diffusivity, unit is d^2/s where d is 
    distance in meters."""

    __t1_interp_method = 'second_order'
    """str: Method used to interpolate T1, either linear or 'second_order'"""

    def __init__(self):
        """"""
        super().__init__()
        # Fixme: Remove the default field, spin concentration etc. Set them as required arguments
        self.field = 348.5
        self.spin_C = 100
        self.T10 = 1.5
        self.T100 = 2.5

    @property
    def t1_interp_method(self):
        """str: Method used to interpolate T1, either `linear` or `second_order`"""
        return self.__t1_interp_method

    @t1_interp_method.setter
    def t1_interp_method(self, value: str):
        if value in ['second_order', 'linear']:
            self.__t1_interp_method = value
        else:
            raise ValueError('t1_interp_method should be either `linear` or `second_order`')
            
    @property
    def smax_model(self):
        """str: Method used to determine smax. Either `tethered` or `free`"""
        return self.__smax_model

    @smax_model.setter
    def smax_model(self, value: str):
    
        if value == 'tethered':
            self.__smax_model = value
        elif value == 'free':
            self.__smax_model = value
        else:
            raise ValueError('smax_model should be either `tethered` or `free`')
    
    # These two function enable dictionary-like getting and setting properties.
    def __getitem__(self, key):
        if key in ['smax_model']:
            return self.smax_model
        elif key in ['t1_interp_method']:
            return self.t1_interp_method
        return self.__dict__[key]

    def __setitem__(self, key, value):
        if key in ['smax_model']:
            self.smax_model = value
        elif key in ['t1_interp_method']:
            self.t1_interp_method = value
        else:
            self.__dict__[key] = value
    

class HydrationResults(AttrDict):
    """Class for handling hydration related quantities

    Attributes:
        uncorrected_Ep (numpy.array)    : Fit of Ep array
        interpolated_T1 (numpy.array)   : T1 values interpolated on E_power,
        ksigma_array (numpy.array)      :
            numpy array that is the result of ~(1-E) / [ (constants*T1) ],
            used in ksigma(E_power) fit,
        ksigma_fit (numpy.array)        : ksig_fit,
        k_sigma (float)                 : k_sigma,
        ksigma_error (float)            : ksigma_error,
        ksigma_bulk_ratio (float)       : k_sigma/ksigma_bulk,
        krho (float)                    : krho,
        klow (float)                    : klow,
        klow_bulk_ratio (float)         : klow / klow_bulk,
        coupling_factor (float)         : coupling_factor,
        tcorr (float)                   : tcorr,
        tcorr_bulk_ratio (float)        : tcorr / tcorr_bulk,
        D_local (float)                 : D_local

    """
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.uncorrected_Ep = None
        self.interpolated_T1 = None
        self.ksigma_array = None
        self.ksigma_fit = None
        self.k_sigma = None
        self.ksigma_error = None
        self.ksigma_bulk_ratio = None
        self.krho = None
        self.klow = None
        self.klow_bulk_ratio = None
        self.coupling_factor = None
        self.tcorr = None
        self.tcorr_bulk_ratio = None
        self.D_local = None
        self.update(*args, **kwargs)

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()


class HydrationCalculator:
    """Hydration Results Calculator

    Attributes:
        T1 (numpy.array): T1 array. Unit: second.
        T1_power (numpy.array): E_power in Watt unit, same length as T1.
        E (numpy.array): Enhancements.
        E_power (numpy.array): E_power in Watt unit, same length as E.
        hp (HydrationParameter): Parameters for calculation, including default
            values.
        results (HydrationResults): Hydration results.

    """
    def __init__(self, T1: np.array, T1_power: np.array, E: np.array,
                 E_power: np.array, hp: HydrationParameter):
        """Class Init

        Args:
            T1 (numpy.array): T1 array. Unit: second.
            T1_power (numpy.array): E_power in Watt unit, same length as T1.
            E (numpy.array): Enhancements.
            E_power (numpy.array): E_power in Watt unit, same length as E.
            hp (HydrationParameter): Parameters for calculation, including
                default values.
        """
        super().__init__()

        if T1.size != T1_power.size:
            raise ValueError('T1 and T1_power must have same length')
        if E.size != E_power.size:
            raise ValueError('E and E_power must have same length')

        self.T1, self.T1_power, self.E, self.E_power = T1, T1_power, E, E_power

        self.hp = hp

        self.results = HydrationResults()

    def run(self):
        """Run calculator"""
        interpolated_T1 = self.interpolate_T1(self.E_power, self.T1_power, self.T1)
        results = self.__calculateODNP(self.E_power, self.E, interpolated_T1)
        self.results = results

    def interpolate_T1(self, E_power: np.array, T1_power: np.array, T1: np.array):
        """Returns the one-dimensional piecewise interpolant to a function with
        given discrete data points (T1_power, T1), evaluated at E_power.

        Points outside the data range will be extrapolated

        Args:
            E_power: The x-coordinates at which to evaluate.
            T1_power: The x-coordinates of the data points, must be increasing.
                Otherwise, T1_power is internally sorted.
            T1: The y-coordinates of the data points, same length as T1_power.

        Returns:
            interplatedT1 (np.array): The evaluated values, same shape as E_power.

        """
        T10, T100 = self.hp.T10, self.hp.T100
        spin_C = self.hp.spin_C / 1e6
        
        t1_interp_method = self.hp.t1_interp_method

        # 2nd order fit, Franck and Han MIE (Eq. 22) and (Eq. 23)
        if t1_interp_method == 'second_order':

            delta_T1_water = T1[-1] - T1[0]
            T1_water = T100
            macro_C = spin_C

            kHH = (1. / T10 - 1. / T1_water) / macro_C
            krp = ((1. / T1) - (1. / (T1_water + delta_T1_water * T1_power)) - (kHH * (macro_C))) / (spin_C)

            p = np.polyfit(T1_power, krp, 2)
            T1_fit_2order = np.polyval(p, E_power)

            interp_T1 = 1./(((spin_C) * T1_fit_2order) + (1. / (T1_water + delta_T1_water * E_power)) + (kHH * (macro_C)))

        # linear fit, Franck et al. PNMRS (Eq. 39)
        elif t1_interp_method == 'linear':

            linear_t1=1./((1. / T1) - (1. / T10) + (1. / T100))

            p = np.polyfit(T1_power, linear_t1, 1)
            T1_fit_linear = np.polyval(p, E_power)

            interp_T1=T1_fit_linear/(1.+(T1_fit_linear/T10)-(T1_fit_linear/T100))

        else:
            raise Exception('Error in T1 interpolation')
        
        return interp_T1

    def __calculateODNP(self, power:np.array, Ep:np.array, T1p:np.array):
        """Returns a HydrationResults object that contains all calculated ODNP values

        Following: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
        equations are labeled (#) for where they appear in the paper, sections are specified in some cases

        Args:
            power (numpy.array): Sequence of powers in Watts. Will be internally
                sorted to be ascending.
            Ep (numpy.array): Array of enhancements. Must be same length as E_power
            T1p (numpy.array): Array of T1. Must be same length as E_power.
        """
        # field and spin label concentration are defined in Hydration Parameter
        field = self.hp.field
        spin_C = self.hp.spin_C / 1e6
        
        T10 = self.hp.T10  # this is the T1 with spin label but at 0 mw E_power, unit is sec
        T100 = self.hp.T100  # this is the T1 without spin label and without mw E_power, unit is sec

        if self.hp.smax_model == 'tethered':
            # Option 1, tether spin label
            s_max = 1  # (section 2.2) maximal saturation factor

        elif self.hp.smax_model == 'free':
            # Option 2, free spin probe
            s_max = 1-(2/(3+(3*(spin_C*198.7)))) # from:
            # M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. &
            # J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

        omega_e = (1.76085963023e-1) * (field / 1000)
        # gamma_e in 1/ps for the tcorr unit, then correct by field in T.
        # gamma_e is from NIST. The field cancels in the following omega_ratio but you
        # need these individually for the spectral density functions later.

        omega_H = (2.6752218744e-4) * (field / 1000)
        # gamma_H in 1/ps for the tcorr unit, then correct by field in T.
        # gamma_H is from NIST. The field cancels in the following omega_ratio but you
        # need these individually for the spectral density functions later.

        omega_ratio = ((omega_e / (2 * np.pi)) / (omega_H / (2 * np.pi)))
        # (Eq. 4-6) ratio of omega_e and omega_H, divide by (2*pi) to get angular
        # frequency units in order to correspond to S_0/I_0, this is also ~= to the
        # ratio of the resonance frequencies for the experiment, i.e. MW freq/RF freq

        ksigma_array = ((1 - Ep) / (spin_C * omega_ratio * T1p))
        # (Eq. 41) this calculates the array of k_sigma*s(p) from the enhancement array,
        # dividing by the T1 array for the "corrected" analysis

        popt, pcov = self.get_ksigma(ksigma_array, power)
        # fit to the right side of Eq. 42 to get (k_sigma*smax) and half of the E_power at s_max, called p_12 here
        ksigma_smax = popt[0]
        p_12 = popt[1]
        ksigma_std = np.sqrt(np.diag(pcov))
        ksigma_error = ksigma_std[0] / s_max

        ksigma_fit = (ksigma_smax * power) / (p_12 + power)
        # (Eq. 42) calculate the "corrected" k_sigma*s(p) array using the fit parameters,
        # this can be used to plot over the data ksigma_array array to assess the quality of fit.
        # This would correspond to the corrected curves in Figure 9
        
        k_sigma = ksigma_smax / s_max
        
        # ksig_uncorr = ((1 - Ep) / (spin_C * omega_ratio * T10)) / s_max
        # (Eq. 44) the "uncorrected" model, this can also be plotted with the corrected
        # curve to determine the severity of heating effects, as in Figure 9.
        # Notice the division by T10 instead of the T1 array
        #################

        ksigma_bulk = self.hp.ksigma_bulk  # unit is s^-1 M^-1
        # "Anomalously Rapid Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc. 2015, 137, 12013−12023. Figure 3 caption

        krho = ((1/T10) - (1/T100)) / spin_C  # (Eq. 36) "self" relaxivity, unit is s^-1 M^-1

        coupling_factor = k_sigma / krho  # (Eq. 3) this is the coupling factor, unitless

        tcorr = self.get_tcorr(coupling_factor, omega_e, omega_H)
        # (Eq. 21-23) this calls the fit to the spectral density functions. The fit
        # optimizes the value of tcorr in the calculation of coupling_factor, the correct tcorr
        # is the one for which the calculation of coupling_factor from the spectral density
        # functions matches the coupling_factor found experimentally. tcorr unit is ps

        tcorr_bulk = self.hp.tcorr_bulk  # (section 2.5), "corrected" bulk tcorr, unit is ps

        D_H2O = self.hp.D_H2O
        # (Eq. 19-20) bulk water diffusivity, unit is d^2/s where d is distance in
        # meters. *didnt use m to avoid confusion with mass

        D_SL = self.hp.D_SL
        # (Eq. 19-20) spin label diffusivity, unit is d^2/s where d is distance in
        # meters. *didnt use m to avoid confusion with mass

        D_local = (tcorr_bulk / tcorr) * (D_H2O + D_SL)
        # (Eq. 19-20) local diffusivity, i.e. diffusivity of the water near the spin label

        ############################################################################
        # This is defined in its most compact form in:
        # Frank, JM and Han, SI;  Chapter Five - Overhauser Dynamic Nuclear Polarization
        # for the Study of Hydration Dynamics, Explained. Methods in Enzymology, Volume 615, 2019
        #
        # But also explained well in:
        # Franck, JM, et. al.; "Anomalously Rapid Hydration Water Diffusion Dynamics
        # Near DNA Surfaces" J. Am. Chem. Soc. 2015, 137, 12013−12023.

        klow = ((5 * krho) - (7 * k_sigma)) / 3
        # section 6, (Eq. 13). this describes the relatively slowly diffusing water
        # near the spin label, sometimes called "bound" water
        ############################################################################

        klow_bulk = self.hp.klow_bulk  # unit is s^-1 M^-1
        # "Anomalously Rapid Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc. 2015, 137, 12013−12023. Figure 3 caption
        
        results = self.get_uncorrected_xi(Ep, power, T10, T100, omega_ratio, s_max)
        xi_unc = results.x[0]
        p_12_unc = results.x[1]
        
        #placeholder
        """
        J = results.jac
        cov = np.linalg.inv(J.T.dot(J))
        results_std = np.sqrt(np.diagonal(cov))
        xi_unc_error = results_std[0]
        """
        
        uncorrected_Ep = 1-((xi_unc*(1-(T10/T100))*omega_ratio)*((power*s_max)/(p_12_unc+power)))
        
        return HydrationResults({
            'uncorrected_Ep'    : uncorrected_Ep,
            'interpolated_T1'   : T1p,
            'ksigma_array'      : ksigma_array,
            'ksigma_fit'        : ksigma_fit,
            'k_sigma'           : k_sigma,
            'ksigma_error'      : ksigma_error,
            'ksigma_bulk_ratio' : k_sigma/ksigma_bulk,
            'krho'              : krho,
            'klow'              : klow,
            'klow_bulk_ratio'   : klow / klow_bulk,
            'coupling_factor'   : coupling_factor,
            'tcorr'             : tcorr,
            'tcorr_bulk_ratio'  : tcorr / tcorr_bulk,
            'D_local'           : D_local
        })

    @staticmethod
    def get_tcorr(coupling_factor: float, omega_e: float, omega_H: float):
        """Returns correlation time tcorr in pico second

        Args:
            coupling_factor (float):
            omega_e (float):
            omega_H (float):

        Returns:
            float: correlation time in pico second

        Raises:
            FitError: If no available root is found.

        """

        def calc_xi(tcorr: float, omega_e: float, omega_H: float):
            '''Returns coupling_factor for any given tcorr

            Args:
                tcorr (float):
                omega_e (float):
                omega_H (float):

            Returns:
                coupling_factor (float):

            '''

            # Using Franck et al. PNMRS (2013)

            zdiff = np.sqrt(1j * (omega_e - omega_H) * tcorr)
            zsum  = np.sqrt(1j * (omega_e + omega_H) * tcorr)
            zH    = np.sqrt(1j * omega_H * tcorr)

            # (Eq. 2)
            Jdiff = (1 + (zdiff / 4)) / (1 + zdiff + ((4 * (zdiff**2)) / 9) + ((zdiff**3) / 9))

            Jsum  = (1 + (zsum / 4)) / (1 + zsum + ((4 * (zsum**2)) / 9) + ((zsum**3) / 9))

            JH    = (1 + (zH / 4)) / (1 + zH + ((4 * (zH**2)) / 9) + ((zH**3) / 9))

            # (Eq. 23) calculation of coupling_factor from the spectral density functions
            xi_tcorr = ((6 * np.real(Jdiff)) - np.real(Jsum)) / ((6 * np.real(Jdiff)) + (3 * np.real(JH)) + np.real(Jsum))

            return xi_tcorr

        # root finding
        # see https://docs.scipy.org/doc/scipy/reference/optimize.html
        result = optimize.root_scalar(
            lambda tcorr: calc_xi(tcorr, omega_e=omega_e, omega_H=omega_H) - coupling_factor,
            method='brentq',
            bracket=[1, 1e5])

        if not result.converged:
            raise FitError("Could not find tcorr")
        return result.root

    @staticmethod
    def get_ksigma(ksig_sp: np.array, power: np.array):
        """Get k_sigma and E_power at half max of ksig

        Args:
            ksig (numpy.array): Array of k_sigma.
            power (numpy.array): Array of E_power.

        Returns:
            popt: fit results
            pcov: covariance matrix

        Asserts:
            k_sigma (popt[0]) is greater than zero

        """

        def calc_ksigma(power: np.array, ksigma_smax: float, p_12: float):
            """Function to calcualte ksig array for any given ksigma and p_12

            Args:
                power (numpy.array): Array of E_power.

            Returns:
                fit to ksig array

            """

            # Again using: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

            # Right side of Eq. 42. This function should fit to ksig_sp
            ksig_fit = (ksigma_smax * power) / (p_12 + power)

            return ksig_fit

        # curve fitting
        # see https://docs.scipy.org/doc/scipy/reference/optimize.html
        popt, pcov = optimize.curve_fit(calc_ksigma, power, ksig_sp,
                                         p0=[50, (max(power) / 2)], method='lm')

        assert popt[0] > 0, 'Unexpected ksigma value: %d < 0' % popt[0]
        return popt, pcov

    @staticmethod
    def get_uncorrected_xi(Ep: np.array, power: np.array, T10: float, T100: float, wRatio: float, s_max: float):
        """Get coupling_factor and E_power at half saturation

        Args:
            Ep (numpy.array): Array of enhancements.
            power (numpy.array): Array of E_power.
            T10 (float): T10
            T100 (float): T100
            wRatio (float): ratio of electron & proton Larmor frequencies
            s_max (float): maximal saturation factor

        Returns:
            A tuple of float (coupling_factor, p_12).

        Raises:
            FitError: If least square fitting is not succeed.

        """

        def residual(x, Ep: np.array, power: np.array, T10: float, T100: float, wRatio: float, s_max: float):
            """Residual function for Ep for any given xi and p_12

            Args:
                x (tuple): length of 2
                Ep (numpy.array): Array of enhancements.
                power (numpy.array): Array of E_power.
                T10 (float): T10
                T100 (float): T100
                wRatio (float): ratio of electron & proton Larmor frequencies
                s_max (float): maximal saturation factor

            Returns:
                Residuals.

            """
            xi_unc, p_12_unc = x[0], x[1]

            # Again using: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

            # Right side of Eq. 42. This function should fit to ksig_sp
            Ep_fit = 1-((xi_unc*(1-(T10/T100))*wRatio)*((power*s_max)/(p_12_unc+power)))

            return Ep - Ep_fit

        # least-squares fitting.
        # see https://docs.scipy.org/doc/scipy/reference/optimize.html
        results = optimize.least_squares(fun=residual,
                                         x0=[0.5, (max(power) / 2)],
                                         args=(Ep, power, T10, T100, wRatio, s_max),
                                         jac='2-point', method='lm')
        if not results.success:
            raise FitError('Could not fit Ep')
        assert results.x[0] > 0, 'Unexpected coupling_factor value: %d < 0' % results.x[0]
        return results


def hydration(ws):

    # Hydration required data
    hyd = ws['hydration']
    T1, T1_power, E, E_power = [hyd[k] for k in ['T1', 'T1_power', 'E', 'E_power']]

    # Create hydration parameter
    hp = HydrationParameter()
    hp.T10, hp.T100, hp.spin_C, hp.field, hp.smax_model, hp.t1_interp_method = \
        [hyd[k] for k in ['T10', 'T100', 'spin_C', 'field', 'smax_model', 't1_interp_method']]

    # compute
    hc = HydrationCalculator(T1, T1_power, E, E_power, hp)
    hc.run()

    return
