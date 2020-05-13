""" Hydration module

This module calculates hydration related quantities using processed ODNP data.

"""

import numpy as np
from scipy import interpolate
from scipy import optimize
from odnpLab.parameter import AttrDict, Parameter


class FitError(Exception):
    """Exception of Failed Fitting"""
    pass


class HydrationParameter(Parameter):
    """Hydration Parameters Getting and Setting"""

    def __init__(self):
        """"""
        super().__init__()

        self.field = 348.5
        """float: Static magnetic field in mT, needed to find omega_e and _H"""

        self.slC = 100
        """float: (Eq. 1-2) unit is M, spin label concentration for scaling 
        relaxations to get "relaxivities" """

        self.__smaxMod = 'tethered'  # either 'tethered' or 'free'
        """str: Method used to determine smax"""

        self.ksig_bulk = 95.4
        """float: unit is s^-1 M^-1

        The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
        Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
        2015, 137, 12013−12023. Figure 3 caption
        """

        self.T10 = 1.5
        """float: T1 with spin label but at 0 mw power, unit is sec"""

        self.T100 = 2.5
        """float: T1 without spin label and without mw power, unit is sec"""

        self.tcorr_bulk = 54
        """float: (section 2.5), "corrected" bulk tcorr, unit is ps"""

        self.dH2O = 2.3e-9
        """float: (Eq. 19-20) bulk water diffusivity, unit is d^2/s where d is 
        distance in meters. *didnt use m to avoid confusion with mass"""

        self.dSL = 4.1e-10
        """float: (Eq. 19-20) spin label diffusivity, unit is d^2/s where d is distance in
        meters. *didnt use m to avoid confusion with mass"""

        self.k_low_bulk = 366
        """float: unit is s^-1 M^-1

        The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
        Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
        2015, 137, 12013−12023. Figure 3 caption"""

        self.__t1InterpMethod = '2ord'
        """str: Method used to interpolate T1, either linear or 2nd order"""

    @property
    def t1InterpMethod(self):
        """str: Method used to interpolate T1, either linear or 2nd order"""
        return self.__t1InterpMethod

    @t1InterpMethod.setter
    def t1InterpMethod(self, value: str):
        if value in ['2ord', 'linear']:
            self.__t1InterpMethod = value
        else:
            raise ValueError('t1InterpMethod should be either `linear` or `2ord`')

    @property
    def smaxMod(self):
        """str: Method used to determine smax"""
        return self.__smaxMod

    @smaxMod.setter
    def smaxMod(self, value: str):
        if value == 'tethered':
            self.__smaxMod = value
        elif value == 'free':
            self.__smaxMod = value
        else:
            raise ValueError('smaxMod should be either `tethered` or `free`')

    # These two function enable dictionary-like getting and setting properties.
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


class HydrationResults(AttrDict):
    """Class for handling hydration related quantities

    Attributes:
        Ep_unc (numpy.array): Fit of Ep array.
        T1interp (numpy.array): Interpolated T1 values on E_power.
        k_sigma_array (numpy.array)     : ksig_sp / sp,
        k_sigma_fit (numpy.array)       : ksigsp_fit,
        k_sigma (float)                 : k_sigma,
        ksigma_kbulk_invratio (float)   : 1/(k_sigma/ksig_bulk),
        k_rho (float)                   : k_rho,
        k_low (float)                   : k_low,
        klow_klow_bulk_ratio (float)    : k_low / k_low_bulk,
        ksi (float)                     : ksi,
        tcorr (float)                   : tcorr,
        tcorr_tcorr_bulk_ratio (float)  : tcorr / tcorr_bulk,
        dLocal (float)                  : dLocal

    """
    def __init__(self, *args, **kwargs):
        super().__init__()
        self.T1interp = None
        self.k_sigma_array = None
        self.k_sigma = None
        self.ksigma_kbulk_invratio = None
        self.k_rho = None
        self.k_low = None
        self.klow_klow_bulk_ratio = None
        self.ksi = None
        self.tcorr = None
        self.tcorr_tcorr_bulk_ratio = None
        self.dLocal = None
        self.update(*args, **kwargs)


class HydrationCalculator:
    """Hydration Results Calculator

    Attributes:
        T1 (numpy.array): T1 array. Unit: second.
        T1_power (numpy.array): power in Watt unit, same length as T1.
        E (numpy.array): Enhancements.
        E_power (numpy.array): power in Watt unit, same length as E.
        hp (HydrationParameter): Parameters for calculation, including default
            values.
        results (HydrationResults): Hydration results.

    """
    def __init__(self, T1: np.array, T1_power: np.array, E: np.array,
                 E_power: np.array, hp: HydrationParameter):
        """Class Init

        Args:
            T1 (numpy.array): T1 array. Unit: second.
            T1_power (numpy.array): power in Watt unit, same length as T1.
            E (numpy.array): Enhancements.
            E_power (numpy.array): power in Watt unit, same length as E.
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
        T1interp = self.interpT1(self.E_power, self.T1_power, self.T1)
        results = self._calcODNP(self.E_power, self.E, T1interp)
        self.results = results

    def interpT1(self, power: np.array, T1power: np.array, T1p: np.array):
        """Returns the one-dimensional piecewise interpolant to a function with
        given discrete data points (T1power, T1p), evaluated at power.

        Points outside the data range will be extrapolated

        Args:
            power: The x-coordinates at which to evaluate.
            T1power: The x-coordinates of the data points, must be increasing.
                Otherwise, T1power is internally sorted.
            T1p: The y-coordinates of the data points, same length as T1power.

        Returns:
            interplatedT1 (np.array): The evaluated values, same shape as power.

        """
        T10, T100, slC = self.hp.T10, self.hp.T100, self.hp.slC

        t1_interp_method = self.hp.t1InterpMethod

        #TODO: sort T1p and T1powers

        if t1_interp_method=='2ord': # 2nd order fit, Franck and Han MIE (Eq. 22) and (Eq. 23)

            delT1w=T1p[-1]-T1p[0]  #Fixme: This requires T1p to be ascending, any better way?
            T1w=T100
            macroC=slC

            kHH = ((1./T10) - (1./(T1w))) / (macroC/1e6)
            krp=((1./T1p)-(1./(T1w + delT1w * T1power))-(kHH*(macroC/1e6))) / (slC/1e6)

            p = np.polyfit(T1power, krp, 2)
            flinear = np.polyval(p, power)

            intT1 = 1./(((slC/1e6) * flinear)+(1./(T1w + delT1w * power))+(kHH * (macroC/1e6)))

        elif t1_interp_method=='linear': # linear fit, Franck et al. PNMRS (Eq. 39)

            linearT1=1./((1./T1p)-(1./T10)+(1./T100))

            p = np.polyfit(T1power, linearT1, 1)
            flinear = np.polyval(p, power)

            intT1=flinear/(1.+(flinear/T10)-(flinear/T100))

        else:
            raise Exception('NotImplemented T1 t1InterpMethod')

        return intT1

    def _calcODNP(self, power:np.array, Ep:np.array, T1p:np.array):
        """Returns a HydrationResults object that contains all calculated ODNP values

        Following: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
        equations are labeled (#) for where they appear in the paper, sections are specified in some cases

        Args:
            power (numpy.array): Sequence of powers in Watts. Will be internally
                sorted to be ascending.
            Ep (numpy.array): Array of enhancements. Must be same length as power
            T1p (numpy.array): Array of T1. Must be same length as power.
        """
        # field and spin label concentration are defined in Hydration Parameter
        field, slC = self.hp.field, self.hp.slC

        T10 = self.hp.T10  # this is the T1 with spin label but at 0 mw power, unit is sec
        T100 = self.hp.T100  # this is the T1 without spin label and without mw power, unit is sec

        if self.hp.smaxMod == 'tethered':
            # Option 1, tether spin label
            s_max = 1  # (section 2.2) maximal saturation factor

        else:
            # Option 2, free spin probe
            s_max = 1-(2/(3+(3*(slC*1e-6*198.7)))) # from:
            # M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. &
            # J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

        omega_e = (1.76085963023e-1) * (field / 1000)
        # gamma_e in 1/ps for the tcorr unit, then correct by field in T.
        # gamma_e is from NIST. The field cancels in the following wRatio but you
        # need these individually for the spectral density functions later.

        omega_H = (2.6752218744e-4) * (field / 1000)
        # gamma_H in 1/ps for the tcorr unit, then correct by field in T.
        # gamma_H is from NIST. The field cancels in the following wRatio but you
        # need these individually for the spectral density functions later.

        wRatio = ((omega_e / (2 * np.pi)) / (omega_H / (2 * np.pi)))
        # (Eq. 4-6) ratio of omega_e and omega_H, divide by (2*pi) to get angular
        # frequency units in order to correspond to S_0/I_0, this is also ~= to the
        # ratio of the resonance frequencies for the experiment, i.e. MW freq/RF freq

        ksig_sp = (1 - Ep) / (slC * wRatio * T1p)
        # (Eq. 41) this calculates the array of k_sigma*s(p) from the enhancement array,
        # dividing by the T1p array for the "corrected" analysis

        ksigma_smax , p_12 = self.getksigsmax(ksig_sp, power)
        # fit to the right side of Eq. 42 to get (k_sigma*smax) and half of the power at s_max, called p_12 here
        
        #################
        # TODO: plot: 'Data', {power, ksig_sp}, and the 'Corrected', {power, ksigsp_fit},
        #  and 'Uncorrected', {power, ksig_sp_uncorr} that are calculated below to assess
        #  the quality of fit as well as compare "corrected" to "uncorrected" analyses to assess heating.

        ksigsp_fit = (ksigma_smax * power) / (p_12 + power)
        # (Eq. 42) calculate the "corrected" k_sigma*s(p) array using the fit parameters,
        # this can be used to plot over the data ksig_sp array to assess the quality of fit.
        # This would correspond to the corrected curves in Figure 9

        ksig_sp_uncorr = (1 - Ep) / (slC * wRatio * T10)
        # (Eq. 44) the "uncorrected" model, this can also be plotted with the corrected
        # curve to determine the severity of heating effects, as in Figure 9.
        # Notice the division by T10 instead of the T1p array
        #################

        k_sigma = ksigma_smax / s_max
        # (Eq. 43) this divides by the s_max to isolate k_sigma, the "cross" relaxivity, unit is s^-1 M^-1

        ksig_bulk = self.hp.ksig_bulk  # unit is s^-1 M^-1
        # The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
        # Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
        # 2015, 137, 12013−12023. Figure 3 caption

        k_rho = ((1/T10) - (1/T100)) / slC  # (Eq. 36) "self" relaxivity, unit is s^-1 M^-1

        ksi = k_sigma / k_rho  # (Eq. 3) this is the coupling factor, unitless

        tcorr = self.getTcorr(ksi, omega_e, omega_H)
        # (Eq. 21-23) this calls the fit to the spectral density functions. The fit
        # optimizes the value of tcorr in the calculation of ksi, the correct tcorr
        # is the one for which the calculation of ksi from the spectral density
        # functions matches the ksi found experimentally. tcorr unit is ps

        tcorr_bulk = self.hp.tcorr_bulk  # (section 2.5), "corrected" bulk tcorr, unit is ps

        dH2O = self.hp.dH2O
        # (Eq. 19-20) bulk water diffusivity, unit is d^2/s where d is distance in
        # meters. *didnt use m to avoid confusion with mass

        dSL = self.hp.dSL
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

        k_low_bulk = self.hp.k_low_bulk  # unit is s^-1 M^-1
        # The only place I can find this is Franck, JM, et. al.; "Anomalously Rapid
        # Hydration Water Diffusion Dynamics Near DNA Surfaces" J. Am. Chem. Soc.
        # 2015, 137, 12013−12023. Figure 3 caption
        
        ksi_unc , p_12_unc = self.getksiunc(Ep, power, T10, T100, wRatio, s_max)
        
        Ep_unc = 1-((ksi_unc*(1-(T10/T100))*wRatio)*((power*s_max)/(p_12_unc+power)))
        
        # this list should be in the Results object,
        # should also include flags, exceptions, etc. related to calculations
        return HydrationResults({
            'Ep_unc'   : Ep_unc,
            'T1interp' : T1p,
            'k_sigma_array' : ksig_sp / s_max,
            'k_sigma_fit' : ksigsp_fit,
            'k_sigma': k_sigma,
            'ksigma_kbulk_invratio' : 1/(k_sigma/ksig_bulk),
            'k_rho'  : k_rho,
            'k_low'  : k_low,
            'klow_klow_bulk_ratio': k_low / k_low_bulk,
            'ksi'    : ksi,
            'tcorr'  : tcorr,
            'tcorr_tcorr_bulk_ratio': tcorr / tcorr_bulk,
            'dLocal' : dLocal
        })

    @staticmethod
    def getTcorr(ksi: float, omega_e: float, omega_H: float):
        """Returns correlation time tcorr in pico second

        Args:
            ksi (float):
            omega_e (float):
            omega_H (float):

        Returns:
            tcorr (float): correlation time in pico second

        Raises:
            FitError: If no available root is found.

        """

        def f_ksi(tcorr: float, omega_e: float, omega_H: float):
            '''Returns ksi for any given tcorr

            Args:
                tcorr (float):
                omega_e (float):
                omega_H (float):

            Returns:
                ksi (float):

            '''

            # Using Franck et al. PNMRS (2013)

            zdiff = np.sqrt(1j * (omega_e - omega_H) * tcorr)
            zsum  = np.sqrt(1j * (omega_e + omega_H) * tcorr)
            zH    = np.sqrt(1j * omega_H * tcorr)

            # (Eq. 2)
            Jdiff = (1 + (zdiff / 4)) / (1 + zdiff + ((4 * (zdiff**2)) / 9) + ((zdiff**3) / 9))

            Jsum  = (1 + (zsum / 4)) / (1 + zsum + ((4 * (zsum**2)) / 9) + ((zsum**3) / 9))

            JH    = (1 + (zH / 4)) / (1 + zH + ((4 * (zH**2)) / 9) + ((zH**3) / 9))

            # (Eq. 23) calculation of ksi from the spectral density functions
            ksi_tcorr = ((6 * np.real(Jdiff)) - np.real(Jsum)) / ((6 * np.real(Jdiff)) + (3 * np.real(JH)) + np.real(Jsum))

            return ksi_tcorr

        # root finding
        # see https://docs.scipy.org/doc/scipy/reference/optimize.html
        result = optimize.root_scalar(
            lambda tcorr: f_ksi(tcorr, omega_e=omega_e, omega_H=omega_H) - ksi,
            method='brentq',
            bracket=[1, 1e5])

        if not result.converged:
            raise FitError("Could not find tcorr")
        return result.root

    @staticmethod
    def getksigsmax(ksig_sp: np.array, power: np.array):
        """Get ksig * smax and power at half ksig

        Args:
            ksig_sp (numpy.array): Array of (k_sigma * s(power)).
            power (numpy.array): Array of power.

        Returns:
            A tuple of float (ksigma_smax, p_12).

        Raises:
            FitError: If least square fitting is not succeed.

        """

        def residual(x, power: np.array, ksig_sp: np.array):
            """Residual function for ksigs_p for any given ksigma_smax and p_12

            Args:
                x (tuple): length of 2
                power (numpy.array): Array of power.
                ksig_sp (numpy.array): Array of (k_sigma * s(power)).

            Returns:
                Residuals.

            """
            ksigma_smax, p_12 = x[0], x[1]

            # Again using: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

            # Right side of Eq. 42. This function should fit to ksig_sp
            ksigsp_fit = (ksigma_smax * power) / (p_12 + power)

            return ksigsp_fit - ksig_sp

        # least-squares fitting. I like this one because it can calculate a jacobian that we can use to get an estimate of the error in k_sigma.
        # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html#scipy.optimize.least_squares
        result = optimize.least_squares(fun=residual,
                                         x0=[50, (max(power) / 2)],
                                         args=(power, ksig_sp),
                                         jac='2-point', method='lm')
        if not result.success:
            raise FitError('Could not fit ksigma ~ power')
        assert result.x[0] > 0, 'Unexpected ksigma value: %d < 0' % result.x[0]
        return result.x

    @staticmethod
    def getksiunc(Ep: np.array, power: np.array, T10: float, T100: float, wRatio: float, s_max: float):
        """Get ksi and power at half saturation

        Args:
            Ep (numpy.array): Array of enhancements.
            power (numpy.array): Array of power.
            T10 (float): T10
            T100 (float): T100
            wRatio (float): ratio of electron & proton Larmor frequencies
            s_max (float): maximal saturation factor

        Returns:
            A tuple of float (ksi, p_12).

        Raises:
            FitError: If least square fitting is not succeed.

        """

        def residual(x, Ep: np.array, power: np.array, T10: float, T100: float, wRatio: float, s_max: float):
            """Residual function for Ep for any given ksiand p_12

            Args:
                x (tuple): length of 2
                Ep (numpy.array): Array of enhancements.
                power (numpy.array): Array of power.
                T10 (float): T10
                T100 (float): T100
                wRatio (float): ratio of electron & proton Larmor frequencies
                s_max (float): maximal saturation factor

            Returns:
                Residuals.

            """
            ksi_unc, p_12_unc = x[0], x[1]

            # Again using: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

            # Right side of Eq. 42. This function should fit to ksig_sp
            Ep_fit = 1-((ksi_unc*(1-(T10/T100))*wRatio)*((power*s_max)/(p_12_unc+power)))

            return Ep - Ep_fit

        # least-squares fitting. I like this one because it can calculate a jacobian that we can use to get an estimate of the error in k_sigma.
        # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html#scipy.optimize.least_squares
        result = optimize.least_squares(fun=residual,
                                         x0=[0.5, (max(power) / 2)],
                                         args=(Ep, power, T10, T100, wRatio, s_max),
                                         jac='2-point', method='lm')
        if not result.success:
            raise FitError('Could not fit Ep')
        assert result.x[0] > 0, 'Unexpected ksi value: %d < 0' % result.x[0]
        return result.x
