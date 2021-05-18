""" dnpHydration module

This module calculates hydration related quantities using processed ODNP data.

"""

import numpy as np
from scipy import interpolate
from scipy import optimize


def calculate_smax(spin_C=False):
    """Returns maximal saturation factor according to: M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. & J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

    .. math::
        s_{max} = 1 - (2 / (3 + (3 * (spin_C * 198.7))))

    Args:
        spin_C: unpaired spin concentration in units of uM

    Returns:
        smax (float): maximal saturation factor
    """

    return 1 - (2 / (3 + (3 * (spin_C * 198.7))))


def interpolate_T1(
    E_powers=False,
    T1_powers=False,
    T1_array=False,
    interpolate_method="linear",
    delta_T1_water=False,
    T1_water=False,
    macro_C=False,
    spin_C=1,
    T10=2.0,
    T100=2.5,
):
    """Returns interpolated T1 data using Eq. 39 of http://dx.doi.org/10.1016/j.pnmrs.2013.06.001 for "linear" or Eq. 22 of https://doi.org/10.1016/bs.mie.2018.09.024 for "second_order"

    Args:
        E_powers: The x-coordinates at which to evaluate
        T1_powers: The x-coordinates of the data points, must be increasing
            Otherwise, T1_power is internally sorted
        T1_array: The y-coordinates of the data points, same length as T1_power
        interpolate_method: "second_order" or "linear"
        delta_T1_water: change in T1 of water at max microwave power
        T1_water: T1 of pure water
        macro_C: concentration of macro molecule
        spin_C: unpaired electron spin concentration in uM
        T10: T1 measured with unpaired electrons
        T100: T1 measured without unpaired electrons

    Returns:
        interpolated_T1 (np.array): The evaluated values, same shape as E_powers.
    """

    spin_C = spin_C / 1e6

    # 2nd order fit, Franck and Han MIE (Eq. 22) and (Eq. 23)
    if interpolate_method == "second_order":
        if not delta_T1_water:
            delta_T1_water = T1_array[-1] - T1_array[0]
        if not T1_water:
            T1_water = T100
        if not macro_C:
            macro_C = spin_C

        kHH = (1.0 / T10 - 1.0 / T1_water) / macro_C
        krp = (
            (1.0 / T1_array)
            - (1.0 / (T1_water + delta_T1_water * T1_powers))
            - (kHH * (macro_C))
        ) / (spin_C)

        p = np.polyfit(T1_powers, krp, 2)
        T1_fit_2order = np.polyval(p, E_powers)

        interpolated_T1 = 1.0 / (
            ((spin_C) * T1_fit_2order)
            + (1.0 / (T1_water + delta_T1_water * E_powers))
            + (kHH * (macro_C))
        )

    # linear fit, Franck et al. PNMRS (Eq. 39)
    elif interpolate_method == "linear":

        linear_t1 = 1.0 / ((1.0 / T1_array) - (1.0 / T10) + (1.0 / T100))

        p = np.polyfit(T1_powers, linear_t1, 1)
        T1_fit_linear = np.polyval(p, E_powers)

        interpolated_T1 = T1_fit_linear / (
            1.0 + (T1_fit_linear / T10) - (T1_fit_linear / T100)
        )

    else:
        raise Exception("invalid interp_method")

    return interpolated_T1


def calculate_ksigma_array(powers=False, ksigma_smax=95.4, p_12=False):
    """Function to calcualte ksig array for any given ksigma and p_12

    Again using: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

    Args:
        powers (numpy.array): Array of powers
        ksigma_smax (float): product of ksigma and smax
        p_12 (float): power at half max for ksigma fit

    Returns:
        ksig_fit (numpy.array): calculated ksigma array

    """

    # Right side of Eq. 42. This function should fit to ksig_sp
    ksig_fit = (ksigma_smax * powers) / (p_12 + powers)

    return ksig_fit


def calculate_ksigma(ksigma_sp=False, powers=False, smax=1):
    """Get ksigma and E_power at half max of ksig

    Args:
        ksig (numpy.array): Array of ksigma.
        powers (numpy.array): Array of E_power.

    Returns:
        ksigma (float): calculated ksigma
        ksigma_stdd (float): standard deviation in calculated ksigma
        p_12 (float): power at half max for ksigma fit

    Asserts:
        ksigma (popt[0]) is greater than zero

    """

    # curve fitting
    # see https://docs.scipy.org/doc/scipy/reference/optimize.html
    popt, pcov = optimize.curve_fit(
        calculate_ksigma_array,
        powers,
        ksigma_sp,
        p0=[95.4 / 2, (max(powers) * 0.1)],
        method="lm",
    )

    assert popt[0] > 0, "Unexpected ksigma value: %d < 0" % popt[0]

    ksigma_smax = popt[0]
    p_12 = popt[1]
    ksigma_std = np.sqrt(np.diag(pcov))
    ksigma_stdd = ksigma_std[0] / smax

    ksigma_fit = calculate_ksigma_array(powers, ksigma_smax, p_12)

    ksigma = ksigma_smax / smax

    return ksigma, ksigma_stdd, ksigma_fit


def calculate_xi(tcorr=54, omega_e=0.0614, omega_H=9.3231e-05):
    """Returns coupling_factor for any given tcorr

    Args:
        tcorr (float): translational diffusion correlation time
        omega_e (float): electron gyromagnetic ratio
        omega_H (float): proton gyromagnetic ratio

    Returns:
        xi (float): coupling factor

    """

    # Using Franck et al. PNMRS (2013)

    zdiff = np.sqrt(1j * (omega_e - omega_H) * tcorr)
    zsum = np.sqrt(1j * (omega_e + omega_H) * tcorr)
    zH = np.sqrt(1j * omega_H * tcorr)

    # (Eq. 2)
    Jdiff = (1 + (zdiff / 4)) / (
        1 + zdiff + ((4 * (zdiff ** 2)) / 9) + ((zdiff ** 3) / 9)
    )

    Jsum = (1 + (zsum / 4)) / (1 + zsum + ((4 * (zsum ** 2)) / 9) + ((zsum ** 3) / 9))

    JH = (1 + (zH / 4)) / (1 + zH + ((4 * (zH ** 2)) / 9) + ((zH ** 3) / 9))

    # (Eq. 23) calculation of coupling_factor from the spectral density functions
    xi = ((6 * np.real(Jdiff)) - np.real(Jsum)) / (
        (6 * np.real(Jdiff)) + (3 * np.real(JH)) + np.real(Jsum)
    )

    return xi


def calculate_tcorr(coupling_factor=0.27, omega_e=0.0614, omega_H=9.3231e-05):
    """Returns correlation time tcorr in pico second

    Args:
        coupling_factor (float): coupling factor
        omega_e (float): electron gyromagnetic ratio
        omega_H (float): proton gyromagnetic ratio

    Returns:
        result.root (float): tcorr, translational diffusion correlation time in pico second

    Raises:
        FitError: If no available root is found.

    """

    # root finding
    # see https://docs.scipy.org/doc/scipy/reference/optimize.html
    result = optimize.root_scalar(
        lambda tcorr: calculate_xi(tcorr, omega_e=omega_e, omega_H=omega_H)
        - coupling_factor,
        method="brentq",
        bracket=[1, 1e5],
    )

    if not result.converged:
        raise FitError("Could not find tcorr")
    return result.root


def calculate_uncorrected_Ep(
    uncorrected_xi=0.33,
    p_12_unc=0,
    E_powers=False,
    T10=2.0,
    T100=2.5,
    omega_ratio=658.5792,
    smax=1,
):
    """Function for Ep for any given xi and p_12

    Again using: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

    Args:
        uncorrected_xi (float): uncorrected coupling factor
        p_12_unc (float): power at half max for uncorrected_xi fit
        E_array (numpy.array): Array of enhancements.
        E_powers (numpy.array): Array of E_power.
        T10 (float): T10
        T100 (float): T100
        omega_ratio (float): ratio of electron & proton Larmor frequencies
        smax (float): maximal saturation factor

    Returns:
        Ep_fit (array): uncorrected Enhancement curve

    """

    # Right side of Eq. 42. This function should fit to ksig_sp
    Ep_fit = 1 - (
        (uncorrected_xi * (1 - (T10 / T100)) * omega_ratio)
        * ((E_powers * smax) / (p_12_unc + E_powers))
    )

    return Ep_fit


def _residual_Ep(
    x,
    E_array: np.array,
    E_powers: np.array,
    T10: float,
    T100: float,
    omega_ratio: float,
    smax: float,
):
    """Function for Ep for any given xi and p_12

    Again using: J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

    Args:
        x (list): (uncorrected coupling factor, power at half max for uncorrected_xi fit)
        E_array (numpy.array): Array of enhancements.
        E_powers (numpy.array): Array of E_power.
        T10 (float): T10
        T100 (float): T100
        omega_ratio (float): ratio of electron & proton Larmor frequencies
        smax (float): maximal saturation factor

    Returns:
        Ep_fit (array): uncorrected Enhancement curve

    """

    return E_array - calculate_uncorrected_Ep(
        uncorrected_xi=x[0],
        p_12_unc=x[1],
        E_powers=E_powers,
        T10=T10,
        T100=T100,
        omega_ratio=omega_ratio,
        smax=smax,
    )


def calculate_uncorrected_xi(
    E_array=False,
    E_powers=False,
    T10=2.0,
    T100=2.5,
    omega_ratio=658.5792,
    smax=1,
):
    """Get coupling_factor and E_power at half saturation

    Args:
        E_array (numpy.array): Array of enhancements
        E_powers (numpy.array): Array of powers
        T10 (float): T10
        T100 (float): T100
        omega_ratio (float): ratio of electron & proton Larmor frequencies
        smax (float): maximal saturation factor

    Returns:
        uncorrected_xi (float): uncorrected coupling factor
        p_12_unc (float): power at half max for uncorrected_xi fit

    Raises:
        FitError: If least square fitting is not succeed.

    """

    # least-squares fitting.
    # see https://docs.scipy.org/doc/scipy/reference/optimize.html
    results = optimize.least_squares(
        fun=_residual_Ep,
        x0=[0.27, (max(E_powers) * 0.1)],
        args=(E_array, E_powers, T10, T100, omega_ratio, smax),
        jac="2-point",
        method="lm",
    )
    if not results.success:
        raise FitError("Could not fit Ep")
    assert results.x[0] > 0, "Unexpected coupling_factor value: %d < 0" % results.x[0]

    uncorrected_xi = results.x[0]
    p_12_unc = results.x[1]

    return uncorrected_xi, p_12_unc


def odnp(
    inputs={},
    constants={
        "ksigma_bulk": 95.4,
        "krho_bulk": 353.4,
        "klow_bulk": 366,
        "tcorr_bulk": 54,
        "D_H2O": 2.3e-9,
        "D_SL": 4.1e-10,
        "delta_T1_water": False,
        "T1_water": False,
        "macro_C": False,
    },
):
    """Function for performing ODNP calculations

    Args:
        inputs (dict): keys and values below
        uncorrected_Ep (numpy.array)    : fit of enhancement array,
        uncorrected_xi (float)          : coupling factor for fit of enhancement array,
        interpolated_T1 (numpy.array)   : T1 values interpolated on enhancement powers,
        ksigma_array (numpy.array)      : ksigma array
        ksigma_fit (numpy.array)        : fit of ksigma_array,
        ksigma (float)                  : ksigma,
        ksigma_stdd (float)             : standard deviation in ksigma,
        ksigma_bulk_ratio (float)       : ratio ksigma / ksigma_bulk,
        krho (float)                    : krho,
        krho_bulk_ratio (float)         : ratio krho / krho_bulk,
        klow (float)                    : klow,
        klow_bulk_ratio (float)         : ratio klow / klow_bulk,
        coupling_factor (float)         : coupling_factor,
        tcorr (float)                   : tcorr,
        tcorr_bulk_ratio (float)        : ratio tcorr / tcorr_bulk,
        Dlocal (float)                  : Dlocal
        constants (dict): keys and values below
        ksigma_bulk (float): bulk ksigma value
        krho_bulk (float): bulk krho value
        klow_bulk (float): bulk klow value
        tcorr_bulk (float): bulk tcorr value
        D_H2O (float): diffusivity of bulk water
        D_SL (float): diffusivity of spin probe in bulk water
        delta_T1_water (float): change in T1 of water of E_powers range
        T1_water (float): T1 of water protons
        macro_C (float): concentration of macromolecule

    Returns:
        results (dict): keys and values below
        uncorrected_Ep (numpy.array): uncorrected_Ep,
        uncorrected_xi (float): xi_unc,
        interpolated_T1 (numpy.array): T1p,
        ksigma_array (numpy.array): ksigma_array,
        ksigma_fit (numpy.array): ksigma_fit,
        ksigma (float): ksigma,
        ksigma_stdd (float): ksigma_stdd,
        ksigma_bulk_ratio (float): ksigma / constants["ksigma_bulk"],
        krho (float): krho,
        krho_bulk_ratio (float): krho / constants["krho_bulk"],
        klow (float): klow,
        klow_bulk_ratio (float): klow / constants["klow_bulk"],
        coupling_factor (float): coupling_factor,
        tcorr (float): tcorr,
        tcorr_bulk_ratio (float): tcorr / constants["tcorr_bulk"],
        Dlocal (float): Dlocal
    """

    if not inputs:
        raise ValueError("Please supply a valid inputs dictionary")

    inputs["spin_C"] /= 1e6

    if inputs["smax_model"] == "tethered":
        # Option 1, tether spin label
        s_max = 1  # (section 2.2) maximal saturation factor

    elif inputs["smax_model"] == "free":
        # Option 2, free spin probe
        s_max = calculate_smax(inputs["spin_C"])  # from:
        # M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. &
        # J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

    if isinstance(inputs["smax_model"], int) or isinstance(inputs["smax_model"], float):
        # Option 3, manual input of smax
        if not (inputs["smax_model"] <= 1 and inputs["smax_model"] > 0):
            raise ValueError("smax must be a number between 0 and 1")
        s_max = inputs["smax_model"]

    omega_e = (1.76085963023e-1) * (inputs["field"] / 1000)
    # gamma_e in 1/ps for the tcorr unit, then correct by field in T.
    # gamma_e is from NIST. The field cancels in the following omega_ratio but you
    # need these individually for the spectral density functions later.

    omega_H = (2.6752218744e-4) * (inputs["field"] / 1000)
    # gamma_H in 1/ps for the tcorr unit, then correct by field in T.
    # gamma_H is from NIST. The field cancels in the following omega_ratio but you
    # need these individually for the spectral density functions later.

    omega_ratio = (omega_e / (2 * np.pi)) / (omega_H / (2 * np.pi))
    # (Eq. 4-6) ratio of omega_e and omega_H, divide by (2*pi) to get angular
    # frequency units in order to correspond to S_0/I_0, this is also ~= to the
    # ratio of the resonance frequencies for the experiment, i.e. MW freq/RF freq

    if "T1_powers" in inputs.keys():
        T1p = interpolate_T1(
            E_powers=inputs["E_powers"],
            T1_powers=inputs["T1_powers"],
            T1_array=inputs["T1_array"],
            interpolate_method=inputs["interpolate_method"],
            delta_T1_water=constants["delta_T1_water"],
            T1_water=constants["T1_water"],
            macro_C=constants["macro_C"],
            spin_C=inputs["spin_C"],
            T10=inputs["T10"],
            T100=inputs["T100"],
        )
    else:
        if len(inputs["T1_array"]) == len(inputs["E_array"]):
            T1p = inputs["T1_array"]
        else:
            raise ValueError(
                "'T1_array' must be equal in length to 'E_array'. Otherwise give 'T1_powers' equal in length to 'T1_array' to interpolate."
            )

    ksigma_array = (1 - inputs["E_array"]) / (inputs["spin_C"] * omega_ratio * T1p)
    # (Eq. 41) this calculates the array of ksigma*s(p) from the enhancement array,
    # dividing by the T1 array for the "corrected" analysis

    ksigma, ksigma_stdd, ksigma_fit = calculate_ksigma(
        ksigma_array, inputs["E_powers"], s_max
    )
    # fit to the right side of Eq. 42 to get (ksigma*smax) and half of the E_power at s_max, called p_12 here

    krho = ((1 / inputs["T10"]) - (1 / inputs["T100"])) / inputs[
        "spin_C"
    ]  # "self" relaxivity, unit is s^-1 M^-1

    coupling_factor = ksigma / krho  # coupling factor, unitless

    tcorr = calculate_tcorr(coupling_factor, omega_e, omega_H)
    # (Eq. 21-23) this calls the fit to the spectral density functions. The fit
    # optimizes the value of tcorr in the calculation of coupling_factor, the correct tcorr
    # is the one for which the calculation of coupling_factor from the spectral density
    # functions matches the coupling_factor found experimentally. tcorr unit is ps

    Dlocal = (constants["tcorr_bulk"] / tcorr) * (
        constants["D_H2O"] + constants["D_SL"]
    )
    # (Eq. 19-20) local diffusivity, i.e. diffusivity of the water near the spin label

    ############################################################################
    # This is defined in its most compact form in:
    # Frank, JM and Han, SI;  Chapter Five - Overhauser Dynamic Nuclear Polarization
    # for the Study of Hydration Dynamics, Explained. Methods in Enzymology, Volume 615, 2019
    #
    # But also explained well in:
    # Franck, JM, et. al.; "Anomalously Rapid Hydration Water Diffusion Dynamics
    # Near DNA Surfaces" J. Am. Chem. Soc. 2015, 137, 12013−12023.

    klow = ((5 * krho) - (7 * ksigma)) / 3
    # section 6, (Eq. 13). this describes the relatively slowly diffusing water
    # near the spin label, sometimes called "bound" water
    ############################################################################

    xi_unc, p_12_unc = calculate_uncorrected_xi(
        inputs["E_array"],
        inputs["E_powers"],
        inputs["T10"],
        inputs["T100"],
        omega_ratio,
        s_max,
    )

    uncorrected_Ep = calculate_uncorrected_Ep(
        xi_unc,
        p_12_unc,
        inputs["E_powers"],
        inputs["T10"],
        inputs["T100"],
        omega_ratio,
        s_max,
    )

    return {
        "uncorrected_Ep": uncorrected_Ep,
        "uncorrected_xi": xi_unc,
        "interpolated_T1": T1p,
        "ksigma_array": ksigma_array,
        "ksigma_fit": ksigma_fit,
        "ksigma": ksigma,
        "ksigma_stdd": ksigma_stdd,
        "ksigma_bulk_ratio": ksigma / constants["ksigma_bulk"],
        "krho": krho,
        "krho_bulk_ratio": krho / constants["krho_bulk"],
        "klow": klow,
        "klow_bulk_ratio": klow / constants["klow_bulk"],
        "coupling_factor": coupling_factor,
        "tcorr": tcorr,
        "tcorr_bulk_ratio": tcorr / constants["tcorr_bulk"],
        "Dlocal": Dlocal,
    }


def hydration(ws):
    """Function for calculating hydration quantities

    Args:
        ws (dnpdata_collection): see function 'odnp' below

    Returns:
        results (dict): see function 'odnp' below
    """

    if "hydration_inputs" in ws.keys():

        odnp_constants = {
            "ksigma_bulk": 95.4,
            "krho_bulk": 353.4,
            "klow_bulk": 366,
            "tcorr_bulk": 54,
            "D_H2O": 2.3e-9,
            "D_SL": 4.1e-10,
            "delta_T1_water": False,
            "T1_water": False,
            "macro_C": False,
        }

        if "hydration_constants" in ws.keys():
            for ky in odnp_constants.keys():
                if ky in ws["hydration_constants"].keys():
                    odnp_constants[ky] = ws["hydration_constants"][ky]

        odnp_inputs = ws["hydration_inputs"]

        results = odnp(odnp_inputs, odnp_constants)

        ws["hydration_results"] = results

        return results

    else:
        raise TypeError("the hydration_inputs dictionary is missing!")
