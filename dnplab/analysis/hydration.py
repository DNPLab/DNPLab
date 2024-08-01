import numpy as _np
from scipy import optimize
import warnings
from ..constants import constants as _const


def calculate_smax(spin_C=False):
    r"""Returns maximal saturation factor.

    Args:
        spin_C (float): unpaired spin concentration (M)

    Returns:
        smax (float): maximal saturation factor (unitless)

    .. math::
        \mathrm{s_{max}} = 1 - (2 / (3 + (3 * (\mathrm{spin\_C} * 198.7))))

    M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. & J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.
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
    """Returns interpolated T1 data.

    Args:
        E_powers (numpy.array): The microwave powers at which to evaluate
        T1_powers (numpy.array): The microwave powers of the T1s to interpolate
        T1_array (numpy.array): The original T1s (s)
        interpolate_method (str): "second_order" or "linear"
        spin_C (float): unpaired electron spin concentration (M)
        T10 (float): T1 measured with unpaired electrons (s)
        T100 (float): T1 measured without unpaired electrons (s)
        delta_T1_water (optional) (float): change in T1 of water at max microwave power (s)
        T1_water (optional) (float): T1 of pure water (s)
        macro_C (optional) (float): concentration of macromolecule (M)

    Returns:
        interpolated_T1 (numpy.array): Array of T1 values same shape as E_powers and E_array

    T1 data is interpolated using Eq. 39 of http://dx.doi.org/10.1016/j.pnmrs.2013.06.001 for "linear" or Eq. 22 of https://doi.org/10.1016/bs.mie.2018.09.024 for "second_order"
    """

    # 2nd order fit, Franck and Han MIE (Eq. 22) and (Eq. 23)
    if interpolate_method == "second_order":
        if not macro_C:
            macro_C = spin_C

        if not delta_T1_water:
            delta_T1_water = T1_array[-1] - T1_array[0]
        if not T1_water:
            T1_water = T100

        kHH = (1.0 / T10 - 1.0 / T1_water) / macro_C
        krp = (
            (1.0 / T1_array)
            - (1.0 / (T1_water + delta_T1_water * T1_powers))
            - (kHH * (macro_C))
        ) / (spin_C)

        p = _np.polyfit(T1_powers, krp, 2)
        T1_fit_2order = _np.polyval(p, E_powers)

        interpolated_T1 = 1.0 / (
            ((spin_C) * T1_fit_2order)
            + (1.0 / (T1_water + delta_T1_water * E_powers))
            + (kHH * (macro_C))
        )

    # linear fit, Franck et al. PNMRS (Eq. 39)
    elif interpolate_method == "linear":
        linear_t1 = 1.0 / ((1.0 / T1_array) - (1.0 / T10) + (1.0 / T100))

        p = _np.polyfit(T1_powers, linear_t1, 1)
        T1_fit_linear = _np.polyval(p, E_powers)

        interpolated_T1 = T1_fit_linear / (
            1.0 + (T1_fit_linear / T10) - (T1_fit_linear / T100)
        )

    else:
        raise Exception("invalid interpolate_method")

    return interpolated_T1


def calculate_ksigma_array(powers=False, ksigma_smax=95.4, p_12=False):
    """Function to calcualte ksig array for any given ksigma and p_12

    Args:
        powers (numpy.array): Array of powers
        ksigma_smax (float): product of ksigma and smax (s^-1 * M^-1)
        p_12 (float): power at half max for ksigma fit

    Returns:
        ksig_fit (numpy.array): calculated ksigma array

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    """

    # Right side of Eq. 42. This function should fit to ksig_sp
    ksig_fit = (ksigma_smax * powers) / (p_12 + powers)

    return ksig_fit


def calculate_ksigma(ksigma_sp=False, powers=False, smax=1):
    """Get ksigma and E_power at half max of ksig

    Args:
        ksig (numpy.array): Array of ksigmas
        powers (numpy.array): Array of E_powers

    Returns:
        ksigma (float): calculated ksigma
        ksigma_stdd (float): standard deviation in ksigma
        p_12 (float): power at half max for ksigma fit

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
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
    ksigma_std = _np.sqrt(_np.diag(pcov))
    ksigma_stdd = ksigma_std[0] / smax

    ksigma_fit = calculate_ksigma_array(powers, ksigma_smax, p_12)

    ksigma = ksigma_smax / smax

    return ksigma, ksigma_stdd, ksigma_fit


def calculate_xi(tcorr=54e-12, omega_e=0.0614, omega_H=9.3231e-05):
    """Returns coupling_factor for any given tcorr

    Args:
        tcorr (float): translational diffusion correlation time (s)
        omega_e (float): electron gyromagnetic ratio
        omega_H (float): proton gyromagnetic ratio

    Returns:
        xi (float): coupling factor

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    """

    # Using Franck et al. PNMRS (2013)
    if tcorr < 0.1:
        tcorr *= 1e12

    zdiff = _np.sqrt(1j * (omega_e - omega_H) * tcorr)
    zsum = _np.sqrt(1j * (omega_e + omega_H) * tcorr)
    zH = _np.sqrt(1j * omega_H * tcorr)

    # (Eq. 2)
    Jdiff = (1 + (zdiff / 4)) / (1 + zdiff + ((4 * (zdiff**2)) / 9) + ((zdiff**3) / 9))

    Jsum = (1 + (zsum / 4)) / (1 + zsum + ((4 * (zsum**2)) / 9) + ((zsum**3) / 9))

    JH = (1 + (zH / 4)) / (1 + zH + ((4 * (zH**2)) / 9) + ((zH**3) / 9))

    # (Eq. 23) calculation of coupling_factor from the spectral density functions
    xi = ((6 * _np.real(Jdiff)) - _np.real(Jsum)) / (
        (6 * _np.real(Jdiff)) + (3 * _np.real(JH)) + _np.real(Jsum)
    )

    return xi


def calculate_tcorr(coupling_factor=0.27, omega_e=0.0614, omega_H=9.3231e-05):
    """Returns translational correlation time (tcorr) in pico second

    Args:
        coupling_factor (float): coupling factor
        omega_e (float): electron gyromagnetic ratio
        omega_H (float): proton gyromagnetic ratio

    Returns:
        tcorr (float): translational diffusion correlation time (s)

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    """

    # root finding
    # see https://docs.scipy.org/doc/scipy/reference/optimize.html
    result = optimize.root_scalar(
        lambda t_corr: calculate_xi(t_corr, omega_e=omega_e, omega_H=omega_H)
        - coupling_factor,
        method="brentq",
        bracket=[1, 1e5],
    )

    if not result.converged:
        raise ValueError("Could not find tcorr")

    tcorr = result.root
    return tcorr * 1e-12


def calculate_uncorrected_Ep(
    uncorrected_xi=0.33,
    p_12_unc=0,
    E_powers=False,
    T10=2.0,
    T100=2.5,
    omega_ratio=658.5792,
    smax=1,
):
    """Function for E(p) for any given xi and p_12

    Args:
        uncorrected_xi (float): uncorrected coupling factor
        p_12_unc (float): power at half max for uncorrected_xi fit
        E_array (numpy.array): Array of enhancements
        E_powers (numpy.array): Array of E_powers
        T10 (float): T1(0), proton T1 with microwave power=0 (s)
        T100 (float): T10(0), proton T1 with spin_C=0 and microwave power=0 (s)
        omega_ratio (float): ratio of electron & proton gyromagnetic ratios
        smax (float): maximal saturation factor

    Returns:
        Ep_fit (numpy.array): uncorrected enhancement curve

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    """

    # Right side of Eq. 42. This function should fit to ksig_sp
    Ep_fit = 1 - (
        (uncorrected_xi * (1 - (T10 / T100)) * omega_ratio)
        * ((E_powers * smax) / (p_12_unc + E_powers))
    )

    return Ep_fit


def _residual_Ep(
    x,
    E_array: _np.array,
    E_powers: _np.array,
    T10: float,
    T100: float,
    omega_ratio: float,
    smax: float,
):
    """Function for residuals between E(p) for any given xi and p_12 and the experimental E_array

    Args:
        x (list): [uncorrected coupling factor, power at half max for uncorrected_xi fit]
        E_array (numpy.array): Array of enhancements
        E_powers (numpy.array): Array of E_power
        T10 (float): T1(0), proton T1 with microwave power=0 (s)
        T100 (float): T10(0), proton T1 with spin_C=0 and microwave power=0 (s)
        omega_ratio (float): ratio of electron & proton gyromagnetic ratios
        smax (float): maximal saturation factor

    Returns:
        Ep_fit (numpy.array): uncorrected enhancement curve

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
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
        T10 (float): T1(0), proton T1 with microwave power=0 (s)
        T100 (float): T10(0), proton T1 with spin_C=0 and microwave power=0 (s)
        omega_ratio (float): ratio of electron & proton gyromagnetic ratios
        smax (float): maximal saturation factor

    Returns:
        uncorrected_xi (float): uncorrected coupling factor
        p_12_unc (float): power at half max for uncorrected_xi fit

    J.M. Franck et al.; Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
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
        raise ValueError("Could not fit Ep")
    assert results.x[0] > 0, "Unexpected coupling_factor value: %d < 0" % results.x[0]

    uncorrected_xi = results.x[0]
    p_12_unc = results.x[1]

    return uncorrected_xi, p_12_unc


def hydration(data={}, constants={}):
    """Function for performing ODNP calculations

    Args:
        data (dict)                   : keys and values are described in the example
        constants (dict)              : (optional) keys and values are described in the example

    Returns:
        (dict)                        : keys and values are described in the example

    J.M. Franck et al.; Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    https://www.sciencedirect.com/science/article/abs/pii/S0079656513000629

    J.M. Franck, S. Han; Methods in Enzymology, Chapter 5, Volume 615, (2019) 131-175
    https://www.sciencedirect.com/science/article/abs/pii/S0076687918303872
    """

    if not data:
        raise ValueError("Please supply a valid data dictionary, see example")

    if "hydration_inputs" in data.keys():
        warnings.warn(
            "The workspace concept is depreciated, see the example in the docs for the new syntax"
        )
        _data = data["hydration_inputs"]
        if "hydration_constants" in data.keys():
            constants = data["hydration_constants"]
        data = _data

    if "tcorr_bulk" in constants.keys() and constants["tcorr_bulk"] > 0.1:
        warnings.warn(
            "tcorr_bulk should be given in seconds, support for picoseconds will be removed in a future release"
        )
        constants["tcorr_bulk"] *= 1e-12

    if "macro_C" in constants.keys() and constants["macro_C"] > 0.1:
        warnings.warn(
            "macro_C should be given in molar, support for micromolar will be removed in a future release"
        )
        constants["macro_C"] *= 1e-6

    if data["spin_C"] > 0.1:
        warnings.warn(
            "spin_C should be given in molar, support for micromolar will be removed in a future release"
        )
        data["spin_C"] *= 1e-6

    if "field" in data.keys():
        warnings.warn(
            "keyword 'field' is depreciated, please use 'magnetic_field' from now on"
        )
        if "magnetic_field" in data.keys():
            warnings.warn(
                "you supplied both 'field' and 'magnetic_field', only 'magnetic_field' will be used"
            )
        else:
            data["magnetic_field"] = data["field"]
        data.pop("field")

    if data["magnetic_field"] > 3:
        warnings.warn(
            "magnetic_field should be given in T, support for mT will be removed in a future release"
        )
        data["magnetic_field"] *= 1e-3

    standard_constants = {
        "ksigma_bulk": 95.4,
        "krho_bulk": 353.4,
        "klow_bulk": 366,
        "tcorr_bulk": 54e-12,
        "D_H2O": 2.3e-9,
        "D_SL": 4.1e-10,
        "delta_T1_water": False,
        "T1_water": False,
        "macro_C": False,
    }
    # these constants have been compiled from the various ODNP literature

    odnp_constants = {**standard_constants, **constants}

    if data["smax_model"] == "tethered":
        # Option 1, tether spin label
        s_max = 1  # (section 2.2) maximal saturation factor

    elif data["smax_model"] == "free":
        # Option 2, free spin probe
        s_max = calculate_smax(data["spin_C"])  # from:
        # M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. &
        # J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

    elif isinstance(data["smax_model"], float):
        # Option 3, manual input of smax
        if not (data["smax_model"] <= 1 and data["smax_model"] > 0):
            raise ValueError(
                "if given directly, smax must be type float between 0 and 1"
            )
        s_max = data["smax_model"]
    else:
        raise ValueError(
            "'smax_model' must be 'tethered', 'free', or a float between 0 and 1"
        )

    omega_e = 1.76085963023e-1 * data["magnetic_field"]
    # gamma_e in 1/ps for the tcorr unit, then correct by magnetic_field in T.
    # gamma_e is from NIST. The magnetic_field cancels in the following omega_ratio but you
    # need these individually for the spectral density functions later.

    omega_H = 2.6752218744e-4 * data["magnetic_field"]
    # gamma_H in 1/ps for the tcorr unit, then correct by magnetic_field in T.
    # gamma_H is from NIST. The magnetic_field cancels in the following omega_ratio but you
    # need these individually for the spectral density functions later.

    omega_ratio = (omega_e / (2 * _const.pi)) / (omega_H / (2 * _const.pi))
    # (Eq. 4-6) ratio of omega_e and omega_H, divide by (2*pi) to get angular
    # frequency units in order to correspond to S_0/I_0, this is also ~= to the
    # ratio of the resonance frequencies for the experiment, i.e. MW freq/RF freq

    if "T1_powers" in data.keys():
        T1p = interpolate_T1(
            E_powers=data["E_powers"],
            T1_powers=data["T1_powers"],
            T1_array=data["T1_array"],
            interpolate_method=data["interpolate_method"],
            delta_T1_water=odnp_constants["delta_T1_water"],
            T1_water=odnp_constants["T1_water"],
            macro_C=odnp_constants["macro_C"],
            spin_C=data["spin_C"],
            T10=data["T10"],
            T100=data["T100"],
        )
    else:
        if len(data["T1_array"]) == len(data["E_array"]):
            T1p = data["T1_array"]
        else:
            raise ValueError(
                "'T1_array' must be equal in length to 'E_array'. Otherwise give 'T1_powers' equal in length to 'T1_array' in order to interpolate."
            )

    ksigma_array = (1 - data["E_array"]) / (data["spin_C"] * omega_ratio * T1p)
    # (Eq. 41) this calculates the array of ksigma*s(p) from the enhancement array,
    # dividing by the T1 array for the "corrected" analysis

    ksigma, ksigma_stdd, ksigma_fit = calculate_ksigma(
        ksigma_array, data["E_powers"], s_max
    )
    # fit to the right side of Eq. 42 to get (ksigma*smax) and half of the E_power at s_max, called p_12 here

    krho = ((1 / data["T10"]) - (1 / data["T100"])) / (
        data["spin_C"]
    )  # (Eq. 36) "self" relaxivity, unit is s^-1 M^-1

    coupling_factor = ksigma / krho  # coupling factor, unitless

    tcorr = calculate_tcorr(coupling_factor, omega_e, omega_H)
    # (Eq. 21-23) this calls the fit to the spectral density functions. The fit
    # optimizes the value of tcorr in the calculation of coupling_factor, the correct tcorr
    # is the one for which the calculation of coupling_factor from the spectral density
    # functions matches the coupling_factor found experimentally. tcorr unit is ps

    Dlocal = (odnp_constants["tcorr_bulk"] / tcorr) * (
        odnp_constants["D_H2O"] + odnp_constants["D_SL"]
    )
    # (Eq. 19-20) local diffusivity, i.e. diffusivity of the water near the spin label

    klow = ((5 * krho) - (7 * ksigma)) / 3
    # section 6, (Eq. 13). this describes the relatively slowly diffusing water
    # near the spin label, sometimes called "bound" water.
    # This is defined in its most compact form in:
    # Frank, JM and Han, SI;  Chapter Five - Overhauser Dynamic Nuclear Polarization
    # for the Study of Hydration Dynamics, Explained. Methods in Enzymology, Volume 615, 2019
    # But also explained well in:
    # Franck, JM, et. al.; "Anomalously Rapid Hydration Water Diffusion Dynamics
    # Near DNA Surfaces" J. Am. Chem. Soc. 2015, 137, 12013−12023.

    xi_unc, p_12_unc = calculate_uncorrected_xi(
        data["E_array"],
        data["E_powers"],
        data["T10"],
        data["T100"],
        omega_ratio,
        s_max,
    )
    # (Eqs. 7 and 44) this calculates the coupling factor using the "uncorrected" analysis

    uncorrected_Ep = calculate_uncorrected_Ep(
        xi_unc,
        p_12_unc,
        data["E_powers"],
        data["T10"],
        data["T100"],
        omega_ratio,
        s_max,
    )
    # (Eqs. 7 and 44) this calculates the "uncorrected" enhancement array using xi_unc

    return {
        "uncorrected_Ep": uncorrected_Ep,
        "uncorrected_xi": xi_unc,
        "interpolated_T1": T1p,
        "ksigma_array": ksigma_array,
        "ksigma_fit": ksigma_fit,
        "ksigma": ksigma,
        "ksigma_stdd": ksigma_stdd,
        "ksigma_bulk_ratio": ksigma / odnp_constants["ksigma_bulk"],
        "krho": krho,
        "krho_bulk_ratio": krho / odnp_constants["krho_bulk"],
        "klow": klow,
        "klow_bulk_ratio": klow / odnp_constants["klow_bulk"],
        "coupling_factor": coupling_factor,
        "tcorr": tcorr,
        "tcorr_bulk_ratio": tcorr / odnp_constants["tcorr_bulk"],
        "Dlocal": Dlocal,
    }
