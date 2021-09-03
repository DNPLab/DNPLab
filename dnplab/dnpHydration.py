import numpy as np
from scipy import interpolate
from scipy import optimize


def calculate_smax(spin_C=False):
    """Returns maximal saturation factor according to: M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. & J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

    .. math::
        \mathrm{s_{max}} = 1 - (2 / (3 + (3 * (\mathrm{spin\_C} * 198.7))))

    Args:
        spin_C (float): unpaired spin concentration in units of uM

    Returns:
        smax (float): maximal saturation factor

    """

    return 1 - (2 / (3 + (3 * (spin_C * 1e-6 * 198.7))))


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
        E_powers (numpy.array): The microwave powers at which to evaluate
        T1_powers (numpy.array): The microwave powers of the T1s to interpolate
        T1_array (numpy.array): The original T1s
        interpolate_method (str): "second_order" or "linear"
        spin_C (float): unpaired electron spin concentration in uM
        T10 (float): T1 measured with unpaired electrons
        T100 (float): T1 measured without unpaired electrons
        delta_T1_water (optional) (float): change in T1 of water at max microwave power
        T1_water (optional) (float): T1 of pure water
        macro_C (optional) (float): concentration of macromolecule in uM

    Returns:
        interpolated_T1 (numpy.array): Array of T1 values same shape as E_powers and E_array

    """

    # 2nd order fit, Franck and Han MIE (Eq. 22) and (Eq. 23)
    if interpolate_method == "second_order":
        spin_C = spin_C / 1e6
        if macro_C:
            macro_C = macro_C / 1e6
        else:
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
        raise Exception("invalid interpolate_method")

    return interpolated_T1


def calculate_ksigma_array(powers=False, ksigma_smax=95.4, p_12=False):
    """Function to calcualte ksig array for any given ksigma and p_12

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

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

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

    Args:
        ksig (numpy.array): Array of ksigmas
        powers (numpy.array): Array of E_powers

    Returns:
        ksigma (float): calculated ksigma
        ksigma_stdd (float): standard deviation in ksigma
        p_12 (float): power at half max for ksigma fit

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

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

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
    """Returns translational correlation time (tcorr) in pico second

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

    Args:
        coupling_factor (float): coupling factor
        omega_e (float): electron gyromagnetic ratio
        omega_H (float): proton gyromagnetic ratio

    Returns:
        t_corr (float): tcorr, translational diffusion correlation time in pico second

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

    t_corr = result.root
    return t_corr


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

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

    Args:
        uncorrected_xi (float): uncorrected coupling factor
        p_12_unc (float): power at half max for uncorrected_xi fit
        E_array (numpy.array): Array of enhancements
        E_powers (numpy.array): Array of E_powers
        T10 (float): T1(0), proton T1 with microwave power=0
        T100 (float): T10(0), proton T1 with spin_C=0 and microwave power=0
        omega_ratio (float): ratio of electron & proton gyromagnetic ratios
        smax (float): maximal saturation factor

    Returns:
        Ep_fit (numpy.array): uncorrected Enhancement curve

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
    """Function for residuals between E(p) for any given xi and p_12 and the experimental E_array

    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

    Args:
        x (list): [uncorrected coupling factor, power at half max for uncorrected_xi fit]
        E_array (numpy.array): Array of enhancements
        E_powers (numpy.array): Array of E_power
        T10 (float): T1(0), proton T1 with microwave power=0
        T100 (float): T10(0), proton T1 with spin_C=0 and microwave power=0
        omega_ratio (float): ratio of electron & proton gyromagnetic ratios
        smax (float): maximal saturation factor

    Returns:
        Ep_fit (numpy.array): uncorrected enhancement curve

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

    J.M. Franck et al.; Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56

    Args:
        E_array (numpy.array): Array of enhancements
        E_powers (numpy.array): Array of powers
        T10 (float): T1(0), proton T1 with microwave power=0
        T100 (float): T10(0), proton T1 with spin_C=0 and microwave power=0
        omega_ratio (float): ratio of electron & proton gyromagnetic ratios
        smax (float): maximal saturation factor

    Returns:
        uncorrected_xi (float): uncorrected coupling factor
        p_12_unc (float): power at half max for uncorrected_xi fit

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


def odnp(inputs={}, constants={}):
    """Function for performing ODNP calculations

    J.M. Franck et al.; Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    http://dx.doi.org/10.1016/j.pnmrs.2013.06.001

    J.M. Franck, S. Han; Methods in Enzymology, Chapter 5, Volume 615, (2019) 131-175
    https://doi.org/10.1016/bs.mie.2018.09.024

    Args:
        inputs (dict)                   : keys and values described in example above
        constants (optional) (dict)     : keys and values described in example above

    Returns:
        hydration_results (dict)        : keys and values described in table above

    """

    if not inputs:
        raise ValueError("Please supply a valid inputs dictionary")

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
    # these constants have been compiled from the various ODNP literature

    if constants:
        for ky in odnp_constants.keys():
            if ky in constants.keys():
                odnp_constants[ky] = constants[ky]

    if inputs["smax_model"] == "tethered":
        # Option 1, tether spin label
        s_max = 1  # (section 2.2) maximal saturation factor

    elif inputs["smax_model"] == "free":
        # Option 2, free spin probe
        s_max = calculate_smax(inputs["spin_C"])  # from:
        # M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. &
        # J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

    if isinstance(inputs["smax_model"], (int, float)):
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
            delta_T1_water=odnp_constants["delta_T1_water"],
            T1_water=odnp_constants["T1_water"],
            macro_C=odnp_constants["macro_C"],
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

    ksigma_array = (1 - inputs["E_array"]) / (
        inputs["spin_C"] * 1e-6 * omega_ratio * T1p
    )
    # (Eq. 41) this calculates the array of ksigma*s(p) from the enhancement array,
    # dividing by the T1 array for the "corrected" analysis

    ksigma, ksigma_stdd, ksigma_fit = calculate_ksigma(
        ksigma_array, inputs["E_powers"], s_max
    )
    # fit to the right side of Eq. 42 to get (ksigma*smax) and half of the E_power at s_max, called p_12 here

    krho = ((1 / inputs["T10"]) - (1 / inputs["T100"])) / (
        inputs["spin_C"] * 1e-6
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
        inputs["E_array"],
        inputs["E_powers"],
        inputs["T10"],
        inputs["T100"],
        omega_ratio,
        s_max,
    )
    # (Eqs. 7 and 44) this calculates the coupling factor using the "uncorrected" analysis

    uncorrected_Ep = calculate_uncorrected_Ep(
        xi_unc,
        p_12_unc,
        inputs["E_powers"],
        inputs["T10"],
        inputs["T100"],
        omega_ratio,
        s_max,
    )
    # (Eqs. 7 and 44) this calculates the "uncorrected" enhnacement array using xi_unc

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


def hydration(workspace):
    """Function for calculating hydration quantities

    J.M. Franck et al.; Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    http://dx.doi.org/10.1016/j.pnmrs.2013.06.001

    J.M. Franck, S. Han; Methods in Enzymology, Chapter 5, Volume 615, (2019) 131-175
    https://doi.org/10.1016/bs.mie.2018.09.024

    Args:
        workspace (dnpdata_collection): workspace or dictionary with 'hydration_inputs', see above

    Returns:
        results (dict)                : 'hydration_results' dictionary, see above

    Raises:
        TypeError: If 'hydration_inputs' dictionary is missing

    """

    if "hydration_inputs" in workspace.keys():

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

        if "hydration_constants" in workspace.keys():
            for ky in odnp_constants.keys():
                if ky in workspace["hydration_constants"].keys():
                    odnp_constants[ky] = workspace["hydration_constants"][ky]

        odnp_inputs = workspace["hydration_inputs"]

        results = odnp(odnp_inputs, odnp_constants)

        workspace["hydration_results"] = results

        return results

    else:
        raise TypeError("the 'hydration_inputs' dictionary is missing!")
