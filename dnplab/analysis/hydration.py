import numpy as _np
from scipy import optimize
import warnings
from ..io.load import load
from ..io.save import save
from ..core.data import DNPData
from ..fitting.general import *
from ..math import *
from scipy.constants import *
import dnplab as _dnp
from ..config.config import DNPLAB_CONFIG


def calculate_smax(spin_C = False, smax_model = 'tethered'):
    """
    Calculate smax using radical concentration and model

    Args:
        spin_C (float): unpaired spin concentration (M)
        smax_model (str): model to determine smax (default: 'tethered')

    Returns:
        smax (float): maximal saturation factor (unitless)

    .. math::
        \mathrm{s_{max}} = 1 - (2 / (3 + (3 * (\mathrm{spin\_C} * 198.7))))

    M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. & J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.

    """
    if smax_model == "tethered":
        # Option 1, tether spin label
        smax = 1  # (section 2.2) maximal saturation factor

    elif smax_model == "free":
        # Option 2, free spin probe
        smax = 1 - (2 / (3 + (3 * (spin_C * 198.7))))  # from:
        # M.T. Türke, M. Bennati, Phys. Chem. Chem. Phys. 13 (2011) 3630. &
        # J. Hyde, J. Chien, J. Freed, J. Chem. Phys. 48 (1968) 4211.
    
    elif isinstance(smax_model, float):
        # Option 3, manual input of smax
        if not (smax_model <= 1 and smax_model > 0):
            raise ValueError(
                "if given directly, smax must be type float between 0 and 1"
            )
        smax = smax_model

    else:
        raise ValueError(
            "'smax_model' must be 'tethered', 'free', or a float between 0 and 1"
        )
    
    return smax

def t1_fit_coefficients(t1_powers, t1_arrays, degree):
    """
    Fit t1 vs. power curve and return t1 at the power of 0

    Args:
        t1_powers (numpy.array): microwave power array for t1s
        t1_arrays (numpy.array): t1s array
        degree (int): polynomials degree

    Returns:
        coeff_array: (numpy array) coefficients of t1 fit
        err_array: (numpy array) errors of t1 fit

    Methods:
        n-degree polynomials fitting
    """

    coeff_array = []
    err_array = []
    for t1_array in t1_arrays:
        coeff, cov = _np.polyfit(t1_powers, t1_array, degree, cov=True)
        err = _np.sqrt(_np.diag(cov))
        coeff_array.append(coeff)
        err_array.append(err)

    return _np.array(coeff_array), _np.array(err_array)

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


def calculate_ksigma(ksigma_sps=False, powers=False, smax=1):
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
    
    ksigma = []
    ksigma_coeff = []
    ksigma_stdd = []
        
    for ksigma_sp in ksigma_sps:
        
        popt, pcov = optimize.curve_fit(
            calculate_ksigma_array,
            powers,
            ksigma_sp,
            p0=[max(_np.real(ksigma_sp)) / 2, (max(powers) * 0.1)],
            method="lm",
        )

        assert popt[0] > 0, "Unexpected ksigma value: %d < 0" % popt[0]
        
        ksigma_coeff.append(popt)
        
        ksigma_smax = popt[0]
        p_12 = popt[1]
        
        ksigma_std = _np.sqrt(_np.diag(pcov))
        ksigma_stdd.append(ksigma_std[0] / smax)

        ksigma_fit = calculate_ksigma_array(powers, ksigma_smax, p_12)

        ksigma.append(ksigma_smax / smax)

    return _np.array(ksigma_coeff), _np.array(ksigma), _np.array(ksigma_stdd), 

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
    Jdiff = (1 + (zdiff / 4)) / (
        1 + zdiff + ((4 * (zdiff**2)) / 9) + ((zdiff**3) / 9)
    )

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
        coupling_factor (float or numpy.array): coupling factor
        omega_e (float): electron gyromagnetic ratio
        omega_H (float): proton gyromagnetic ratio

    Returns:
        tcorr (numpy.array): translational diffusion correlation time (s)


    J.M. Franck et al. / Progress in Nuclear Magnetic Resonance Spectroscopy 74 (2013) 33–56
    """

    # root finding
    # see https://docs.scipy.org/doc/scipy/reference/optimize.html

    coupling_factor = _np.array(coupling_factor)

    tcorr = []
    for index in range(len(coupling_factor)):
        if coupling_factor[index] > 0.21 and coupling_factor[index] < 0.37:
            result = optimize.root_scalar(
                lambda t_corr: calculate_xi(t_corr, omega_e=omega_e, omega_H=omega_H)
                - coupling_factor[index],
                method="brentq",
                bracket=[1, 1e5],
            )

            if not result.converged:
                raise ValueError("Could not find tcorr")
            
            root = result.root * 1e-12
            
        else:
            root = -99.            
        
        tcorr.append(root)
    return _np.array(tcorr)

def calculate_omega_ratio(field):
    """
    Calculate omega ratio using magnetic field

    Args:
        field (float): using magnetic field (T)

    Returns:
        omega_ratio (float)

    """

    if field > 3:
        warnings.warn(
            "Field (magnetic_field) should be given in T, support for mT will be removed in a future release"
        )
        field *= 1e-3

    omega_e = 1.76085963023e-1 * field
    # gamma_e in 1/ps for the tcorr unit, then correct by magnetic_field in T.
    # gamma_e is from NIST. The magnetic_field cancels in the following omega_ratio but you
    # need these individually for the spectral density functions later.

    omega_H = 2.6752218744e-4 * field
    # gamma_H in 1/ps for the tcorr unit, then correct by magnetic_field in T.
    # gamma_H is from NIST. The magnetic_field cancels in the following omega_ratio but you
    # need these individually for the spectral density functions later.

    omega_ratio = (omega_e / (2 * pi)) / (omega_H / (2 * pi))
    # (Eq. 4-6) ratio of omega_e and omega_H, divide by (2*pi) to get angular
    # frequency units in order to correspond to S_0/I_0, this is also ~= to the
    # ratio of the resonance frequencies for the experiment, i.e. MW freq/RF freq

    return omega_ratio, omega_e, omega_H

def calculate_leakage_factor(t1n0, t10n, t1n0_err=None, t10n_err=None):
    """
    Calculating leakage factor and its error

    Args:
        t1n0: (numpy array) t1 at the power of 0
        t10n: (numpy array) t10 (no radical)
        t1n0_err: (numpy array) the error of t1 at the power of 0
        t10n_err: (numpy array) errors of t10 (no radical)

    Returns:
        f: (numpy array), leakage factor
        f_err: (numpy array), leakage factor errors

    ..math::
        f = t1n0 / t10n
        f_err = leakage_factor * ((t1n0_err / t1n0) ** 2 +  (t10n_err / t10n) ** 2) ** (1 / 2)

    Reference:
    [1] 10.1016/j.pnmrs.2013.06.001
    [2] Wikipedia: 'Propagation of uncertainty'
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty

    """

    # if t1n0_err == None:
    #     t1n0_err = _np.zeros(_np.size(t1n0))
    #     print("t1n0 does not have errors, assigning 0 to all errors")

    # if t10n_err == None:
    #     t10n_err = _np.zeros(_np.size(t10n))
    #     print("t10n does not have errors, assigning 0 to all errors")

    f = [1 - x / y for x, y in zip(t1n0, t10n)]
    f_err = [
        (1 - x / y) * _np.sqrt((x_err / x) ** 2 + (y_err / y) ** 2)
        for x, x_err, y, y_err in zip(t1n0, t1n0_err, t10n, t10n_err)
    ]

    return _np.array(f), _np.array(f_err)

def calculate_krho(t1n0, t10n, radical_concentration, t1n0_err=None, t10n_err=None):
    """
    Calculating k rho and its error

    Args:
        t1n0: (numpy array) t1 at the power of 0
        t10n: (numpy array) t10 (no radical)
        radical_concentration (float): radical concentration, unit M
        t1n0_err: (numpy array) the error of t1 at the power of 0
        t10n_err: (numpy array) errors of t10 (no radical)

    Returns:
        krho: (numpy array): krho
        krho_err: (numpy array): errors of krho

    ..math::
        krho = (t1n0 ^ -1 - t10n ^ -1) / radical_conc
        krho_err = ((t1n0_err / (t1n0 ** 2)) ** 2 +  (t10n_err / (t10n ** 2)) ** 2) ** (1 / 2)


    Reference:
    [1] 10.1016/j.pnmrs.2013.06.001
    [2] Wikipedia: 'Propagation of uncertainty'
    https://en.wikipedia.org/wiki/Propagation_of_uncertainty

    """

    c = radical_concentration

    krho = [(x**-1 - y**-1) / c for x, y in zip(t1n0, t10n)]
    krho_err = [
        _np.sqrt((x_err / x**2) ** 2 + (y_err / y**2) ** 2) / c
        for x, x_err, y, y_err in zip(t1n0, t1n0_err, t10n, t10n_err)
    ]

    return _np.array(krho), _np.array(krho_err)

def calculate_coupling_factor(ksig, krho, ksig_err, krho_err):
    """
    Calculating coupling factor xi and its error

    Args:
        ksig_smax_array: (numpy array) ksigma*smax
        krho: (numpy array): krho
        ksig_smax_err_array: (numpy array) the errors of ksigma*smax
        krho_err: (numpy array): errors of krho

    Returns:
        xi_smax: (numpy array) coupling factor
        xi_smax_err: (numpy array): coupling factor errors

    ..math::

        xi_smax = ksigma_smax / krho
        xi_smax_err = xi_smax * ((ksigma_smax_err / ksigma_smax) ** 2 +  (krho_err / krho) ** 2) ** (1 / 2)

    Reference:
        [1] 10.1016/j.pnmrs.2013.06.001
        [2] Wikipedia: 'Propagation of uncertainty'
        https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    """

    xi = [x / y for x, y in zip(ksig, krho)]
    xi_err = [
        x / y * _np.sqrt((x_err / x) ** 2 + (y_err / y) ** 2)
        for x, x_err, y, y_err in zip(ksig, ksig_err, krho, krho_err)
    ]

    return _np.array(xi), _np.array(xi_err)

def hydration_analysis(path = None, smax_model = 'tethered', t10_method = 'interpolate', interpolation_degree = 1, 
                       field = 350, update_constant = {}, save_file=False, show_result=True, verbose = False):
    """
    Perform hydration analysis and return leakage factor, coupling factor, ksigma, krho, etc.

    Args:
        path (str): the path to h5 file that contains a dictionary
        save (binary): True for saving new h5 file
        print (binary): True for print hydration analysis results

    Returns:
        hydration_info(dict): DNPData objects and hydration analysis results (dict)

    ..Math::
        see examples of each hydration function
    """
    hydration_results = {}

    default_constants = {
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

    odnp_constants = {**default_constants, **update_constant}

    # load dictionary file
    hydration_info = load(path)

    # extract radical concentration
    if 'radical_concentration' in hydration_info["sample_information"]:
        radical_concentration = hydration_info["sample_information"]["radical_concentration"]
    else:
        radical_concentration = 1
        warnings.warn("'radical_concentration' was not found in 'sample_information'. The radical_concentration is set to 1 M")
    
    # calculate smax
    smax = calculate_smax(radical_concentration, smax_model)

    # calculate omega ratio
    omega_ratio, omega_e, omega_H = calculate_omega_ratio(field)
    
    # check datasets in hydration info
    check_list = ["odnp_enhancements", "odnp_ir_coefficients"] # check enhancements data and t1 vs. power data
    for key in check_list:
        if key not in hydration_info:
            raise ValueError(
                "key: %s is missing, please specify the type of radical" % key
            )
    # get enhancements data 
    enh_arrays = hydration_info["odnp_enhancements"].values.T
    enh_powers = hydration_info["odnp_enhancements"].coords['power']

    if 'odnp_fit_coefficients' in hydration_info:
        enh_max = hydration_info['odnp_fit_coefficients'].values[0]
        p_half = hydration_info['odnp_fit_coefficients'].values[1]

        if 'odnp_fit_errors' in hydration_info: # if odnp_fit_coefficients doesn't exist, it is meaningless to get errors
            enh_max_err = hydration_info['odnp_fit_errors'].values[0]
            p_half_err = hydration_info['odnp_fit_errors'].values[1]

        else:
            enh_max_err = _np.zeros(len(enh_max))
            p_half_err = _np.zeros(len(p_half))  

    # get T1s data
    t1_arrays = hydration_info["odnp_ir_coefficients"].values[0].T
    t1_powers = hydration_info["odnp_ir_coefficients"].coords['power'] 
    m0 = hydration_info["odnp_ir_coefficients"].values[1].T
    m_inf = hydration_info["odnp_ir_coefficients"].values[2].T

    if 'odnp_ir_errors' in hydration_info: # get errors for T1 
        t1_arrays_err = hydration_info["odnp_ir_errors"].values[0].T
        m0_err = hydration_info["odnp_ir_errors"].values[1].T
        m_inf_err = hydration_info["odnp_ir_errors"].values[2].T
    
    else:
        t1_arrays_err = _np.full(_np.shape(t1_arrays), 0)
        m0_err = _np.full(_np.shape(m0), 0)
        m_inf_err = _np.full(_np.shape(m_inf), 0)

    # check if enhancements data and t1 vs. power data have the same number of peak
    if _np.shape(t1_arrays)[0] != _np.shape(enh_arrays)[0]:
        raise ValueError("Peaks number is not same in 'odnp_enhancements' and 'odnp_ir_coefficients'")
    else:
        number_of_peaks = _np.shape(t1_arrays)[0]

    # get T100 
    if "T1_Water" in update_constant: # force T100 to a user-defined value
        t100 = _np.array(update_constant["T1_Water"]) 
        if len(t100) != number_of_peaks: # make the T100 length is equal as the number of peaks
            t100 = _np.append(t100, _np.full(number_of_peaks - len(t100), 2.5))
            t100_message = "aren't given in updated constant completely, not given t100s are forces to 2.5 s."
        else:
            t100_message = "are given in updated constant completely."

    elif "ir_coefficients" in hydration_info: # T100 from dictionary
        t100 = hydration_info["ir_coefficients"].values[0]
        t100_message = 'are calculated from dictionary.' # a message to indicate where T100 comes
    else: # when T100 is not given
        t100 = _np.full(number_of_peaks, 2.5) # create default an array for T100 
        t100_message = 'are not given and forced to 2.5 s.'

    # get T100 errors
    if "ir_errors" not in hydration_info:
        t100_err = _np.zeros(number_of_peaks) # force errors to 0
    else:
        t100_err = hydration_info["ir_errors"].values[0] # create default an array for T100 errors 
    
    # interpolate T1 vs. power curves
    t1_coeff, t1_coeff_err = t1_fit_coefficients(t1_powers, t1_arrays, interpolation_degree)

    # interpolate T10 (T1 at 0 dBm) from previous curves
    interp_t10 = t1_coeff[:, 1][:]
    interp_t10_err = t1_coeff_err[:, 1][:]

    # get t10s from dictionary if t10s exist
    if _np.min(t1_powers) == 0:
        raw_t10 = t1_arrays[:, _np.where(t1_powers == 0)]
        raw_t10_err = t1_arrays_err[:, _np.where(t1_powers == 0)]
        raw_t10_flag = True
    else:
        raw_t10_flag = False

    # select which t10s will be used for hydration analysis
    if t10_method.lower() == 'interpolate':
        t10 = interp_t10
        t10_err = interp_t10_err
    elif t10_method.lower() == 'raw' and raw_t10_flag:
        t10 = raw_t10
        t10_err = raw_t10_err
    elif t10_method.lower() == 'raw' and raw_t10_flag != False:
        raise ValueError("T10 method is 'raw', but T10s cannot be found in T1 data.")
    else:
        raise ValueError("T10 methods are 'interpolate' or 'raw'.")        

    # find ksigma*s(p)
    # (Eq. 41) this calculates the array of ksigma*s(p) from the enhancement array,
    # dividing by the T1 array for the "corrected" analysis
    ksigma_array = _np.array([(1 - enh_arrays[x]) / (radical_concentration * omega_ratio * _np.poly1d(t1_coeff[x])(enh_powers)) for x in range(number_of_peaks)])

    # calculate ksigma
    ksigma_sp_coefficients, ksigma, ksigma_err = calculate_ksigma(
        ksigma_array, enh_powers, smax
    )

    # calcualte leakage factor
    f, f_err = calculate_leakage_factor(t10, t100, t10_err, t100_err)

    # calculate krho
    krho, krho_err = calculate_krho(
        t10, t100, radical_concentration, t10_err, t100_err
    )

    # calculate coupling factors
    xi, xi_err = calculate_coupling_factor(ksigma, krho, ksigma_err, krho_err)

    # calculate translational correction time (tcorr) in pico second
    tcorr= calculate_tcorr(xi, omega_e, omega_H)

    # calculate D_local
    Dlocal = (odnp_constants["tcorr_bulk"] / tcorr) * (
        odnp_constants["D_H2O"] + odnp_constants["D_SL"]
    )

    if show_result == True:
        print("----------------------------------------------------------")
        print("Hydration Analysis Results:")
        print("----------------------------------------------------------")
        print("Number of Peaks: %d" %number_of_peaks)
        print("T1 Intepolation Degree: %d" %interpolation_degree)
        print("Radical Concentration: %0.03f M" %radical_concentration)
        print("Smax Model: %s" %smax_model)
        print("Smax: %0.03f" %smax)
        print("T100 %s" %t100_message) 

        print("----------------------------------------------------------")

        for index, integrals in enumerate(
            hydration_info["odnp_enhancements"].coords["integrals"]
        ):
            peak_number = index + 1
            print("Peak NO.%i:" % peak_number)
            print('Emax: %0.03f +/- %0.03f s' % (enh_max[index], enh_max_err[index]))
            print('p_half: %0.03f +/- %0.03f s' % (p_half[index], p_half_err[index]))
            if raw_t10_flag:
                print("Raw T10: %0.03f +/- %0.03f s" % (raw_t10[index], raw_t10_err[index]))
            else:
                print('Raw T10 is not found')
            print("Interpolated T10: %0.03f +/- %0.03f s" % (interp_t10[index], interp_t10_err[index]))
            print("T100 %0.03f +/- %0.03f s" %(t100[index], t100_err[index]))
            print("Leakage Factor: %0.03f +/- %0.03f" % (f[index], f_err[index]))
            print(
                "ksigma: %0.03f +/- %0.03f"
                % (ksigma[index], ksigma_err[index])
            )
            print("krho: %0.03f +/- %0.03f" % (krho[index], krho_err[index]))
            print(
                "Coupling Factor: %0.03f +/- %0.03f"
                % (xi[index], xi_err[index])
            )
            if tcorr[index] < 0:
                print('No Translation Correction Time (tcorr) because coupling factor is out of range.')
            else:
                print(
                "Estimated Translation Correlation Time: %0.03f ps"
                % (tcorr[index]*1e12)
            )
                
            if Dlocal[index] < 0:
                print('No Dlocal.')
            else:
                print(
                "Dlocal: %0.03fe-9 m^2s-1"
                % (Dlocal[index]*1e9)
            )

            print("----------------------------------------------------------")
    
    # assign values in hydration results    
    hydration_results["t1_coefficients"] = t1_coeff
    hydration_results["t1_coefficients_errors"] = t1_coeff_err

    hydration_results["t10"] = t10
    hydration_results["t10_errors"] = t10_err

    hydration_results["ksigma_sp_coefficients"] = ksigma_sp_coefficients
    hydration_results["ksigma"] = ksigma
    hydration_results["ksigma_errors"] = ksigma_err

    hydration_results["leakage_factors"] = f
    hydration_results["leakage_factor_errors"] = f_err

    hydration_results["krho"] = krho
    hydration_results["krho_errors"] = krho_err

    hydration_results["coupling_factor"] = xi
    hydration_results["coupling_factor_errors"] = xi_err

    hydration_results["translation_correlation_time"] = tcorr

    hydration_results["Dlocal"] = Dlocal

    hydration_info["hydration_results"] = hydration_results

    if verbose == True:
        color_cycler_list = [
            DNPLAB_CONFIG.get("COLORS", color_key)
            for color_key in DNPLAB_CONFIG["COLORS"].keys()
        ]

        # plot enhancements curves
        _dnp.plt.figure("odnp_enhancements")
        hydration_info["odnp_enhancements"].attrs['experiment_type'] = 'enhancements_PW'
        for index in range(number_of_peaks):
            _dnp.fancy_plot(hydration_info["odnp_enhancements"]['integrals', index], marker = 'o', ls = '', color = color_cycler_list[index], label = 'Peak No.%d: $E_{max}$ = %0.03f +/- %0.03f' %(index+1, enh_max[index], enh_max_err[index])) 
            _dnp.plt.plot(_np.linspace(0, _np.max(enh_powers)), _dnp.relaxation.buildup_function(_np.linspace(0, _np.max(enh_powers)), enh_max[index], p_half[index]), color = color_cycler_list[index])
        _dnp.plt.legend()
        _dnp.plt.tight_layout()

        # plot t1 curves
        hydration_info['t1_integrals'].attrs['experiment_type'] = 'inversion_recovery'
        for index in range(number_of_peaks):
            _dnp.plt.figure("inversion_recovery, peak No.%d" %(index+1))
            for i, power in enumerate(t1_powers):
                _dnp.fancy_plot(hydration_info['t1_integrals']['integrals', index, 'power', power], marker = 'o', ls = '', color = color_cycler_list[i], label = 'MW: %0.01f dBm, $T_1$ = %0.03f +/- %0.03f s' %(_dnp.w2dBm(power), t1_arrays[index,i], t1_arrays_err[index,i])) 
                _dnp.plt.plot(_np.linspace(0, _np.max(hydration_info['t1_integrals'].coords['t1'])), _dnp.relaxation.t1(_np.linspace(0, _np.max(hydration_info['t1_integrals'].coords['t1'])), t1_arrays[index, i], m0[index, i], m_inf[index, i]), color = color_cycler_list[i])
            _dnp.plt.legend()
            _dnp.plt.tight_layout()

        # plot t1 vs.power
        _dnp.plt.figure("t1_interpolation")
        for index in range(number_of_peaks):
            _dnp.plt.plot(t1_powers, t1_arrays[index], color = color_cycler_list[index], marker = 'o', ls = '', fillstyle = 'none')
            _dnp.plt.plot(0, interp_t10[index], color = color_cycler_list[index], marker = 'o', ls = '',  label = "Peak No.%d: Interpolated $T_{10}$ =  %0.03f +/- %0.03f s" %(index+1, interp_t10[index], interp_t10_err[index]))
            if raw_t10_flag:
                _dnp.plt.plot(0, raw_t10[index], color = color_cycler_list[index], marker = 's', ls = '',  label = "Peak No.%d: raw $T_{10}$ =  %0.03f +/- %0.03f s" %(index+1, raw_t10[index], raw_t10_err[index]))
            _dnp.plt.plot(_np.linspace(0, _np.max(enh_powers)), _np.poly1d(t1_coeff[index])(_np.linspace(0, _np.max(enh_powers))), color = color_cycler_list[index])
        _dnp.plt.title('$T_1$ Interpolation, degree = %d' %interpolation_degree)
        _dnp.plt.xlabel('Microwave Power (W)')
        _dnp.plt.ylabel('$T_1$ (s)')
        _dnp.plt.grid(ls = ':')
        _dnp.plt.legend()
        _dnp.plt.tight_layout()
          

        # plot ksigma curve
        _dnp.plt.figure("ksigma")
        # no fancy plot for ksigma_array
        for index in range(number_of_peaks):
            _dnp.plt.plot(enh_powers, ksigma_array[index], color = color_cycler_list[index], marker = 'o', ls = '', fillstyle = 'none', label = 'Peak No.%d: $k_{\sigma}$ = %0.03f +/- %0.03f' %(index +1, ksigma[index], ksigma_err[index]))
            _dnp.plt.plot(_np.linspace(0, _np.max(enh_powers)), calculate_ksigma_array(_np.linspace(0, _np.max(enh_powers)), ksigma_sp_coefficients[index][0], ksigma_sp_coefficients[index][1]))
        _dnp.plt.title('$k_{\sigma}$ vs. Power')
        _dnp.plt.xlabel('Microwave Power (W)')
        _dnp.plt.ylabel('$k_{\sigma}$S(p) (a.u.)')
        _dnp.plt.grid(ls = ':')
        _dnp.plt.legend()
        _dnp.plt.tight_layout()
        _dnp.plt.show()  

    if save_file == True:
        
        save(hydration_info, path, overwrite=True)

    return hydration_info



## ready to be depreciarted.
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
    http://dx.doi.org/10.1016/j.pnmrs.2013.06.001

    J.M. Franck, S. Han; Methods in Enzymology, Chapter 5, Volume 615, (2019) 131-175
    https://doi.org/10.1016/bs.mie.2018.09.024
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

    omega_ratio = (omega_e / (2 * pi)) / (omega_H / (2 * pi))
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

