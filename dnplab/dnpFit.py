from . import dnpdata as dnpdata, dnpdata_collection
import numpy as _np
from scipy.optimize import curve_fit


def return_data(all_data):

    is_workspace = False
    if isinstance(all_data, dnpdata):
        data = all_data.copy()
    elif isinstance(all_data, dict):
        raise ValueError("Type dict is not supported")
    elif isinstance(all_data, dnpdata_collection):
        is_workspace = True
        if all_data.processing_buffer in all_data.keys():
            data = all_data[all_data.processing_buffer]
        else:
            raise ValueError("No data in processing buffer")
    else:
        raise ValueError("Data type not supported")

    return data, is_workspace


def t1_function(t_axis, T1, M_0, M_inf):
    return M_0 - M_inf * _np.exp(-1.0 * t_axis / T1)


def t2_function_stretch(x_axis, M_0, T2, p):
    return M_0 * _np.exp(-2.0 * (x_axis / T2) ** p)


def t2_function_nostretch(x_axis, M_0, T2):
    return M_0 * _np.exp(-2.0 * (x_axis / T2) ** 1.0)


def exp_fit_func_1(x_axis, C1, C2, tau):
    return C1 + C2 * _np.exp(-1.0 * x_axis / tau)


def exp_fit_func_2(x_axis, C1, C2, tau1, C3, tau2):
    return C1 + C2 * _np.exp(-1.0 * x_axis / tau1) + C3 * _np.exp(-1.0 * x_axis / tau2)


def exponential_fit(
    all_data, type="mono", stretched=False, dim="f2", indirect_dim=None
):
    """Fits various forms of exponential functions

    .. math::

        f(t) = M_0 - M_{\infty} e^{-t/T_1}
        f(t) = M_0 e^{(-2(t/T_2)^p}
        f(t) = M_0 e^{(-2(t/T_2)}
        f(t) = C1 + C2 e^{-x/tau}
        f(t) = C1 + C2 e^{-x/tau1} + C3 e^{-x/tau2}

    Args:
        all_data (dnpdata, dict): data container after processing inversion recovery data, after integration with dnpNMR.integrate
        type (str) : "T1" for inversion recovery fit, "T2" for stretched exponential, "mono", or "bi"
        stretch (boolean) : if False "p" is set to 1, if True "p" is a fit parameter
        dim (str) : dimension to fit down

    Returns:
        dnpdata_collection or dnpdata: Processed data in container, updated with fit data
        attributes: "T1" value and "T1_stdd" standard deviation for type="T1", "T2" value and "T2_stdd" standard deviation for type="T2", or "tau" and "tau_stdd" for type="mono", "tau1", "tau1_stdd", "tau2", and "tau2_stdd" for type="bi"

    """

    data, isDict = return_data(all_data)

    if not indirect_dim:
        if len(data.dims) == 2:
            ind_dim = list(set(data.dims) - set([dim]))[0]
        elif len(data.dims) == 1:
            ind_dim = data.dims[0]
        else:
            raise ValueError(
                "you must specify the indirect dimension, use argument indirect_dim= "
            )
    else:
        ind_dim = indirect_dim

    x_axis = data.coords[ind_dim]
    new_axis = _np.r_[_np.min(x_axis) : _np.max(x_axis) : 100j]
    inputData = _np.real(data.values)

    if type == "T1":

        x0 = [1.0, inputData[-1], inputData[-1]]
        out, cov = curve_fit(t1_function, x_axis, inputData, x0, method="lm")
        stdd = _np.sqrt(_np.diag(cov))
        fit = t1_function(new_axis, out[0], out[1], out[2])

        fitData = dnpdata(fit, [new_axis], [ind_dim])
        fitData.attrs["T1"] = out[0]
        fitData.attrs["T1_stdd"] = stdd[0]
        fitData.attrs["M_0"] = out[1]
        fitData.attrs["M_inf"] = out[2]

    elif type == "T2":

        if stretched:
            x0 = [inputData[0], 1.0, 1.0]
            out, cov = curve_fit(
                t2_function_stretch, x_axis, inputData, x0, method="lm"
            )
            stdd = _np.sqrt(_np.diag(cov))
            fit = t2_function_stretch(new_axis, out[0], out[1], out[2])
        else:
            x0 = [inputData[0], 1.0]
            out, cov = curve_fit(
                t2_function_nostretch, x_axis, inputData, x0, method="lm"
            )
            stdd = _np.sqrt(_np.diag(cov))
            fit = t2_function_nostretch(new_axis, out[0], out[1])

        fitData = dnpdata(fit, [new_axis], [ind_dim])
        fitData.attrs["T2"] = out[1]
        fitData.attrs["T2_stdd"] = stdd[1]
        fitData.attrs["M_0"] = out[0]
        if stretched:
            fitData.attrs["p"] = out[2]

    elif type == "mono":

        x0 = [inputData[-1], 1.0, 100]
        out, cov = curve_fit(exp_fit_func_1, x_axis, inputData, x0, method="lm")
        stdd = _np.sqrt(_np.diag(cov))
        fit = exp_fit_func_1(new_axis, out[0], out[1], out[2])

        fitData = dnpdata(fit, [new_axis], [ind_dim])
        fitData.attrs["tau"] = out[2]
        fitData.attrs["tau_stdd"] = stdd[2]
        fitData.attrs["C1"] = out[0]
        fitData.attrs["C2"] = out[1]

    elif type == "bi":

        x0 = [inputData[-1], 1.0, 100, 1.0, 100]
        out, cov = curve_fit(exp_fit_func_2, x_axis, inputData, x0, method="lm")
        stdd = _np.sqrt(_np.diag(cov))
        fit = exp_fit_func_2(new_axis, out[0], out[1], out[2], out[3], out[4])

        fitData = dnpdata(fit, [new_axis], [ind_dim])
        fitData.attrs["tau1"] = out[2]
        fitData.attrs["tau1_stdd"] = stdd[2]
        fitData.attrs["tau2"] = out[4]
        fitData.attrs["tau2_stdd"] = stdd[4]
        fitData.attrs["C1"] = out[0]
        fitData.attrs["C2"] = out[1]
        fitData.attrs["C3"] = out[3]

    else:
        raise TypeError("Invalid fit type")

    if isDict:
        all_data["fit"] = fitData
        return all_data
    else:
        return fitData


def enhancement_function(power_array, E_max, power_half):
    return E_max * power_array / (power_half + power_array)


def enhancement_fit(dataDict):
    """Fits enhancement curves to return Emax and power and one half maximum saturation

    .. math::

        f(p) = E_{max} p / (p_{1/2} + p)

    Args:
        workspace

    Returns:
        all_data (dnpdata, dict): Processed data in container, updated with fit data
        attributes: Emax value and Emax standard deviation

                    p_one_half value and p_one_half standard deviation

    Example::

        ### INSERT importing and processing ###
        dnplab.dnpNMR.integrate(workspace, {})

        workspace.new_dim('power', power_list)

        dnplab.dnpFit.enhancementFit(workspace)

        Emax_value = workspace['fit'].attrs['E_max']
        Emax_standard_deviation = workspace['fit'].attrs['E_max_stdd']
        p_one_half_value = workspace['fit'].attrs['p_half']
        p_one_half_standard_deviation = workspace['fit'].attrs['p_half_stdd']
        Emax_fit = workspace['fit'].values
        Emax_fit_xaxis = workspace['fit'].coords

    """

    data, isDict = return_data(all_data)

    power_axes = data.coords["power"]

    inputData = _np.real(data.values)

    x0 = [inputData[-1], 0.1]

    out, cov = curve_fit(enhancement_function, power_axes, inputData, x0, method="lm")
    stdd = _np.sqrt(_np.diag(cov))

    fit = enhancement_function(power_axes, out[0], out[1])

    fitData = dnpdata(fit, [power_axes], ["power"])
    fitData.attrs["E_max"] = out[0]
    fitData.attrs["E_max_stdd"] = stdd[0]
    fitData.attrs["power_half"] = out[1]
    fitData.attrs["power_half_stdd"] = stdd[1]

    if isDict:
        dataDict["fit"] = fitData
        return dataDict
    else:
        return fitData


def interpolate_T1(
    E_powers=False,
    T1_powers=False,
    T1_array=False,
    interp_method="linear",
    spin_C=100,
    T10=2.0,
    T100=2.5,
):
    """Returns the one-dimensional piecewise interpolant to a function with
    given discrete data points (T1_powers, T1), evaluated at E_powers.

    Points outside the data range will be extrapolated

    Args:
        E_powers: The x-coordinates at which to evaluate.
        T1_powers: The x-coordinates of the data points, must be increasing.
            Otherwise, T1_power is internally sorted.
        T1_array: The y-coordinates of the data points, same length as T1_power.
        interp_method: "second_order" or "linear".
        spin_C: unpaired electron spin concentration in uM.
        T10: T1 measured with unpaired electrons.
        T100: T1 measured without unpaired electrons.

    Returns:
        interp_T1 (np.array): The evaluated values, same shape as E_powers.

    """

    spin_C = spin_C / 1e6

    # 2nd order fit, Franck and Han MIE (Eq. 22) and (Eq. 23)
    if interp_method == "second_order":

        delta_T1_water = T1_array[-1] - T1_array[0]
        T1_water = T100
        macro_C = spin_C

        kHH = (1.0 / T10 - 1.0 / T1_water) / macro_C
        krp = (
            (1.0 / T1_array)
            - (1.0 / (T1_water + delta_T1_water * T1_powers))
            - (kHH * (macro_C))
        ) / (spin_C)

        p = _np.polyfit(T1_powers, krp, 2)
        T1_fit_2order = _np.polyval(p, E_powers)

        interp_T1 = 1.0 / (
            ((spin_C) * T1_fit_2order)
            + (1.0 / (T1_water + delta_T1_water * E_powers))
            + (kHH * (macro_C))
        )

    # linear fit, Franck et al. PNMRS (Eq. 39)
    elif interp_method == "linear":

        linear_t1 = 1.0 / ((1.0 / T1_array) - (1.0 / T10) + (1.0 / T100))

        p = _np.polyfit(T1_powers, linear_t1, 1)
        T1_fit_linear = _np.polyval(p, E_powers)

        interp_T1 = T1_fit_linear / (
            1.0 + (T1_fit_linear / T10) - (T1_fit_linear / T100)
        )

    else:
        raise Exception("invalid interp_method")

    return interp_T1
