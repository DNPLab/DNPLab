from . import dnpdata as _dnpdata, dnpdata_collection
import numpy as _np

from scipy.optimize import curve_fit


def t1Function(t, T1, M_0, M_inf):

    return M_0 - M_inf * _np.exp(-1.0 * t / T1)


def t1Fit(dataDict):
    """Fits inversion recovery data to extract T1 value in seconds

    .. math::

        f(t) = M_0 - M_{\infty} e^{-t/T_1}

    Args:
        workspace after processing inversion recovery data, after integration with dnpNMR.integrate

    Returns:
        all_data (dnpdata, dict): Processed data in container, updated with fit data
        attributes: T1 value and T1 standard deviation

    Example:

    .. code-block:: python

        ### INSERT importing and processing ###
        dnplab.dnpNMR.integrate(workspace, {})

        dnplab.dnpFit.t1Fit(workspace)

        T1_value = workspace['fit'].attrs['t1']
        T1_standard_deviation = workspace['fit'].attrs['t1_stdd']
        T1_fit = workspace['fit'].values
        T1_fit_xaxis = workspace['fit'].coords

    """

    isDict = False
    if isinstance(dataDict, (dict, dnpdata_collection)):
        data = dataDict["proc"].copy()
        isDict = True
    elif isinstance(dataDict, _dnpdata):
        data = dataDict.copy()
    else:
        print("Incompatible data type:")
        print(type(dataDict))
        return

    t1_axes = data.coords["t1"]

    inputData = _np.real(data.values)

    x0 = [1.0, inputData[-1], inputData[-1]]

    out, cov = curve_fit(t1Function, t1_axes, inputData, x0, method="lm")
    stdd = _np.sqrt(_np.diag(cov))

    new_axes = _np.r_[_np.min(t1_axes) : _np.max(t1_axes) : 100j]
    fit = t1Function(new_axes, out[0], out[1], out[2])

    fitData = _dnpdata(fit, [new_axes], ["t1"])
    fitData.attrs["t1"] = out[0]
    fitData.attrs["t1_stdd"] = stdd[0]
    fitData.attrs["M_0"] = out[1]
    fitData.attrs["M_inf"] = out[2]

    if isDict:
        dataDict["fit"] = fitData
        return dataDict
    else:
        return fitData


def enhancementFunction(powerArray, E_max, power_half):

    return E_max * powerArray / (power_half + powerArray)


def enhancementFit(dataDict):
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

    isDict = False
    if isinstance(dataDict, (dict, dnpdata_collection)):
        data = dataDict["proc"].copy()
        isDict = True
    elif isinstance(dataDict, _dnpdata):
        data = dataDict.copy()
    else:
        print("Incompatible data type:")
        print(type(dataDict))
        return

    power_axes = data.coords["power"]

    inputData = _np.real(data.values)

    x0 = [inputData[-1], 0.1]

    out, cov = curve_fit(enhancementFunction, power_axes, inputData, x0, method="lm")
    stdd = _np.sqrt(_np.diag(cov))

    fit = enhancementFunction(power_axes, out[0], out[1])

    fitData = _dnpdata(fit, [power_axes], ["power"])
    fitData.attrs["E_max"] = out[0]
    fitData.attrs["E_max_stdd"] = stdd[0]
    fitData.attrs["power_half"] = out[1]
    fitData.attrs["power_half_stdd"] = stdd[1]

    if isDict:
        dataDict["fit"] = fitData
        return dataDict
    else:
        return fitData
