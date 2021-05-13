from . import dnpdata as dnpdata, dnpdata_collection
from . import dnpMath
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


def exponential_fit(
    all_data,
    type="mono",
    stretched=False,
    dim="t1",
    indirect_dim=None,
    ws_key="integrals",
):
    """Fits various forms of exponential functions

    .. math::

        f(t) &= M_0 - M_{\infty} e^{-t/T_1} \\
        
             &= M_0 e^{(-2(t/T_2)^p} \\
        
             &= C1 + C2 e^{-x/tau} \\
        
             &= C1 + C2 e^{-x/tau1} + C3 e^{-x/tau2}


    Args:
        all_data (dnpdata, dict): data container, after integration with dnpTools.integrate
        
    +---------------+------+-------------+-----------------------------------------------------------------------------------+
    | parameter     | type | default     | description                                                                       |
    +===============+======+=============+===================================================================================+
    | dim           | str  | "f2"        | dimension to fit down                                                             |
    +---------------+------+-------------+-----------------------------------------------------------------------------------+
    | type          | str  | "mono"      | "T1" for inversion recovery fit, "T2" for stretched exponential, "mono", or "bi"  |
    +---------------+------+-------------+-----------------------------------------------------------------------------------+
    | stretch       | bool | False       | if False "p" is set to 1, if True "p" is a fit parameter                          |
    +---------------+------+-------------+-----------------------------------------------------------------------------------+
    | indirect_dim  | str  | None        | indirect dimension                                                                |
    +---------------+------+-------------+-----------------------------------------------------------------------------------+
    | ws_key        | str  | "integrals" | if False "p" is set to 1, if True "p" is a fit parameter                          |
    +---------------+------+-------------+-----------------------------------------------------------------------------------+


    Returns:
        dnpdata_collection or dnpdata: Processed data in container, updated with fit data
        attributes: "T1" value and "T1_stdd" standard deviation for type="T1", "T2" value and "T2_stdd" standard deviation for type="T2", "tau" and "tau_stdd" for type="mono", or "tau1", "tau1_stdd", "tau2", and "tau2_stdd" for type="bi"
    """

    data, isDict = return_data(all_data)

    if ws_key in all_data.keys():
        x_axis = all_data[ws_key].coords[dim]
        new_axis = _np.r_[_np.min(x_axis) : _np.max(x_axis) : 100j]
        input_data = _np.real(all_data[ws_key].values)
        ind_dim = dim
    else:
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
        input_data = _np.real(data.values)

    if type == "T1":

        x0 = [1.0, input_data[-1], input_data[-1]]
        out, cov = curve_fit(dnpMath.t1_function, x_axis, input_data, x0, method="lm")
        stdd = _np.sqrt(_np.diag(cov))
        fit = dnpMath.t1_function(new_axis, out[0], out[1], out[2])

        fit_data = dnpdata(fit, [new_axis], [ind_dim])
        fit_data.attrs["T1"] = out[0]
        fit_data.attrs["T1_stdd"] = stdd[0]
        fit_data.attrs["M_0"] = out[1]
        fit_data.attrs["M_inf"] = out[2]

    elif type == "T2":

        if stretched:
            x0 = [input_data[0], 1.0, 1.0]
            out, cov = curve_fit(
                dnpMath.t2_function, x_axis, input_data, x0, method="lm"
            )
            stdd = _np.sqrt(_np.diag(cov))
            fit = dnpMath.t2_function(new_axis, out[0], out[1], out[2])
        else:
            x0 = [input_data[0], 1.0]
            out, cov = curve_fit(
                dnpMath.t2_function, x_axis, input_data, x0, method="lm"
            )
            stdd = _np.sqrt(_np.diag(cov))
            fit = dnpMath.t2_function(new_axis, out[0], out[1])

        fit_data = dnpdata(fit, [new_axis], [ind_dim])
        fit_data.attrs["T2"] = out[1]
        fit_data.attrs["T2_stdd"] = stdd[1]
        fit_data.attrs["M_0"] = out[0]
        if stretched:
            fit_data.attrs["p"] = out[2]

    elif type == "mono":

        x0 = [input_data[-1], 1.0, 100]
        out, cov = curve_fit(dnpMath.monoexp_fit, x_axis, input_data, x0, method="lm")
        stdd = _np.sqrt(_np.diag(cov))
        fit = dnpMath.monoexp_fit(new_axis, out[0], out[1], out[2])

        fit_data = dnpdata(fit, [new_axis], [ind_dim])
        fit_data.attrs["tau"] = out[2]
        fit_data.attrs["tau_stdd"] = stdd[2]
        fit_data.attrs["C1"] = out[0]
        fit_data.attrs["C2"] = out[1]

    elif type == "bi":

        x0 = [input_data[-1], 1.0, 100, 1.0, 100]
        out, cov = curve_fit(dnpMath.biexp_fit, x_axis, input_data, x0, method="lm")
        stdd = _np.sqrt(_np.diag(cov))
        fit = dnpMath.biexp_fit(new_axis, out[0], out[1], out[2], out[3], out[4])

        fit_data = dnpdata(fit, [new_axis], [ind_dim])
        fit_data.attrs["tau1"] = out[2]
        fit_data.attrs["tau1_stdd"] = stdd[2]
        fit_data.attrs["tau2"] = out[4]
        fit_data.attrs["tau2_stdd"] = stdd[4]
        fit_data.attrs["C1"] = out[0]
        fit_data.attrs["C2"] = out[1]
        fit_data.attrs["C3"] = out[3]

    else:
        raise TypeError("Invalid fit type")

    if isDict:
        all_data["fit"] = fit_data
    else:
        return fit_data


def enhancement_fit(dataDict):
    """Fits enhancement curves to return Emax and power and one half maximum saturation

    .. math::

        f(p) = \mathrm{E_{max} * p} / (\mathrm{p_{1/2}} + \mathrm{p})

    Args:
        all_data (dnpdata, dict): data container

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

    if "enhancements" not in all_data.keys():
        raise TypeError("please use dnpNMR.calculate_enhancement() first")

    power_axes = all_data["enhancements"].coords["power"]

    input_data = _np.real(all_data["enhancements"].values)

    x0 = [input_data[-1], 0.1]

    out, cov = curve_fit(
        dnpMath.buildup_function, power_axes, input_data, x0, method="lm"
    )
    stdd = _np.sqrt(_np.diag(cov))

    fit = dnpMath.buildup_function(power_axes, out[0], out[1])

    fit_data = dnpdata(fit, [power_axes], ["power"])
    fit_data.attrs["E_max"] = out[0]
    fit_data.attrs["E_max_stdd"] = stdd[0]
    fit_data.attrs["power_half"] = out[1]
    fit_data.attrs["power_half_stdd"] = stdd[1]

    if isDict:
        dataDict["fit"] = fit_data
    else:
        return fit_data
