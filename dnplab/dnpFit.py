from . import return_data, dnpdata, dnpdata_collection
from . import dnpMath
import numpy as _np
from scipy.optimize import curve_fit


def fit(f, data, dim = None, p0 = None, sigma = None, absolute_sigma = False, check_finite = True, bounds = (-1*_np.inf, _np.inf), method = None, jac = None, **kwargs):
    """ Fit data to generic function and return the fitting parameters and fit
    """

    '''
    proc_parameters = {"dim": dim}
    original_order = data.dims  # Original order of dims
    data.reorder([dim])  # Move dim to first dimension
    all_values = data.values  # Extract Data Values for alignment
    '''



    ydata = data.values
    xdata = data.coords[dim]
    
    out = curve_fit(f, xdata, ydata, p0 = p0, sigma = sigma, absolute_sigma = absolute_sigma, check_finite = check_finite, bounds = bounds, method = method, jac = jac, **kwargs)
    
    fit = f(xdata, *out[0])

    return fit

def exponential_fit(
    all_data,
    type="mono",
    stretched=False,
    bounds=None,
    p0=None,
    dim="t1",
    ws_key="integrals",
):
    """Fits various forms of exponential functions

    .. math::

        f(t) &= M_{0} - M_{\infty} e^{-t/T_{1}} \\
        
             &= M_{0} e^{(-2(t/T_{2})^{p}} \\
        
             &= C1 + C2 e^{-t/tau} \\
        
             &= C1 + C2 e^{-t/tau1} + C3 e^{-t/tau2}


    Args:
        all_data (dnpdata, dict): data container, after integration with dnpTools.integrate
        
    +---------------+-------+-------------+-----------------------------------------------------------------------------------+
    | parameter     | type  | default     | description                                                                       |
    +===============+=======+=============+===================================================================================+
    | type          | str   | "mono"      | "T1" for inversion recovery fit, "T2" for stretched exponential, "mono", or "bi"  |
    +---------------+-------+-------------+-----------------------------------------------------------------------------------+
    | stretch       | bool  | False       | if False "p" is set to 1, if True "p" is a fit parameter                          |
    +---------------+-------+-------------+-----------------------------------------------------------------------------------+
    | bounds        | tuple | None        | bounds on fit parameters                                                          |
    +---------------+-------+-------------+-----------------------------------------------------------------------------------+
    | p0            | list  | None        | initial guesses for fit parameters                                                |
    +---------------+-------+-------------+-----------------------------------------------------------------------------------+
    | dim           | str   | "t1"        | direct dimension                                                                  |
    +---------------+-------+-------------+-----------------------------------------------------------------------------------+
    | ws_key        | str   | "integrals" | if False "p" is set to 1, if True "p" is a fit parameter                          |
    +---------------+-------+-------------+-----------------------------------------------------------------------------------+


    Returns:
        dnpdata_collection or dnpdata: Processed data in container, updated with fit data
        attributes: "T1" value and "T1_stdd" standard deviation for type="T1", "T2" value and "T2_stdd" standard deviation for type="T2", "tau" and "tau_stdd" for type="mono", or "tau1", "tau1_stdd", "tau2", and "tau2_stdd" for type="bi"
    """

    data, isDict = return_data(
        all_data,
    )

    if isDict:
        x_axis = all_data[ws_key].coords[dim]
        input_data = all_data[ws_key].real.values
    else:
        x_axis = data.coords[dim]
        input_data = data.real.values

    new_axis = _np.r_[_np.min(x_axis) : _np.max(x_axis) : 100j]

    if type == "T1":
        if p0 is None:
            x0 = [1.0, input_data[-1], input_data[-1]]
        elif isinstance(p0, (list, tuple)) and len(p0) == 3:
            x0 = p0
        else:
            raise TypeError(
                "p0 must be a list or tuple of length=3, see the T1 function"
            )

        if bounds:
            out, cov = curve_fit(
                dnpMath.t1_function, x_axis, input_data, x0, bounds=bounds, method="trf"
            )
        else:
            out, cov = curve_fit(
                dnpMath.t1_function, x_axis, input_data, x0, method="lm"
            )

        stdd = _np.sqrt(_np.diag(cov))
        fit = dnpMath.t1_function(new_axis, out[0], out[1], out[2])

        fit_data = dnpdata(fit, [new_axis], [dim])
        fit_data.attrs["T1"] = out[0]
        fit_data.attrs["T1_stdd"] = stdd[0]
        fit_data.attrs["M_0"] = out[1]
        fit_data.attrs["M_inf"] = out[2]

    elif type == "T2":
        if p0 is None:
            x0 = [input_data[0], 1.0, 1.0]
        elif isinstance(p0, (list, tuple)) and len(p0) == 3:
            x0 = p0
        else:
            raise TypeError(
                "p0 must be a list or tuple of length=3, see the T2 function"
            )

        if stretched:
            if bounds:
                out, cov = curve_fit(
                    dnpMath.t2_function,
                    x_axis,
                    input_data,
                    x0,
                    bounds=bounds,
                    method="trf",
                )
            else:
                out, cov = curve_fit(
                    dnpMath.t2_function, x_axis, input_data, x0, method="lm"
                )
        else:
            if bounds is None:
                bounds = (
                    [float("-inf"), float("-inf"), 0.99999],
                    [float("inf"), float("inf"), 1.00001],
                )

            out, cov = curve_fit(
                dnpMath.t2_function,
                x_axis,
                input_data,
                x0,
                bounds=bounds,
                method="trf",
            )
            out[2] = 1.0

        stdd = _np.sqrt(_np.diag(cov))
        fit = dnpMath.t2_function(new_axis, out[0], out[1], out[2])

        fit_data = dnpdata(fit, [new_axis], [dim])
        fit_data.attrs["T2"] = out[1]
        fit_data.attrs["T2_stdd"] = stdd[1]
        fit_data.attrs["M_0"] = out[0]
        fit_data.attrs["p"] = out[2]

    elif type == "mono":
        if p0 is None:
            x0 = [input_data[-1], 1.0, 100]
        elif isinstance(p0, (list, tuple)) and len(p0) == 3:
            x0 = p0
        else:
            raise TypeError(
                "p0 must be a list or tuple of length=3, see the mono-exponential function"
            )

        if bounds:
            out, cov = curve_fit(
                dnpMath.monoexp_fit, x_axis, input_data, x0, bounds=bounds, method="trf"
            )
        else:
            out, cov = curve_fit(
                dnpMath.monoexp_fit, x_axis, input_data, x0, method="lm"
            )

        stdd = _np.sqrt(_np.diag(cov))
        fit = dnpMath.monoexp_fit(new_axis, out[0], out[1], out[2])

        fit_data = dnpdata(fit, [new_axis], [dim])
        fit_data.attrs["tau"] = out[2]
        fit_data.attrs["tau_stdd"] = stdd[2]
        fit_data.attrs["C1"] = out[0]
        fit_data.attrs["C2"] = out[1]

    elif type == "bi":
        if p0 is None:
            x0 = [input_data[-1], 1.0, 100, 1.0, 100]
        elif isinstance(p0, (list, tuple)) and len(p0) == 5:
            x0 = p0
        else:
            raise TypeError(
                "p0 must be a list or tuple of length=5, see the bi-exponential function"
            )

        if bounds:
            out, cov = curve_fit(
                dnpMath.biexp_fit, x_axis, input_data, x0, bounds=bounds, method="trf"
            )
        else:
            out, cov = curve_fit(dnpMath.biexp_fit, x_axis, input_data, x0, method="lm")

        stdd = _np.sqrt(_np.diag(cov))
        fit = dnpMath.biexp_fit(new_axis, out[0], out[1], out[2], out[3], out[4])

        fit_data = dnpdata(fit, [new_axis], [dim])
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


def enhancement_fit(all_data, dim="power", bounds=None, p0=None):
    """Fits enhancement curves to return Emax and power and one half maximum saturation

    .. math::

        f(p) = \mathrm{E_{max} * p} / (\mathrm{p_{1/2}} + \mathrm{p})

    Args:
        all_data (dnpdata, dict): data container
        dim (str): name of power dimension
        bounds (tuple): bounds on fit parameters
        p0 (list): initial guesses for fit parameters

    Returns:
        fit_data (dnpdata): Processed data in container, updated with fit data
        attrs (dict): Emax, Emax standard deviation, p_one_half, and p_one_half standard deviation

    Example::

        ### INSERT importing and processing ###
        dnplab.dnpNMR.integrate(workspace)

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

    if isDict:
        if "enhancements" not in all_data.keys():
            raise TypeError("please use dnpNMR.calculate_enhancement() first")
        else:
            power_axes = all_data["enhancements"].coords[dim]
            input_data = _np.real(all_data["enhancements"].values)
    else:
        power_axes = data.coords[dim]
        input_data = _np.real(data.values)

    if p0 is None:
        x0 = [input_data[-1], 0.1]
    elif isinstance(p0, (list, tuple)) and len(p0) == 2:
        x0 = p0
    else:
        raise TypeError("p0 must be a list or tuple of length=2, see the emax function")

    if bounds:
        out, cov = curve_fit(
            dnpMath.buildup_function,
            power_axes,
            input_data,
            x0,
            bounds=bounds,
            method="trf",
        )
    else:
        out, cov = curve_fit(
            dnpMath.buildup_function, power_axes, input_data, x0, method="lm"
        )

    stdd = _np.sqrt(_np.diag(cov))
    fit = dnpMath.buildup_function(power_axes, out[0], out[1])

    fit_data = dnpdata(fit, [power_axes], [dim])
    fit_data.attrs["E_max"] = out[0]
    fit_data.attrs["E_max_stdd"] = stdd[0]
    fit_data.attrs["power_half"] = out[1]
    fit_data.attrs["power_half_stdd"] = stdd[1]

    if isDict:
        dataDict["fit"] = fit_data
    else:
        return fit_data
