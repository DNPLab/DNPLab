import numpy as np
from scipy.optimize import curve_fit
from ..core.data import DNPData


def fit(
    f,
    data,
    dim=None,
    p0=None,
    sigma=None,
    absolute_sigma=False,
    check_finite=True,
    bounds=(-1 * np.inf, np.inf),
    method=None,
    jac=None,
    **kwargs
):
    """
    proc_parameters = {"dim": dim}
    original_order = data.dims  # Original order of dims
    data.reorder([dim])  # Move dim to first dimension
    all_values = data.values  # Extract Data Values for alignment
    """

    fit = data.copy()
    fit.unfold(dim)

    ydata = fit.values
    xdata = fit.coords[dim]

    popt_list = []
    perr_list = []
    for ix in range(fit.shape[1]):
        ydata = fit.values[:, ix]
        out = curve_fit(
            f,
            xdata,
            ydata,
            p0=p0,
            sigma=sigma,
            absolute_sigma=absolute_sigma,
            check_finite=check_finite,
            bounds=bounds,
            method=method,
            jac=jac,
            **kwargs
        )
        fit_values = f(xdata, *out[0])
        fit.values[:, ix] = fit_values
        popt = out[0]
        pcov = out[1]
        perr = np.sqrt(np.diag(pcov))
        popt_list.append(popt)
        perr_list.append(perr)

    p_shape = list(fit.attrs["folded_shape"])
    p_shape[0] = len(p0)
    popt_array = np.array(popt_list).T.reshape(p_shape)
    perr_array = np.array(perr_list).T.reshape(p_shape)

    fit.fold()

    pdims = fit.dims
    pdims[0] = "popt"

    pcoords = [x for x in fit.coords]

    pcoords[0] = np.array(range(0, len(p0)))

    popt_data = DNPData(popt_array, pdims, pcoords)
    perr_data = DNPData(perr_array, pdims, pcoords)

    out = {
        "fit": fit,
        "popt": popt_data,
        "err": perr_data,
    }

    return out
