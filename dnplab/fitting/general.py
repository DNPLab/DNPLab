import numpy as _np
from scipy.optimize import curve_fit
from ..core.data import DNPData


def fit(
    f,
    data,
    dim,
    p0,
    fit_points=50,
    sigma=None,
    absolute_sigma=False,
    check_finite=True,
    bounds=(-1 * _np.inf, _np.inf),
    method=None,
    jac=None,
    **kwargs
):
    """Fitting function for DNPData

    Args:
        f (func): Function used in scipy.curve_fit
        data (DNPData): Data for fit
        dim (str): Dimension to perform fit along
        p0 (tuple): Initial guess for fit
        fit_points (int): Number of points to use in the fit. If None (default), the number of points is the same as the data.
        kwargs: Additional parameters for scipy.curve_fit

    Returns:
        out (dict): Dictionary of fit, fitting parameters, and error
    """

    fit = data.copy()

    index = fit.index(dim)
    dims = fit.dims
    coords = list(fit.coords.copy())
    shape = list(fit.shape)

    coord = coords[index]
    if fit_points is not None:
        new_coord = _np.r_[_np.min(coord) : _np.max(coord) : 1j * fit_points]
    else:
        new_coord = coord

    shape[index] = len(new_coord)
    coords[index] = new_coord
    fit_out = DNPData(_np.zeros(shape), dims, coords)

    fit_out.unfold(dim)
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
        fit_values = f(new_coord, *out[0])
        fit_out.values[:, ix] = fit_values
        popt = out[0]
        pcov = out[1]
        perr = _np.sqrt(_np.diag(pcov))
        popt_list.append(popt)
        perr_list.append(perr)

    p_shape = list(fit.attrs["folded_shape"])
    p_shape[0] = len(p0)
    popt_array = _np.array(popt_list).T.reshape(p_shape)
    perr_array = _np.array(perr_list).T.reshape(p_shape)

    fit.fold()
    fit_out.fold()

    pdims = list(fit.dims)
    pdims[0] = "popt"

    pcoords = list(fit.coords)

    pcoords[0] = _np.array(range(0, len(p0)))

    popt_data = DNPData(popt_array, pdims, pcoords)
    perr_data = DNPData(perr_array, pdims, pcoords)

    out = {
        "fit": fit_out,
        "popt": popt_data,
        "err": perr_data,
    }

    return out
