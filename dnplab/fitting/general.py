import numpy as _np
from scipy.optimize import curve_fit


def fit(
    f,
    data,
    dim=None,
    p0=None,
    sigma=None,
    absolute_sigma=False,
    check_finite=True,
    bounds=(-1 * _np.inf, _np.inf),
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

    fit.fold()

    return fit
