import numpy as _np
from ..fitting import *
from ..math import relaxation


def remove_background(
    data, dim="t2", deg=0, regions=None, func: callable = None, **kwargs
):
    """Remove polynomial background from data

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform background fit
        deg (int): Polynomial degree
        regions (None, list): Background regions, by default the entire region is used to calculate the background correction. Regions can be specified as a list of tuples [(min, max), ...]
        func (optional callable): The fitting function to fit the background
        **kwargs: arguments for fitting function

    Returns:
        data (DNPData): Background corrected data

    Examples:

        0th-order background removal (DC offset)

        >>> data = dnp.remove_background(data)


        Background removal with a given fit function

        >>> data = dnp.remove_background(data, dim = 'tau', func= dnp.relaxation.general_exp, p0=(1,-1,900))

    """

    out = data.copy()

    proc_parameters = {
        "dim": dim,
        "deg": deg,
        "regions": regions,
    }
    proc_parameters["func"] = func.__name__ if func else None
    proc_parameters = {**proc_parameters, **kwargs}

    bg = background(out, dim=dim, deg=deg, regions=regions, func=func, **kwargs)
    out = out - bg

    proc_attr_name = "remove_background"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def background(data, dim="t2", deg=0, regions=None, func: callable = None, **kwargs):
    """Remove background from data

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform background fit
        deg (int): Polynomial degree
        regions (None, list): Background regions, by default entire region is background corrected. Regions can be specified as a list of tuples [(min, max), ...]
        func (optional callable): The fitting function to fit the background
        **kwargs: arguments for fitting function

    Returns:
        DNPData: Background fit

    Examples:

        0th-order background fit (DC offset)

        >>> bg = dnp.background(data)


        Background with a given fit function

        >>> bg = dnp.background(data, dim = 'tau', func= dnp.relaxation.general_exp, p0=(1,-1,900))

    """

    out = data.copy()

    proc_parameters = {
        "dim": dim,
        "deg": deg,
        "regions": regions,
    }
    proc_parameters["func"] = func.__name__ if func else None
    proc_parameters = {**proc_parameters, **kwargs}

    out.unfold(dim)

    coord = out.coords[dim]

    if regions == None:
        fit_points = [True for x in coord]
    else:
        fit_points = [False for x in coord]
        for region in regions:
            fit_points = [
                fit_points[ix] or ((coord[ix] >= region[0]) & (coord[ix] <= region[1]))
                for ix in range(len(coord))
            ]
    if not func:
        for ix in range(out.shape[1]):
            if _np.iscomplexobj(out.values[:, ix]):
                out_real = out.values[:, ix].real
                out_imag = out.values[:, ix].imag
                p_real = _np.polyfit(coord[fit_points], out_real[fit_points], deg=deg)
                p_imag = _np.polyfit(coord[fit_points], out_imag[fit_points], deg=deg)
                bg_real = _np.polyval(p_real, coord)
                bg_imag = _np.polyval(p_imag, coord)
                out.values[:, ix] = bg_real + 1j * bg_imag
            else:
                p = _np.polyfit(
                    coord[fit_points], out.values[:, ix][fit_points], deg=deg
                )
                bg = _np.polyval(p, coord)
                out.values[:, ix] = bg

        out.fold()
    else:
        if _np.iscomplexobj(data.values):
            out_real = fit(func, data.real, dim=dim, **kwargs)["fit"]
            out_imag = fit(func, data.imag, dim=dim, **kwargs)["fit"]
            out = out_real + 1j * out_imag
        else:
            out = fit(func, data, dim=dim, **kwargs)["fit"]

    proc_attr_name = "background"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out
