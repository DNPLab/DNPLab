import numpy as _np
from ..fitting import *
from ..math import relaxation

def remove_background(data, dim="t2", deg=0, regions=None, f:callable = None, **kwargs):
    """Remove polynomial background from data

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform background fit
        deg (int): Polynomial degree
        regions (None, list): Background regions, by default the entire region is used to calculate the background correction. Regions can be specified as a list of tuples [(min, max), ...]

    Returns:
        data (DNPData): Background corrected data

    Examples:

        0th-order background removal (DC offset)

        >>> data = dnp.remove_background(data)

    """

    out = data.copy()

    proc_parameters = {
        "dim": dim,
        "deg": deg,
        "regions": regions,
    }
    proc_parameters['f'] = f.__name__ if f else None 
    proc_parameters = {**proc_parameters, **kwargs}
    
    bg = background(out, dim=dim, deg=deg, regions=regions, f = f, **kwargs)
    out = out - bg

    proc_attr_name = "remove_background"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def background(data, dim="t2", deg=0, regions=None, f:callable = None, **kwargs):
    """Remove background from data

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform background fit
        deg (int): Polynomial degree
        regions (None, list): Background regions, by default entire region is background corrected. Regions can be specified as a list of tuples [(min, max), ...]

    Returns:
        DNPData: Background fit
    """

    out = data.copy()

    proc_parameters = {
        "dim": dim,
        "deg": deg,
        "regions": regions,
    }
    proc_parameters['f'] = f.__name__ if f else None 
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
    if not f:
        for ix in range(out.shape[1]):
            p = _np.polyfit(coord[fit_points], out.values[:, ix][fit_points], deg=deg)
            bg = _np.polyval(p, coord)
            out.values[:, ix] = bg

        out.fold()
    else:
        out = fit(f, data, dim = dim, **kwargs)['fit']

    proc_attr_name = "background"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out
