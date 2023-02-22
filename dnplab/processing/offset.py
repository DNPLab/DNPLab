import numpy as _np


def remove_background(data, dim="t2", deg=0, regions=None):
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

    proc_parameters = {
        "dim": dim,
        "deg": deg,
        "regions": regions,
    }

    fit = background(data, dim=dim, deg=deg, regions=regions)
    data = data - fit

    proc_attr_name = "remove_background"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data


def background(data, dim="t2", deg=0, regions=None):
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

    for ix in range(out.shape[1]):
        p = _np.polyfit(coord[fit_points], out.values[:, ix][fit_points], deg=deg)
        fit = _np.polyval(p, coord)
        out.values[:, ix] = fit

    out.fold()

    return out
