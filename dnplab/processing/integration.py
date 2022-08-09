import numpy as np
from ..core.data import DNPData
from ..core.util import concat

from scipy.integrate import trapz, cumtrapz


def cumulative_integrate(data, dim="f2", regions=None):
    """Cumulative integration

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform cumulative integration
        regions (None, list): List of tuples to specify range of integration [(min, max), ...]

    Returns:
        data: cumulative sum of data

    Examples:
        Example showing cumulative integration of lorentzian function

        >>> import numpy as np
        >>> from matplotlib.pylab import *
        >>> import dnplab as dnp
        >>> x = np.r_[-10:10:1000j]
        >>> y = dnp.math.lineshape.lorentzian(x,0,1)
        >>> data = dnp.DNPData(y, ['f2'], [x])
        >>> data_int = dnp.cumulative_integrate(data)
        >>> figure()
        >>> dnp.plot(data)
        >>> dnp.plot(data_int)
        >>> show()


    """

    data = data.copy()

    if regions == None:
        index = data.index(dim)
        data.values = cumtrapz(data.values, data.coords[dim], axis=index, initial=0)
        return data
    else:
        data_list = []
        for region in regions:
            data_list.append(cumulative_integrate(data[dim, region], dim))

        return data_list


def integrate(data, dim="f2", regions=None):
    """Integrate data down given dimension

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform integration
        regions (None, list): List of tuples, by default entire dimension is integrated

    Returns:
        data: integrals of data
    """

    data = data.copy()
    data.attrs["experiment_type"] = "integrals"

    index = data.index(dim)
    if regions == None:
        data.values = trapz(data.values, data.coords[dim], axis=index)
        data.coords.pop(dim)

    else:
        data_list = []
        for region in regions:
            data_list.append(integrate(data[dim, region], dim))

        x = np.array(list(range(len(data_list))))
        dim_name = "integrals"
        data = concat(data_list, dim_name, coord=x)

    return data
