import numpy as _np
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
    """Integrate data along given dimension. If no region is given, the integral is calculated over the entire range.

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform integration. Default is "f2"
        regions (None, list): List of tuples defining the region to integrate

    Returns:
        data (DNPData): Integrals of data. If multiple regions are given the first value corresponds to the first region, the second value corresponds to the second region, etc.

    Examples:
        Integrated entire data region:

            >>> data = dnp.integrate(data)

        Integrate single peak/region:

            >>> data = dnp.integrate(data, regions=[(4, 5)])

        Integrate two regions:

            >>> data = dnp.integrate(data, regions=[(1.1, 2.1), (4.5, 4.9)])

    """
    data = data.copy()
    data.attrs["experiment_type"] = "integrals"

    index = data.index(dim)
    if regions == None:
        data.values = np.trapz(data.values, data.coords[dim], axis=index)
        data.coords.pop(dim)

        # if error_regions == None:
        #     data.error = np.zeros(data.shape)
        #     print("add errors")

        # else:
        #     signal = max(data.values)
        #     noise = np.trapz(data.)

    else:
        data_list = []
        for region in regions:
            data_list.append(integrate(data[dim, region], dim))

        x = _np.array(list(range(len(data_list))))
        dim_name = "integrals"
        data = concat(data_list, dim_name, coord=x)

    return data
