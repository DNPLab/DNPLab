import numpy as np
from ..core.data import DNPData
from ..core.util import concat

from scipy.integrate import trapz, cumtrapz

def _check_regions(regions):
    """Check integration regions"""
    # Regions is None
    if regions is None:
        return regions

    # Regions is tuple
    elif isinstance(regions, tuple):
        if len(regions) != 2:
            raise ValueError('If specifying single region, tuple must have two elements')
    
    # Regions is list
    elif isinstance(regions, list):
        for region in regions:
            if not isinstance(region, tuple):
                raise ValueError('Each region specified must be a tuple with two elements specifying the maximum and minimum of the integration region')

            elif len(region) != 2:
                raise ValueError('If specifying single region, tuple must have two elements')


    else:
        raise ValueError('Invalid regions')


    return regions

def cumulative_integrate(data, dim="f2", regions=None):
    """Cumulative integration

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform cumulative integration
        regions (None, list): List of tuples to specify range of integration [(min, max), ...]

    Returns:
        data: cumulative sum of data
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


    if regions == None:
        index = data.index(dim)
        data.values = trapz(data.values, data.coords[dim], axis=index)
        data.coords.pop(dim)
        return data
    
    else:
        regions = _check_regions(regions)
        if isinstance(regions, tuple):
            data = integrate(data[dim, regions], dim)
        else:
            data_list = []
            for region in regions:
                data_list.append(integrate(data[dim, region], dim))

            x = np.array(list(range(len(data_list))))
            dim_name = "integrals"
            data = concat(data_list, dim_name, coord=x)

    return data
