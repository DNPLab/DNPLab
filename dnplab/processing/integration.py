import numpy as np
from ..core.data import DNPData
from ..core.data import concat

from scipy.integrate import trapz, cumtrapz


def cumulative_integrate(
    data,
    dim = "f2",
):
    """ Cumulative integration

    Args:
        data (DNPData): Data object

    Returns:
        data: cumulative sum of data
    """

    data = data.copy()

    index = data.index(dim)
    data.values = cumtrapz(data.values, data.coords[dim], axis = index, initial = 0)

    return data

def integrate(
    data,
    dim="f2",
    regions = None
):
    """Integrate data down given dimension

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform integration
        regions (list): List of tuples, by default entire dimension is integrated

    +------------------+---------------+----------+-------------------------------+
    | parameter        | type          | default  | description                   |
    +==================+===============+==========+===============================+
    | dim              | str           | 'f2'     | dimension to integrate        |
    +------------------+---------------+----------+-------------------------------+
    | regions          | list          | None     | Regions for integration       |
    +------------------+---------------+----------+-------------------------------+

    Returns:
        data: integrals of data
    """

    data = data.copy()

    index = data.index(dim)
    if regions == None:
        data.values = trapz(data.values, data.coords[dim], axis = index)
        data = data.sum(dim)

    else:
        data_list = []
        for region in regions:
            data_list.append(integrate(data[dim, region], dim))
        
        x = np.array(list(range(len(data_list))))
        dim_name = 'integrals'
        data = concat(data_list, dim_name, coord = x)

    return data
