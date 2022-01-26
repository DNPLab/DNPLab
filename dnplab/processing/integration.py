import numpy as np
from ..core.data import DNPData
from ..core.data import concat

from scipy.integrate import trapz

def cumtrapz(

):
    """ Cummulative sum down given dimension

    """

def integrate(
    data,
    dim="f2",
    regions = None
):
    """Integrate data down given dimension

    Args:
        all_data (dnpdata,dict): Data container

    +------------------+---------------+----------+-------------------------------+
    | parameter        | type          | default  | description                   |
    +==================+===============+==========+===============================+
    | dim              | str           | 'f2'     | dimension to integrate        |
    +------------------+---------------+----------+-------------------------------+
    | regions          | list          | None     | Regions for integration       |
    +------------------+---------------+----------+-------------------------------+

    Returns:
        dnpdata: integrals of data
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
