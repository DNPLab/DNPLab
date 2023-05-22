from warnings import warn

import numpy as _np

def average(data, dim = "Average"):
    """Average a dimension

    Args:
        data (object) : DNPData object
        dim (str) : Dimension to average

    Returns:
        DNPData: Averaged data

    Examples:

        >>> data_averaged = dnp.average(data)
    """

    temp = data.copy()
    index = temp.index(dim)
    temp.values = temp.values.sum(index)/_np.size(temp.coords[dim])
    
    temp.coords.pop(dim)

    return temp