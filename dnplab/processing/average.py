import numpy as _np


def average(data, dim="Average"):
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
    temp = _np.mean(temp, axis=dim)

    return temp
