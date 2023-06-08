import numpy as _np


def average(data, axis="Average"):
    """Average a dimension using numpy.mean

    Args:
        data (object) : DNPData object
        dim (str) : Dimension to average

    Returns:
        DNPData: Averaged data

    Examples:

        >>> data_averaged = dnp.average(data)
    """

    data = data.copy()

    data = _np.mean(data, axis=axis)  # it will automatically assign proc_attrs

    return data
