from .data import DNPData
import numpy as _np
from .base import _SPECIAL_NP_HANDLED

def implements_np(np_function):
    "register an numpy function for special handling in SPECIAL_NO_HANDLED"
    def decorator(someFunction):
        _SPECIAL_NP_HANDLED[np_function]=someFunction
        return someFunction
    return decorator

def concat(data_list, dim, coord=None):
    """Concatenates list of data objects down another dimension

    args:
        data_list (list): List of DNPData objects to concatentate
        dim (str): new dimension name
        coord: coords for new dimension

    Returns:
        data (DNPData): concatenated data object

    """

    shape = data_list[0].shape
    values_list = [data.values for data in data_list]

    for values in values_list:
        this_shape = values.shape
        if this_shape != shape:
            raise IndexError(
                "Cannot concatenate data objects. Array shapes do not match.",
                this_shape,
                shape,
            )

    dims = data_list[0].dims
    coords = data_list[0].coords.coords
    attrs = data_list[0].attrs

    values = _np.stack(values_list, axis=-1)

    dims.append(dim)

    if coord is None:
        coords.append(_np.arange(len(data_list)))
    else:
        coords.append(coord)

    data = DNPData(values, dims, coords, attrs)

    return data
