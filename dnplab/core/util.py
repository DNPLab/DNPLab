from .data import DNPData
import numpy as _np


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
    dnplab_attrs = data_list[0].dnplab_attrs

    values = _np.stack(values_list, axis=-1)

    dims.append(dim)

    if coord is None:
        coords.append(_np.arange(len(data_list)))
    else:
        coords.append(coord)

    data = DNPData(values, dims, coords, attrs, dnplab_attrs)

    return data
