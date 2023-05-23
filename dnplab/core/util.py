from .data import DNPData
import numpy as _np
from .base import _SPECIAL_NP_HANDLED


def implements_np(np_function):
    "register an numpy function for special handling in SPECIAL_NO_HANDLED"

    def decorator(someFunction):
        _SPECIAL_NP_HANDLED[np_function] = someFunction
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
    dnplab_attrs = data_list[0].dnplab_attrs

    values = _np.stack(values_list, axis=-1)

    dims.append(dim)

    if coord is None:
        coords.append(_np.arange(len(data_list)))
    else:
        coords.append(coord)

    data = DNPData(values, dims, coords, attrs, dnplab_attrs)

    return data


def update_axis(data, new_dims, start_stop, spacing="lin", verbose=False):
    """Update axis

    Update dimensions (dims) and axis (coords) of a dnpDate object. The name of the dims will be replaced with the name giving in new_dims. The variable start_stop defines the values of the new coords. This can be either a tuple (start values, stop value) or a vector with values. If the start and stop value is provided, either a linear axis (spacing = "lin", default) or a logarithmically space (spacing = "log") will be created. The new axis will replace the coords in the dnpdata object.

    The function is currently implemented for 1D objects only.

    Args:
        data (DNPData): dnpData objetc
        new_dims (str): Name of the new dimension. If None the name will not be changed.
        start_stop(tuple or vector): Coords for new dimension
        spacing (str): "lin" for linear spaced axis or "log" for logarithmically spaced axis

    Returns:
        data (DNPData): concatenated data object
    """

    out = data.copy()

    if new_dims == None:
        new_dims = out.dims[0]

    out.rename(out.dims[0], new_dims)

    if verbose == True:
        print("New dims name:", new_dims)

    if len(start_stop) == 2:
        if spacing == "lin":
            out.coords[new_dims] = _np.linspace(
                start_stop[0], start_stop[1], len(out.values)
            )

        elif spacing == "log":
            out.coords[new_dims] = _np.logspace(
                start_stop[0], start_stop[1], len(out.values)
            )

    elif len(start_stop) > 2:
        out.coords[new_dims] = start_stop

    return out
