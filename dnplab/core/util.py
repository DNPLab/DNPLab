from .data import DNPData
import numpy as _np
from .base import _SPECIAL_NP_HANDLED


def implements_np(np_function):
    "register an numpy function for special handling in SPECIAL_NO_HANDLED"

    def decorator(someFunction):
        _SPECIAL_NP_HANDLED[np_function] = someFunction
        return someFunction

    return decorator


def concat(data_list, dim, coord=None, casting="same_kind"):
    """Concatenates list of data objects down another dimension

    args:
        data_list (list): List of DNPData objects to concatentate
        dim (str): new dimension name
        coord: coords for new dimension

    Returns:
        data (DNPData): concatenated data object

    """

    shape_list = [data.shape for data in data_list]
    length_list = [len(shape) for shape in shape_list]
    values_list = [data.values for data in data_list]

    dims = data_list[0].dims
    coords = data_list[0].coords.coords

    attrs = data_list[0].attrs
    dnplab_attrs = data_list[0].dnplab_attrs

    match casting:
        case "same_kind":
            if not all(shape == shape_list[0] for shape in shape_list):
                raise IndexError(
                    "Cannot concatenate data objects. Array shapes do not match.",
                    shape_list,
                )

        case "unsafe":
            if not all(length == length_list[0] for length in length_list):
                raise IndexError(
                    "Cannot concatenate data objects. Data shapes length do not match.",
                    shape_list,
                )
            transposed_shape_list = _np.array(shape_list).T
            new_shape = []
            for i in range(len(transposed_shape_list)):
                new_shape.append(max(transposed_shape_list[i]))
            final_shape = tuple(_np.array(new_shape).T)

            # add nan data to the non-consistent dataset
            for i, (values, shape) in enumerate(zip(values_list, shape_list)):
                if shape != final_shape:
                    appending_shape = [1] * len(
                        final_shape
                    )  # generate new shape for appending
                    for pointer in range(len(final_shape)):
                        appending_length = final_shape[pointer] - shape[pointer]
                        if appending_length > 0:
                            appending_shape[pointer] = appending_length
                            e = _np.empty(tuple(appending_shape))
                            e[:] = _np.nan
                            values = _np.append(values, e, axis=pointer)
                        appending_shape[pointer] = final_shape[pointer]
                    values_list[i] = values

            # make coords
            for data in data_list:
                for i in range(len(coords)):
                    if len(coords[i]) < len(data.coords.coords[i]):
                        coords[i] = data.coords.coords[i]

        case _:
            raise ValueError("Currently 'same_kind' and 'unsafe' are available")

    values = _np.stack(values_list, axis=-1)
    dims.append(dim)

    if coord is None:
        coords.append(_np.arange(len(data_list)))
    else:
        coords.append(coord)

    data = DNPData(values, dims, coords, attrs, dnplab_attrs)

    return data


def update_axis(data, start_stop, dim=0, new_dims=0, spacing="lin", verbose=False):
    """Update axis

    Update dimensions (dims) and axis (coords) of a dnpDate object. The name of the dims will be replaced with the name giving in new_dims. The variable start_stop defines the values of the new coords. This can be either a tuple (start values, stop value) or a vector with values. If the start and stop value is provided, either a linear axis (spacing = "lin", default) or a logarithmically space (spacing = "log") will be created. The new axis will replace the coords in the dnpdata object.

    The function is currently implemented for 1D objects only.

    Args:
        data (DNPData):                 dnpData object
        start_stop(tuple or vector):    Coords for new dimension
        dim (int):                      Dimension to act on
        new_dims (str):                 Name of the new dimension. If None the name will not be changed.
        spacing (str):                  "lin" for linear spaced axis or "log" for logarithmically spaced axis

    Returns:
        data (DNPData): concatenated data object
    """

    out = data.copy()
    data_shape = _np.shape(out.values)

    if new_dims == None:
        new_dims = out.dims[dim]

    out.rename(out.dims[dim], new_dims)

    if verbose == True:
        print("New dims name:", new_dims)

    if len(start_stop) == 2:
        if spacing == "lin":
            out.coords[new_dims] = _np.linspace(
                start_stop[0], start_stop[1], data_shape[dim]
            )

        elif spacing == "log":
            out.coords[new_dims] = _np.logspace(
                start_stop[0], start_stop[1], data_shape[dim]
            )

    elif len(start_stop) > 2:
        out.coords[new_dims] = start_stop

    return out


def get_slice(data, dim, slice_index):
    """
    Get data slice of DNPData object

    Args:
        data (DNPData):     Input data object
        dim (str):          Selected dimension
        slice_index (int):  Index of slice to be returned

    Returns:
        data (DNPData):     DNPData object with selected slice

    """

    out = data.copy()
    out = out[dim, slice_index]

    return out
