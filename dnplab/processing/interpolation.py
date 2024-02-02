"""
DNPLab uses numpy.interp function to interpolate DNPData object

Author: Yen-Chun Huang
"""

import numpy as _np
from ..core.data import DNPData


def interp(data, dim, new_coord, left=None, right=None):
    """
    Interpolate DNPData object

    Args:
        data (DNPData): Data object
        dim (str): Dimension to interpolate
        new_coord (list or numpy.arrays): New collection of numpy.ndarrays defining the axes
        left (optional: float or complex): Corresponding to data value, see numpy.interp for more details
        right (optional: float or complex): Corresponding to data value, see numpy.interp for more details

    Returns:
        data (DNPData): interpolated data object

    Example:
        data = dnp.interp(data, 'f2', new_coords = np.r_[-10:10:1000j])

    """

    out = data.copy()

    proc_parameters = {"dim": dim, "new_coord": new_coord, "left": left, "right": right}

    if len(_np.shape(new_coord)) != 1:
        raise ValueError("The input coord can only be one dimension")

    if isinstance(new_coord, list):
        new_coord = _np.array(new_coord)

    out.unfold(dim)
    coord = out.coords[dim]
    new_length = len(new_coord)

    folded_order = out.attrs["folded_order"]
    dim_index = folded_order.index(dim)

    folded_shape = out.attrs["folded_shape"]
    out.attrs["folded_shape"] = tuple(
        [
            folded_shape[i] if i != dim_index else new_length
            for i in range(len(folded_shape))
        ]
    )

    all_values = out.values.T
    new_values = []
    for values in all_values:
        new_values.append(_np.interp(new_coord, coord, values, left, right))

    new_values = _np.array(new_values).T

    out.values = new_values
    out.coords[dim] = new_coord

    out.fold()

    proc_attr_name = "interp"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out
