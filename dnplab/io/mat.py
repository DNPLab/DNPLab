"""
Save or load mat file 

Author: Yen-Chun Huang, Timothy Keller

"""

from scipy.io import savemat, loadmat
import numpy as np
from numpy import array, uint32
import os
import dnplab as _dnp


def save_mat(data, path, **kwargs):
    """Save DNPLab data object as Matlab .mat file

    Args:
        data (DNPData): Dictionary of data
        filename (str): directory to save file

    """

    matlab_dict = {}

    matlab_dict["values"] = data.values
    matlab_dict["dims"] = data.dims
    matlab_dict["coords"] = data.coords.coords
    matlab_dict["attrs"] = str(data.attrs)
    matlab_dict["dnplab_attrs"] = str(data.dnplab_attrs)
    matlab_dict["proc_attrs"] = str(data.proc_attrs).replace(
        "None", "'__PYTHON_NONE__'"
    )

    savemat(path, matlab_dict, **kwargs)


def import_mat(path):
    """
    Import mat file and return DNPData object

    Args:
        path (str):                 Path to .mat file

    Returns:
        data (DNPData):             DNPData objects
    """

    if path[-1] == os.sep:
        path = path[:-1]

    if ".mat" not in path:
        raise TypeError("Incorrect file type, must be .mat")

    data = loadmat("test.mat")
    values = _get_values(data)
    shape = np.shape(values)
    dims = _get_dims(data, shape)
    coords = _get_coords(data, shape)
    # attrs = _get_attrs(data, "attrs")
    # dnplab_attrs = _get_attrs(data, "dnplab_attrs")
    # proc_attrs = eval(
    #     str(_get_attrs(data, "proc_attrs")).replace("__PYTHON_NONE__", "None")
    # )
    attrs = {}
    dnplab_attrs = {}
    proc_attrs = []

    data = _dnp.DNPData(
        values=values,
        dims=dims,
        coords=coords,
        attrs=attrs,
        dnplab_attrs=dnplab_attrs,
        proc_attrs=proc_attrs,
    )

    return data


def _get_values(data):
    """
    Return data if values in 'values' in data

    Args:
        data (dict)                 Mat data

    Returns:
        values (numpy.array)        Numpy array of data
    """
    values = data["values"]
    return values


def _get_dims(data, shape):
    """
    Get dims in data

    Args:
        data (dict)                 Mat data
        shape (tuple)               Data shape

    Returns:
        dims (list)                 List of dims
    """

    dims = data.get("dims", [])
    if len(dims):
        dims = dims.tolist()
    else:
        for i in range(len(shape)):
            dims.append("X%i" % (i + 1))
    return dims


def _get_coords(data, shapes):
    """
    Get dims in data

    Args:
        data (dict)                 Mat data
        shape (tuple)               Data shape

    Returns:
        coords (list)                 List of dims
    """

    coords = data.get("coords", [])
    if len(coords):
        coords = coords.tolist()
    else:
        for i in range(len(shapes)):
            coords.append(list(range(shapes[i])))
    return coords


# def _get_attrs(data, key="attrs"):
#     """
#     Get attrs in data

#     Args:
#         data (dict)                 Mat data
#         key (str)                   Key to attributes

#     Returns:
#         attrs (dict, list)          Saved attributes in .mat file
#     """

#     attrs = data.get(key, {})
#     if attrs:
#         attrs = eval(attrs[0])

#     return attrs
