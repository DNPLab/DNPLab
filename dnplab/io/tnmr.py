from .. import DNPData
import numpy as _np
import struct
import re


def import_tnmr(path, squeeze=True):
    """Import tnmr data and return DNPData object

    Args:
        path (str) : Path to .jdf file
        squeeze (bool) : automatically remove length 1 dimensions

    Returns:
        tnmr_data (object) : DNPData object containing tnmr data
    """

    attrs = import_tnmr_pars(path)
    values, dims, coords = import_tnmr_data(path)

    tnmr_data = DNPData(values, dims, coords, attrs)

    if squeeze:
        tnmr_data.squeeze()

    return tnmr_data


def import_tnmr_pars(path):
    """Import parameter fields of tnmr data

    Args:
        path (str) : Path to .tnt file

    Returns:
        params (dict) : dictionary of parameter fields and values
    """

    params = {}

    with open(path, "rb") as f:
        params["version"] = f.read(8).decode("utf-8")

    return params


def import_tnmr_data(path):
    """Import spectrum or spectra of tnmr data

    Args:
        path (str) : Path to .tnt file

    Returns:
        data (ndarray) : spectrum or spectra if >1D
        abscissa (list) : coordinates of axes
        dims (list) : axes names
    """

    with open(path, "rb") as f:
        raw = f.read()

    tmag_attrs = {}
    tmag_attrs["version"] = str(raw[0:8])

    len_tecmag_struct = int.from_bytes(raw[16:20], "little")
    tecmag_struct = raw[20 : 20 + len_tecmag_struct]

    offset = 28 + len_tecmag_struct
    len_data = raw[offset : offset + 4]
    offset += 4

    len_data = int.from_bytes(len_data, "little")

    data = raw[offset : offset + (len_data)]

    data = _np.frombuffer(data, dtype="<f")

    data = data[0::2] + 1j * data[1::2]

    # Parse tecmag struct

    # points in x0, x1, x2, x3
    npts = _np.frombuffer(tecmag_struct[0:16], dtype="<i")

    # Array of points in each dimension (x0, x1, x2, x3)
    actual_npts = _np.frombuffer(tecmag_struct[16:32], dtype="<i")

    acq_pts = int.from_bytes(tecmag_struct[32:36], byteorder="little")
    scans = int.from_bytes(tecmag_struct[36:40], byteorder="little")
    # actual_scans = int.from_bytes(tecmag_struct[40:44], byteorder = 'little')
    # dummy_scans = int.from_bytes(tecmag_struct[44:48], byteorder = 'little')

    sw = struct.unpack("<4d", tecmag_struct[240:272])
    dwell_time = struct.unpack("<4d", tecmag_struct[272:304])
    # dwell_time = float.from_bytes(tecmag_struct[272:304], byteorder = 'little')

    data = data.reshape(npts, order="F")

    coords = []

    for ix, pts in enumerate(npts):
        coord = _np.r_[0:pts] * dwell_time[ix]
        coords.append(coord)

    dims = ["t2", "t1", "t3", "t4"]  # t2 dim is first

    return data, dims, coords
