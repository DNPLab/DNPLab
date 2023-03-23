from .. import DNPData
import numpy as _np
import struct
import re


def import_tnmr(path):
    """Import tnmr data and return DNPData object

    Args:
        path (str) : Path to .jdf file

    Returns:
        tnmr_data (object) : DNPData object containing tnmr data
    """

    attrs = import_tnmr_pars(path)
    values, dims, coords = import_tnmr_data(path)

    tnmr_data = DNPData(values, dims, coords, attrs)

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
        version = f.read(8).decode("utf-8")

        section = None

        while section != "":
            section = f.read(4).decode("utf-8")
            section = str(section)

            if section == "TMAG":
                flag = bool(f.read(4))
                if flag:
                    bytes_to_read = f.read(4)
                    bytes_to_read = struct.unpack("<i", bytes_to_read)[0]

                    header = f.read(bytes_to_read)

                    ### Deal With Header Here ###

            elif section == "DATA":
                flag = bool(f.read(4))
                if flag:
                    bytes_to_read = f.read(4)
                    bytes_to_read = struct.unpack("<i", bytes_to_read)[0]

                    raw_data = f.read(bytes_to_read)

                    raw_data = struct.unpack("%if" % (bytes_to_read / 4), raw_data)

                    raw_data = _np.array(raw_data)

                    data = raw_data[::2] + 1j * raw_data[1::2]

            else:
                flag = bool(f.read(4))
                if flag:
                    bytes_to_read = f.read(4)
                    bytes_to_read = struct.unpack("<i", bytes_to_read)[0]

                    unsupported_bytes = f.read(bytes_to_read)

    abscissa = _np.array(range(0, len(data)))

    dims = ["t2"]
    coords = [abscissa]

    return data, dims, coords
