import numpy as np
from struct import unpack
from .. import dnpdata


def import_delta(path):
    """
    Import Delta data and return dnpdata object

    Args:
        path (str) : Path to .jdf file

    Returns:
        delta_data (object) : dnpdata object containing Delta data
    """

    pars = import_delta_pars(path)
    values, coords, dims, attrs = import_delta_data(path, pars)

    delta_data = dnpdata(values, coords, dims, attrs)

    return delta_data


def import_delta_pars(path):
    """
    Import parameter fields of Delta data

    Args:
        path (str) : Path to .jdf file

    Returns:
        params (dict) : dictionary of parameter fields and values
    """

    file_opened = open(path, "rb")
    file_contents = file_opened.readlines()
    file_opened.close()
    params = {}
    for ix in range(len(file_contents)):
        try:
            file_contents[ix] = str(file_contents[ix], "utf8")
            if "=" in file_contents[ix]:
                new_line = file_contents[ix].split("=")
                new_name = new_line[0].replace("/", "").strip()
                if "@" in new_name or new_name == "":
                    pass
                else:
                    params[new_name] = (
                        new_line[1]
                        .replace(">", "")
                        .replace(";", "")
                        .replace("?", "")
                        .strip()
                    )
        except UnicodeDecodeError:
            pass

    return params


def import_delta_data(path, params):
    """
    Import spectrum or spectra of Delta data

    Args:
        path (str) : Path to .jdf file
        params (dict) : dictionary of parameters

    Returns:
        y_data (ndarray) : spectrum or spectra if >1D
        abscissa (list) : coordinates of axes
        dims (list) : axes names
        params (dict) : updated dictionary of parameters
    """

    if not params:
        params = {}

    file_opened = open(path, "rb")
    file_contents = file_opened.read(1296)
    file_opened.close()

    num_dims = [
        unpack(">B", file_contents[12 + ix : 13 + ix])[0] for ix in range(0, 1, 1)
    ][0]

    params["nmr_frequency"] = [
        unpack(">d", file_contents[1064 + ix : 1072 + ix])[0] for ix in range(0, 64, 8)
    ][0]
    axes_units = [
        unpack(">B", file_contents[32 + ix : 33 + ix])[0] for ix in range(1, 16, 2)
    ][:num_dims]
    params["units"] = []
    for ix in range(num_dims):
        if axes_units[ix] == 1:
            params["units"].append("abundance")
        elif axes_units[ix] == 13:
            params["units"].append("Hz")
        elif axes_units[ix] == 26:
            params["units"].append("ppm")
        elif axes_units[ix] == 27:
            params["units"].append("rad")
        elif axes_units[ix] == 28:
            params["units"].append("s")
        else:
            params["units"].append("indexed")

    endian = [unpack(">B", file_contents[8 + ix : 9 + ix])[0] for ix in range(0, 1, 1)][
        0
    ]
    if endian == 0:
        endian = ">d"
    elif endian == 1:
        endian = "<d"
    else:
        raise UnicodeTranslateError("Failed to determine endianness")

    num_pts = [
        unpack(">I", file_contents[176 + ix : 180 + ix])[0] for ix in range(0, 32, 4)
    ]
    axis_type = [
        unpack(">B", file_contents[24 + ix : 25 + ix])[0] for ix in range(0, 8, 1)
    ][:num_dims]
    axes_start = [
        unpack(">d", file_contents[272 + ix : 280 + ix])[0] for ix in range(0, 64, 8)
    ][:num_dims]
    axes_stop = [
        unpack(">d", file_contents[336 + ix : 344 + ix])[0] for ix in range(0, 64, 8)
    ][:num_dims]
    abscissa = []
    for ix in range(num_dims):
        abscissa.append(np.linspace(axes_start[ix], axes_stop[ix], num_pts[ix]))

    data_start = [
        unpack(">I", file_contents[1284 + ix : 1288 + ix])[0] for ix in range(0, 4, 4)
    ][0]

    file_opened = open(path, "rb")
    file_opened.seek(data_start)
    if num_dims == 2 and axis_type[0] == 3 and axis_type[1] == 3:
        read_pts = np.prod(num_pts) * 4
    else:
        read_pts = np.prod(num_pts) * 2
    data = np.fromfile(file_opened, endian, read_pts)
    file_opened.close()

    if num_dims == 1:
        if axis_type[0] == 1:
            y_data = data
        elif axis_type[0] == 3 or axis_type[0] == 4:
            y_data = np.split(data, 2)[0] - 1j * np.split(data, 2)[1]
        else:
            raise TypeError("Data format not recognized")

        dims = ["t2"]

    elif num_dims == 2:
        if axis_type[0] == 4 or (axis_type[0] == 3 and axis_type[1] == 1):
            data_folded = np.split(data, 2)[0] - 1j * np.split(data, 2)[1]
            data_shaped = np.reshape(
                data_folded, [int(num_pts[0] / 4), int(num_pts[1] / 4), 4, 4], order="F"
            )
            y_data = np.concatenate(np.concatenate(data_shaped, 1), 1)
        elif axis_type[0] == 3 and axis_type[1] == 3:
            data_folded = [
                np.split(data, 4)[0] - 1j * np.split(data, 4)[1],
                np.split(data, 4)[2] - 1j * np.split(data, 4)[3],
            ]
            for idx in enumerate(data_folded):
                data_shaped[idx] = np.reshape(
                    data_folded[idx],
                    [int(num_pts[0] / 32), int(num_pts[1] / 32), 32, 32],
                    order="F",
                )
                y_data[idx] = np.concatenate(np.concatenate(data_shaped[idx], 1), 1)
        else:
            raise ValueError("Data format not recognized")

        dims = ["t2", "t1"]

    else:
        raise TypeError("Only 1D or 2D are supported")

    return y_data, abscissa, dims, params
