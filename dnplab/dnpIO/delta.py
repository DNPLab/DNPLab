import numpy as np
from struct import unpack
from dnplab import dnpdata


def import_delta(path, sep_phase_cycle=False):
    """
    Import Delta data and return dnpdata object

    Args:
        path (str) : Path to .jdf file

    Returns:
        delta_data (object) : dnpdata object containing Delta data
    """
    # pars = import_delta_pars(path)
    values, coords, dims, attrs = import_delta_data(
        path, sep_phase_cycle=sep_phase_cycle
    )

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


def import_delta_data(path, sep_phase_cycle=False):
    """
    Import spectrum or spectra of Delta data

    Args:
        path (str) : Path to .jdf file

    Returns:
        y_data (ndarray) : spectrum or spectra if >1D
        abscissa (list) : coordinates of axes
        dims (list) : axes names
        params (dict) : dictionary of parameters
    """
    file_opened = open(path, "rb")
    file_contents = file_opened.read(1296)
    file_opened.close()

    params = {}
    exp_type = [unpack(">B", file_contents[24 + k : 25 + k])[0] for k in range(0, 8, 1)]
    num_dims = [
        unpack(">B", file_contents[12 + k : 13 + k])[0] for k in range(0, 1, 1)
    ][0]
    params["nmr_frequency"] = [
        unpack(">d", file_contents[1064 + k : 1072 + k])[0] for k in range(0, 64, 8)
    ][0]
    axes_units = [
        unpack(">B", file_contents[32 + k : 33 + k])[0] for k in range(0, 16, 1)
    ]

    if axes_units[1] == 13:
        params["units"] = "Hz"
    elif axes_units[1] == 26:
        params["units"] = "ppm"
    elif axes_units[1] == 28:
        params["units"] = "s"

    axes_coords_start = [
        unpack(">d", file_contents[272 + k : 280 + k])[0] for k in range(0, 64, 8)
    ]
    axes_coords_stop = [
        unpack(">d", file_contents[336 + k : 344 + k])[0] for k in range(0, 64, 8)
    ]
    num_pts = [
        unpack(">I", file_contents[176 + k : 180 + k])[0] for k in range(0, 32, 4)
    ]

    file_opened = open(path, "rb")
    if num_dims == 1:
        num_pts = int(num_pts[0])
        load_pts = num_pts
        file_opened.seek(load_pts)
        if exp_type[0] == 3:
            data = np.fromfile(file_opened, "<d", load_pts * 2)
            y_data = np.split(data, 2)[0] - 1j * np.split(data, 2)[1]
        else:
            y_data = np.fromfile(file_opened, "<d", load_pts)
        abscissa = [np.linspace(axes_coords_start[0], axes_coords_stop[0], num_pts)]
        dims = ["t2"]

    elif num_dims == 2:
        num_pts_x = int(num_pts[0])
        num_pts_y = int(num_pts[1])
        load_pts = num_pts_x * num_pts_y
        file_opened.seek(load_pts)
        data = np.fromfile(file_opened, "<d", load_pts)
        abscissa = []
        if exp_type[0] == 4:
            data_folded = np.split(data, 2)[0] - 1j * np.split(data, 2)[1]
            data_shaped = np.reshape(data_folded, [int(num_pts_x / 2), int(num_pts_y)])
            abscissa.append(
                np.linspace(
                    axes_coords_start[0], axes_coords_stop[0], int(num_pts_x / 8)
                )
            )
            abscissa.append(
                np.linspace(axes_coords_start[1], axes_coords_stop[1], num_pts_y)
            )
            dims = ["t2", "t1"]
            y_data = []
            if sep_phase_cycle:
                for ix in range(num_pts_y):
                    y_data.append(np.split(data_shaped[:, ix], 4))
            else:
                for ix in range(num_pts_y):
                    y_data.append(np.sum(np.split(data_shaped[:, ix], 4), axis=0))
    else:
        raise ValueError("Only 1D or 2D are supported")

    file_opened.close()

    return y_data, abscissa, dims, params
