import numpy as _np
import re as _re

from .. import dnpdata as _dnpdata

import os as _os

_dspfvs_table_10 = {
    2: 44.7500,
    3: 33.5000,
    4: 66.6250,
    6: 59.0833,
    8: 68.5625,
    12: 60.3750,
    16: 69.5313,
    24: 61.0208,
    32: 70.0156,
    48: 61.3438,
    64: 70.2578,
    96: 61.5052,
    128: 70.3789,
    192: 61.5859,
    256: 70.4395,
    384: 61.6263,
    512: 70.4697,
    1024: 70.4849,
    1536: 61.6566,
    2048: 70.4924,
}

_dspfvs_table_11 = {
    2: 46.0000,
    3: 36.5000,
    4: 48.0000,
    6: 50.1667,
    8: 53.2500,
    12: 69.5000,
    16: 72.2500,
    24: 70.1667,
    32: 72.7500,
    48: 70.5000,
    64: 73.0000,
    96: 70.6667,
    128: 72.5000,
    192: 71.3333,
    256: 72.2500,
    384: 71.6667,
    512: 72.1250,
    1024: 72.0625,
    1536: 71.9167,
    2048: 72.0313,
}

_dspfvs_table_12 = {
    2: 46.311,
    3: 36.530,
    4: 47.870,
    6: 50.229,
    8: 53.289,
    12: 69.551,
    16: 71.600,
    24: 70.184,
    32: 72.138,
    48: 70.528,
    64: 72.348,
    96: 70.700,
    128: 72.524,
    192: 71.3333,
    256: 72.2500,
    384: 71.6667,
    512: 72.1250,
    1024: 72.0625,
    1536: 71.9167,
    2048: 72.0313,
}

_dspfvs_table_13 = {
    2: 2.750,
    3: 2.833,
    4: 2.875,
    6: 2.917,
    8: 2.938,
    12: 2.958,
    16: 2.969,
    24: 2.979,
    32: 2.984,
    48: 2.989,
    64: 2.992,
    96: 2.995,
}


def find_group_delay(attrs_dict):
    """
    Determine group delay from tables

    Args:
        attrs_dict (dict): dictionary of topspin acquisition parameters

    Returns:
        float: Group delay. Number of points FID is shifted by DSP. The ceiling of this number (group delay rounded up) is the number of points should be removed from the start of the FID.
    """

    group_delay = 0
    if attrs_dict["DSPFIRM"] != 0 and "GRPDLY" in attrs_dict.keys():
        group_delay = attrs_dict["GRPDLY"]
    elif attrs_dict["DECIM"] == 1.0:
        pass
    else:
        if attrs_dict["DSPFVS"] == 10:
            group_delay = _dspfvs_table_10[int(attrs_dict["DECIM"])]
        elif attrs_dict["DSPFVS"] == 11:
            group_delay = _dspfvs_table_11[int(attrs_dict["DECIM"])]
        elif attrs_dict["DSPFVS"] == 12:
            group_delay = _dspfvs_table_12[int(attrs_dict["DECIM"])]
        elif attrs_dict["DSPFVS"] == 13:
            group_delay = _dspfvs_table_13[int(attrs_dict["DECIM"])]
        else:
            print(
                "GRPDLY and DSPFVS parameters not found in acqus file, setting group delay to 0"
            )

    return group_delay


def import_topspin(path, phase_cycle=[0, 90, 180, 270]):
    """
    Import topspin data and return dnpdata object

    Args:
        path (str): Directory of data
        phase_cycle (list): list of phases used for phase cycling (deg, multiples of 90)

    Returns:
        dnpdata: topspin data
    """
    dir_list = _os.listdir(path)

    if "fid" in dir_list and "ser" not in dir_list:
        data = load_fid_ser(path, type="fid", phase_cycle=phase_cycle)
    elif "ser" in dir_list:
        data = load_fid_ser(path, type="ser", phase_cycle=phase_cycle)
    else:
        raise ValueError("Could Not Identify Data Type in File")

    return data


def load_acqu_proc(path="1", param_filename="acqus", proc_num=1):
    """
    Import topspin acqu or proc files

    Args:
        path (str): directory of acqu or proc file
        param_filename (str): parameters filename
        proc_num (int): number of the folder inside the pdata folder

    Returns:
        dict: Dictionary of acqusition parameters
    """

    if "proc" in param_filename:
        path_filename = _os.path.join(path, "pdata", str(procNum), param_filename)
    else:
        path_filename = _os.path.join(path, param_filename)

    # Import parameters
    with open(path_filename, "r") as f:
        raw_params = f.read()

    # Split parameters by line
    lines = raw_params.strip("\n").split("\n")
    attrs_dict = {}

    # Parse Parameters
    for line in lines:
        if line[0:3] == "##$":
            line_split = line[3:].split("= ")
            try:
                attrs_dict[line_split[0]] = float(line_split[1])
            except:
                attrs_dict[line_split[0]] = line_split[1]
            # if line_split[0] in ["TD", "NS"]:
            # print(line_split[0] + ": " + str(attrs_dict[line_split[0]]))

    needed_params = [
        "SW_h",
        "RG",
        "DECIM",
        "DSPFIRM",
        "DSPFVS",
        "BYTORDA",
        "TD",
        "SFO1",
    ]
    needed_params_2 = ["SW_h", "TD", "SFO1"]
    if param_filename in ["acqu", "acqus"]:
        if not all(
            map(
                attrs_dict.keys().__contains__,
                needed_params,
            )
        ):
            raise KeyError(
                "Unable to find all needed fields in the " + param_filename + " file"
            )
        else:
            attrs_dict = {x: attrs_dict[x] for x in needed_params}

    elif param_filename in ["acqu2", "acqu2s"]:
        if not all(map(attrs_dict.keys().__contains__, needed_params_2)):
            raise KeyError(
                "Unable to find all needed fields in the " + param_filename + " file"
            )
        else:
            attrs_dict = {x + "_2": attrs_dict[x] for x in needed_params_2}

    elif param_filename in ["acqu3", "acqu3s"]:
        if not all(map(attrs_dict.keys().__contains__, needed_params_2)):
            raise KeyError(
                "Unable to find all needed fields in the " + param_filename + " file"
            )
        else:
            attrs_dict = {x + "_3": attrs_dict[x] for x in needed_params_2}

    return attrs_dict


def load_fid_ser(path, type="fid", phase_cycle=None):
    """
    Import topspin fid or ser file

    Args:
        path (str): Directory of data
        type (str): "fid" for 1D, "ser" or "serPhaseCycle" for 2D
        phase_cycle (list): list of phases used for phase cycling (deg, multiples of 90)

    Returns:
        dnpdata: Topspin data
    """
    if isinstance(phase_cycle, list):
        for indx, x in enumerate(phase_cycle):
            if x in [0, 360, 720, 1080]:
                phase_cycle[indx] = 1
            elif x in [90, 450, 810, 1170]:
                phase_cycle[indx] = 1j
            elif x in [180, 540, 900, 1260]:
                phase_cycle[indx] = -1
            elif x in [270, 630, 990, 1350]:
                phase_cycle[indx] = -1j

    dir_list = _os.listdir(path)

    attrs_dict_list = [
        load_acqu_proc(path, x) for x in ["acqus", "acqu2s", "acqu3s"] if x in dir_list
    ]
    attrs_dict = {}
    for a in attrs_dict_list:
        attrs_dict.update(a)

    important_params_dict = {
        "nmr_frequency": attrs_dict["SFO1"] * 1e6,
        "SW_h": attrs_dict["SW_h"],
        "TD": attrs_dict["TD"],
    }

    higher_dim_pars = {
        x: attrs_dict[x]
        for x in ["SW_h_2", "TD_2", "SFO1_2", "SW_h_3", "TD_3", "SFO1_3"]
        if x in attrs_dict.keys()
    }
    important_params_dict.update(higher_dim_pars)

    if attrs_dict["BYTORDA"] == 0:
        endian = "<"
    else:
        endian = ">"

    if type == "fid":
        raw = _np.fromfile(_os.path.join(path, "fid"), dtype=endian + "i4")
    else:
        raw = _np.fromfile(_os.path.join(path, "ser"), dtype=endian + "i4")

    data = raw[0::2] + 1j * raw[1::2]  # convert to complex

    group_delay = find_group_delay(attrs_dict)
    group_delay = int(_np.ceil(group_delay))

    t = (
        1.0
        / attrs_dict["SW_h"]
        * _np.arange(0, int(attrs_dict["TD"] / 2) - group_delay)
    )

    if "vdlist" in dir_list:
        important_params_dict.update({"VDLIST": topspin_vdlist(path)})

    if type == "fid":
        data = data[group_delay : int(attrs_dict["TD"] / 2)] / attrs_dict["RG"]
        coords = [t]
        dims = ["t2"]
    elif type == "ser":
        ser_data = data.reshape(int(attrs_dict["TD_2"]), -1).T
        ser_data = (
            ser_data[group_delay : int(attrs_dict["TD"] / 2), :] / attrs_dict["RG"]
        )
        coords = [t, range(int(attrs_dict["TD_2"]))]
        dims = ["t2", "t1"]
        if "acqu2s" in dir_list and "acqu3s" not in dir_list:
            if "VDLIST" in important_params_dict.keys():
                if len(important_params_dict["VDLIST"]) == int(attrs_dict["TD_2"]):
                    coords = [t, important_params_dict["VDLIST"]]
            else:
                if isinstance(phase_cycle, list):
                    length1d = int((_np.ceil(attrs_dict["TD"] / 256.0) * 256) / 2)
                    ser_data = data.reshape(-1, int(length1d)).T
                    ser_data = ser_data[group_delay : int(attrs_dict["TD"] / 2), :]
                    phs_facs = _np.tile(
                        _np.array(phase_cycle),
                        int(attrs_dict["TD_2"]) / len(phase_cycle),
                    )
                    for indx in range(int(attrs_dict["TD_2"])):
                        ser_data[:, indx] = phs_facs[indx] * ser_data[:, indx]
                    ser_data = ser_data.sum(axis=1) / attrs_dict["RG"]
                    coords = [t]
                    dims = ["t2"]

        elif "acqu2s" in dir_list and "acqu3s" in dir_list:
            t = (
                1.0
                / attrs_dict["SW_h"]
                * _np.arange(
                    0, int(attrs_dict["TD"] / 2 / int(attrs_dict["TD_3"])) - group_delay
                )
            )
            if "VDLIST" in important_params_dict.keys() and len(
                important_params_dict["VDLIST"]
            ) == int(attrs_dict["TD_2"]):
                coords = [
                    t,
                    important_params_dict["VDLIST"],
                    range(int(attrs_dict["TD_3"])),
                ]
            else:
                coords = [
                    t,
                    range(int(attrs_dict["TD_2"])),
                    range(int(attrs_dict["TD_3"])),
                ]
            dims = ["t2", "t1", "t0"]
            ser_data = ser_data.reshape(
                int(attrs_dict["TD"] / 2 / int(attrs_dict["TD_3"])),
                int(attrs_dict["TD_2"]),
                int(attrs_dict["TD_3"]),
            )
        data = ser_data

    return _dnpdata(data, coords, dims, important_params_dict)


def topspin_vdlist(path):
    """
    Return topspin vdlist

    Args:
        path (str): Directory of data

    Returns:
        numpy.ndarray: vdlist as numpy array
    """
    fullPath = _os.path.join(path, "vdlist")

    with open(fullPath, "r") as f:
        raw = f.read()

    lines = raw.rstrip().rsplit()

    unitDict = {
        "n": 1.0e-9,
        "u": 1.0e-6,
        "m": 1.0e-3,
        "k": 1.0e3,
    }
    vdlist = []
    for line in lines:
        if line[-1] in unitDict:
            value = float(line[0:-1]) * unitDict[line[-1]]
            vdlist.append(value)
        else:
            value = float(line)
            vdlist.append(value)

    vdlist = _np.array(vdlist)
    return vdlist


def load_title(
    path="1", title_path=_os.path.join("pdata", "1"), title_filename="title"
):
    """
    Import Topspin Experiment Title File

    Args:
        path (str): Directory of title
        title_path (str): Path within experiment of title
        title_filename (str): filename of title

    Returns:
        str: Contents of experiment title file
    """

    path_filename = _os.path.join(path, title_path, title_filename)

    with open(path_filename, "r") as f:
        rawTitle = f.read()
    title = rawTitle.rstrip()

    return title


def topspin_jcamp_dx(path):
    """
    Return the contents of topspin JCAMP-DX file as dictionary

    Args:
        path: Path to file

    Returns:
        dict: Dictionary of JCAMP-DX file
    """

    attrs = {}

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip()

            if line[0:3] == "##$":
                key, value = tuple(line[3:].split("= ", 1))

                # Test for array
                if value[0] == "(":
                    x = _re.findall("\([0-9]+\.\.[0-9]+\)", value)

                    start, end = tuple(x[0][1:-1].split("..", 1))

                    array_size = int(end) + 1

                    same_line_array = value.split(")", 1)[-1]

                    array = []
                    if same_line_array != "":
                        same_line_array = same_line_array.split(" ")
                        same_line_array = [
                            float(x) if "." in x else int(x) for x in same_line_array
                        ]

                        array += same_line_array

                    while len(array) < array_size:
                        array_line = f.readline().rstrip().split(" ")

                        array_line = [
                            float(x) if "." in x else int(x) for x in array_line
                        ]

                        array += array_line

                    array = _np.array(array)

                    attrs[key] = array

                elif value[0] == "<":
                    value = value[1:-1]
                    if "." in value:
                        try:
                            value = float(value)
                        except:
                            pass
                    else:
                        try:
                            value = int(value)
                        except:
                            pass

                    attrs[key] = value
                else:
                    if "." in value:
                        try:
                            value = float(value)
                        except:
                            pass
                    else:
                        try:
                            value = int(value)
                        except:
                            pass

                    attrs[key] = value

            elif line[0:2] == "##":
                try:
                    key, value = tuple(line[2:].split("= ", 1))
                except:
                    pass

    return attrs
