import numpy as _np
import re
from warnings import warn

from .. import DNPData

import os

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

_required_params = {
    "acqus": [
        "SW_h",
        "RG",
        "DECIM",
        "DSPFIRM",
        "DSPFVS",
        "BYTORDA",
        "TD",
        "SFO1",
    ],
    "acqu2s": ["TD", "SW_h"],
    "acqu3s": ["TD", "SW_h"],
}


def find_group_delay(attrs_dict):
    """Determine group delay from tables

    Args:
        attrs_dict (dict): dictionary of topspin acquisition parameters

    Returns:
        float: Group delay. Number of points FID is shifted by DSP. The ceiling of this number (group delay rounded up) is the number of points should be removed from the start of the FID.
    """

    # This must be revisited
    group_delay = 0
    if attrs_dict["DSPFVS"] >= 13 and "GRPDLY" in attrs_dict.keys():
        group_delay = attrs_dict["GRPDLY"]
        if group_delay > 0:
            return group_delay

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


# This function does too much, should be broken into smaller functions
def import_topspin(
    path, assign_vdlist=False, remove_digital_filter=False, verbose=False
):
    """Import topspin data and return dnpdata object

    Args:
        path (str): Directory of data
        assign_vdlist: False, or the name of dimension to assign topspin vdlist
        remove_digital_filter (bool): option to remove group delay
        verbose (bool): Print additional output for troubleshooting

    Returns:
        dnpdata: topspin data
    """
    dir_list = os.listdir(path)  # All files and folders in directory
    if verbose:
        print("Files in directory:")
        for each in dir_list:
            print(" ", each)

    # Load Acquisition Parameters
    if verbose:
        print("Loading acqus")

    acqus_params = load_acqu(os.path.join(path, "acqus"), verbose=verbose)

    dims = [
        "t2"
    ]  # this may cause issues if the first dimension is not the direct time dimension

    if "fid" in dir_list:
        bin_filename = "fid"
    else:
        bin_filename = "ser"

    if verbose:
        print("Binary File:", bin_filename)

    if acqus_params["BYTORDA"] == 0:
        endian = "<"
    else:
        endian = ">"

    if verbose:
        print("Endian", endian)

    topspin_major_version = int(acqus_params["topspin"].split(".")[0])

    # Is this incorrect?
    # Most topspin data I've seen is i4, however, later versions seem to have i8
    # float is also possible
    if topspin_major_version >= 4:
        data_bytes = 8
    else:
        data_bytes = 4

    if acqus_params["DTYPA"] == 0:
        data_type = "i"
    elif acqus_params["DTYPA"] == 2:
        data_type = "f"

    data_type = data_type + str(data_bytes)

    if verbose:
        print("Data type:", endian + data_type)

    raw = load_bin(os.path.join(path, bin_filename), dtype=endian + data_type)

    if verbose:
        print("Raw data shape:", _np.shape(raw))

    # Convert to complex data depending on AQ_mod
    if acqus_params["AQ_mod"] == 0 or acqus_params["AQ_mod"] == 2:
        if verbose:
            print("Data is complex")

        values = raw

    elif acqus_params["AQ_mod"] == 1 or acqus_params["AQ_mod"] == 3:
        if verbose:
            print("Data is not complex")

        values = raw[0::2] + 1j * raw[1::2]  # convert to complex

    else:
        raise ValueError("Unknown format")

    # Option to remove group delay
    if remove_digital_filter == True:
        group_delay = find_group_delay(acqus_params)
        group_delay = int(_np.floor(group_delay))  # should this be floor or ceil?

    else:
        group_delay = 0

    if verbose:
        print("Group Delay", group_delay)

    # why is dividing by 2 required?
    t2 = 1.0 / acqus_params["SW_h"] * _np.arange(0, int(acqus_params["TD"] / 2))

    if verbose:
        print("Points in FID:", acqus_params["TD"] / 2)

    # Handle t2 group delay
    # t2 = t2[slice(group_delay, int(acqus_params["TD"] / 2))] # Alternative method
    t2 = t2[group_delay:]
    coords = [t2]

    # This will not work for vdlist data
    if "acqu2s" in dir_list:
        if verbose:
            print("Loading acqu2s")

        acqu2s_params = load_acqu(os.path.join(path, "acqu2s"), verbose=verbose)
        dims.insert(0, "t1")
        t1 = 1.0 / acqu2s_params["SW_h"] * _np.arange(0, int(acqu2s_params["TD"]))
        coords.insert(0, t1)

    # 3d data must be tested
    if "acqu3s" in dir_list:
        if verbose:
            print("Loading acqu3s")

        acqu3s_params = load_acqu(os.path.join(path, "acqu3s"), verbose=verbose)
        dims.insert(0, "t3")
        t3 = 1.0 / acqu3s_params["SW_h"] * _np.arange(0, int(acqu3s_params["TD"]))
        coords.insert(0, t3)

    new_shape = [len(coords[ix]) if dims[ix] != "t2" else -1 for ix in range(len(dims))]
    if verbose:
        print("Raw Data Shape:", _np.shape(values))
        print("Reshaping data to:", new_shape)

    if assign_vdlist:
        if verbose:
            print("Assigning vdlist to %s dim" % assign_vdlist)

        vdlist = topspin_vdlist(path)
        if assign_vdlist in dims:
            index = dims.index(assign_vdlist)
            coords[index] = vdlist

        else:
            raise ValueError("Could not identify dimension to assign vdlist")

    # reshape values
    values = values.reshape(new_shape)

    # Handle group delay
    values = values[..., slice(group_delay, int(acqus_params["TD"] / 2))]

    # create data object
    topspin_data = DNPData(values, dims, coords, attrs=acqus_params)

    # Add import path to attributes
    topspin_data.attrs["import_path"] = path

    # Add NMR Frequency to attrs
    topspin_data.attrs["nmr_frequency"] = acqus_params["SFO1"] * 1e6

    # Assign data/spectrum type
    topspin_data.attrs["experiment_type"] = "nmr_spectrum"

    # Assign number of scans
    if "acqu2s" in dir_list:
        topspin_data.attrs["scans"] = int(acqu2s_params["TD"])

    # reorder so that 't2' is first
    topspin_data.reorder(["t2"])

    return topspin_data


# Load topspin should also handle this
def load_pdata(path, verbose=False):
    """Import prospa processed data

    Args:
        path (str): Directory of pdata
        verbose (bool): If true, print output for troubleshooting

    Returns:
        DNPData: Topspin processed data
    """

    # Directory
    dir_list = os.listdir(path)  # All files and folders in directory
    if verbose:
        print("Files in directory:")
        for each in dir_list:
            print(" ", each)

    proc_params = load_acqu(os.path.join(path, "procs"), verbose=verbose)

    if proc_params["BYTORDP"] == 0:
        endian = "<"
    else:
        endian = ">"

    if verbose:
        print("endian", endian)

    real_raw = load_bin(os.path.join(path, "1r"), dtype=endian + "i4")
    imag_raw = load_bin(os.path.join(path, "1i"), dtype=endian + "i4")

    SW = proc_params["SW_p"]
    offset = proc_params["OFFSET"]  # Reference Offset in ppm
    td_eff = proc_params["TDeff"]
    SI = proc_params["SI"]  # What does SI stand for?
    spectrometer_frequency = proc_params["SF"]  # spectrometer frequency in MHz
    phase_0 = proc_params["PHC0"]  # Phase correction, zeroth order phase
    phase_1 = proc_params["PHC1"]  # Phase correction, first order phase

    f2 = (
        -1 * SW * _np.linspace(0, 1, num=SI, endpoint=False) / spectrometer_frequency
        + offset
    )

    raw = real_raw + 1j * imag_raw

    data = DNPData(raw, ["f2"], [f2], attrs=proc_params)
    data.attrs["nmr_frequency"] = spectrometer_frequency
    data.attrs["phase_0"] = phase_0
    data.attrs["phase_1"] = phase_1

    return data


# Is this function obsolete?
def load_acqu(path, required_params=None, verbose=False):
    """Import topspin acqu or proc files

    Args:
        path (str): directory of acqu or proc file
        required_params (list): Only return parameters given
        verbose (bool): If true, print output for troubleshooting

    Returns:
        dict: Dictionary of acqusition parameters
    """

    raw_params = load_topspin_jcamp_dx(path, verbose=False)

    if required_params is not None:
        acqus_params = {}
        for key in required_params:
            acqus_params[key] = raw_params[key]
    else:
        acqus_params = raw_params

    return acqus_params


def topspin_vdlist(path):
    """
    Return topspin vdlist

    Args:
        path (str): Directory of data

    Returns:
        numpy.ndarray: vdlist as numpy array
    """
    fullPath = os.path.join(path, "vdlist")

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


# Legacy import ser
def load_ser(path, dtype=">i4"):
    """Depreciated. Use load bin. Import Topspin Ser file

    Args:
        path (str): Directory of data
        dtype (str): data format for import

    returns:
        raw (np.ndarray): Data from ser file
    """
    warn(
        "This function is deprecated, use load_bin instead",
        DeprecationWarning,
        stacklevel=2,
    )

    raw = _np.fromfile(os.path.join(path), dtype=dtype)

    return raw


def load_bin(path, dtype=">i4"):
    """Import Topspin Ser file

    Args:
        path (str): Directory of data
        dtype (str): data format for import

    returns:
        raw (np.ndarray): Data from ser file
    """

    raw = _np.fromfile(os.path.join(path), dtype=dtype)

    return raw


def load_title(path="1", title_path=os.path.join("pdata", "1"), title_filename="title"):
    """Import Topspin Experiment Title File

    Args:
        path (str): Directory of title
        title_path (str): Path within experiment of title
        title_filename (str): filename of title

    Returns:
        str: Contents of experiment title file
    """

    path_filename = os.path.join(path, title_path, title_filename)

    with open(path_filename, "r") as f:
        rawTitle = f.read()
    title = rawTitle.rstrip()

    return title


def load_topspin_jcamp_dx(path, verbose=False):
    """Return the contents of topspin JCAMP-DX file as dictionary

    Args:
        path (str): Path to file
        verbose (bool): If true, print output for troubleshooting

    Returns:
        dict: Dictionary of JCAMP-DX file parameters
    """

    attrs = {}

    with open(path, "r") as f:
        for line in f:
            if verbose:
                print(line)
            line = line.rstrip()

            if line[0:3] == "##$":
                key, value = tuple(line[3:].split("= ", 1))

                # Test for array
                if value[0] == "(":
                    x = re.findall("\([0-9]+\.\.[0-9]+\)", value)

                    start, end = tuple(x[0][1:-1].split("..", 1))

                    array_size = int(end) + 1

                    same_line_array = value.split(")", 1)[-1]

                    array = []
                    if same_line_array != "":
                        same_line_array = same_line_array.split(" ")

                        try:
                            same_line_array = [
                                float(x) if "." in x else int(x)
                                for x in same_line_array
                            ]
                        except:
                            pass  # Needed in case where "<>" found in some arrays

                        array += same_line_array

                    while len(array) < array_size:
                        array_line = f.readline().rstrip().split(" ")

                        try:
                            array_line = [
                                float(x) if "." in x else int(x) for x in array_line
                            ]
                        except:
                            pass  # Needed in case where "<>" found in some arrays

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
                # Extract Title and TopSpin Version, needed for data type determination
                #                if 'TOPSPIN' in line.upper():
                if "TITLE" in line:
                    version = line.split(" ")[-1]
                    attrs["topspin"] = version
                try:
                    key, value = tuple(line[2:].split("= ", 1))
                except:
                    pass

    return attrs
