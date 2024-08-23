import numpy as _np
import os
import dnplab as _dnp
import re

scale_dict = {
    "p": 1e-12,
    "n": 1e-9,
    "u": 1e-6,
    "m": 1e-3,
    "1": 1,
    "k": 1e3,
    "M": 1e6,
    "G": 1e9,
    "T": 1e12,
}


def import_specman(
    path, autodetect_coords: bool = False, autodetect_dims: bool = False
):
    """Import SpecMan data and return DNPData object

    DNPLab function to import SpecMan4EPR data (https://specman4epr.com/). The function returns a DNPdata object with the spectral data.

    The structure of the DNPdata object can be complex and the variables saved by SpecMan depend on the individual spectrometer configuration. Therefore, the import function returns a numpy array with the dimension "x0", "x1", "x2", "x3", "x4". In any case, the dimension "x0" corresponds to the variables stored in the data file. The spectroscopic data is stored in "x1" to "x4", depending on how many dimensions were recorded.
    The import function will require a parser script to properly assign the spectroscopic data and proper coordinates.

    Args:
        path (str):                 Path to either .exp file
        autodetect_coords(bool):    Autodetect coords based on attrs
        autodetect_dims(bool):      Autodetect dims based on attrs
    Returns:
        data (DNPData):         DNPData object containing SpecMan EPR data
    """

    if path[-1] == os.sep:
        path = path[:-1]
    if path[-1] == "p":
        file_name_d01 = path.replace("exp", "d01")
        file_name_exp = path
    elif path[-1] == "1":
        file_name_exp = path.replace("d01", "exp")
        file_name_d01 = path
    else:
        raise TypeError("Incorrect file type, must be .d01 or .exp")

    attrs = load_specman_exp(file_name_exp)
    data, dims, coords, attrs = load_specman_d01(file_name_d01, attrs)

    if autodetect_coords or autodetect_dims:
        attrs = analyze_attrs(attrs)

    if autodetect_dims:
        new_dims = generate_dims(attrs)
        dims = new_dims
    else:
        new_dims = None

    if autodetect_coords:
        coords = calculate_specman_coords(attrs, new_dims)

    # Add import path
    attrs["import_path"] = path

    # Assign data/spectrum type
    attrs["experiment_type"] = "epr_spectrum"

    specman_data = _dnp.DNPData(data, dims, coords, attrs)

    return specman_data


def load_specman_exp(path):
    """Import SpecMan parameters

    DNPLab function to read and import the SpecMan exp file. The .exp file is a text file that stores the experimental data, the pulse program, and other spectrometer configuration files.

    Args:
        path (str):     Path to either .d01 or .exp file

    Returns:
        attrs (dict):   Dictionary of parameter fields and values (DNPLab attributes)

    """
    exp_file_opened = open(path, encoding="utf8", errors="ignore")
    file_contents = exp_file_opened.read().splitlines()
    exp_file_opened.close()

    attrs = {}

    c = ""

    for i in range(0, len(file_contents)):
        exp_content = str(file_contents[i])
        splt_exp_content = exp_content.split(" = ")
        if "[" in exp_content and "]" in exp_content and "=" not in exp_content:
            c = splt_exp_content[0].replace("[", "").replace("]", "")
        elif exp_content == "":
            c = "param"
        elif len(splt_exp_content) > 1:
            attrs[c + "_" + splt_exp_content[0]] = splt_exp_content[1]
        elif len(splt_exp_content) == 1 and exp_content != "":
            attrs[c + "_" + str(i)] = splt_exp_content
        else:
            pass

    return attrs


def load_specman_d01(path, attrs, verbose=False):
    """Import SpecMan d01 data file

    DNPLab function to import the SpecMan d01 data file. The format of the SpecMan data file is described here:


    Args:
        path (str):         Path to either .d01 or .exp file

    Returns:
        data (ndarray):     SpecMan data as numpy array
        params (dict):      Dictionary with import updated parameters dictionary
    """

    file_opened = open(path, "rb")
    uint_read = _np.fromfile(file_opened, dtype=_np.uint32)
    file_opened.close()

    file_opened = open(path, "rb")
    # "<f4" means float32, little-endian
    float_read = _np.fromfile(file_opened, dtype="<f4")
    file_opened.close()

    # Number of recorded variables stored
    attrs["numberOfVariables"] = uint_read[0]

    # 0 - data stored in double format, 1 -data stored in float format
    attrs["dataFormat"] = uint_read[1]

    # Number of dimensions of stored data. Note: Not clear why this is repeated for every dimension. This seems to be the same number for all dimensions, but is repeated.
    attrs["dims"] = uint_read[2]

    dataShape = uint_read[2 : 2 + uint_read[0] * 6]

    attrs["dataStreamShape"] = dataShape

    attrs["dataStartIndex"] = 2 + attrs["numberOfVariables"] * 6

    if verbose == True:
        print("** Data paramters **")
        print("numberOfVariables : ", attrs["numberOfVariables"])
        print("dataFormat        : ", attrs["dataFormat"])
        print("dims              : ", attrs["dims"])
        print("dataStreamShape   : ", attrs["dataStreamShape"])
        print("dataStartIndex    : ", attrs["dataStartIndex"])

    data = float_read[attrs["dataStartIndex"] :]

    if attrs["dims"] == 1:
        data = _np.reshape(data, (uint_read[0], uint_read[3]), order="C")

    elif attrs["dims"] == 2:
        data = _np.reshape(data, (uint_read[0], uint_read[4], uint_read[3]), order="C")

    elif attrs["dims"] == 3:
        data = _np.reshape(
            data, (uint_read[0], uint_read[5], uint_read[4], uint_read[3]), order="C"
        )

    elif attrs["dims"] == 4:
        data = _np.reshape(
            data,
            (uint_read[0], uint_read[4], uint_read[5], uint_read[6], uint_read[3]),
            order="C",
        )

    elif attrs["dims"] >= 4:
        print("Maximum dimensionality for SpecMan data is 4D")
        return None

    # Swap first axis with last
    data = _np.swapaxes(data, 0, -1)

    dims_full = ["x0", "x1", "x2", "x3", "x4"]
    dims = dims_full[0 : dataShape[0] + 1]

    coords = []
    shape = _np.shape(data)
    for index in range(data.ndim):
        coords.append(_np.arange(0.0, shape[index]))

    # SpecMan data can have a maximum of four dimensions

    return data, dims, coords, attrs


def analyze_attrs(attrs):
    """
    Analyze the attrs and add some important attrs to existing dictionary

    Args:
        attrs (dict): Dictionary of specman acqusition parameters

    Returns:
        attrs (dict): The dictionary of specman acqusition parameters and added parameters

    """

    temp = {}
    for key, val in attrs.items():
        if "params_" in key:
            new_key = key.split("params_")[1]  # get key value for temp dictionary
            val = val.split(";")[0]  # remove non value related information
            val_list = val.split(" ")  # split value string for further analyze

            val = val_list[0].strip(",")
            val_unit = val_list[1] if len(val_list) == 5 else None
            temp[new_key] = int(val) if "." not in val else float(val)
            temp[new_key] *= _convert_unit(val_unit)

            if "step" in val_list:  # when it indicate the step
                step_index = (
                    val_list.index("step") + 1
                )  # the index of the value of 'step' is equal to the index of string 'index' + 1
                step_unit = (
                    val_list[val_list.index("step") + 2] if len(val_list) == 5 else None
                )
                step = float(val_list[step_index]) * _convert_unit(step_unit)
                temp[new_key + "_step"] = step

            if "to" in val_list:  # when it indicate the stop
                stop_index = (
                    val_list.index("to") + 1
                )  # the index of the value of 'stop' is equal to the index of string 'index' + 1
                stop_unit = (
                    val_list[val_list.index("to") + 2] if len(val_list) == 5 else None
                )
                stop = float(val_list[stop_index]) * _convert_unit(stop_unit)
                temp[new_key + "_stop"] = stop

        if "sweep_" in key:
            val_list = val.split(",")
            val = val_list[1]  # get value
            new_key = "sweep_" + val_list[0]
            temp[new_key + "_length"] = int(val)
            # new_key += '_dim' # last item is the key to the parameters, such as t, p...
            temp[new_key + "_dim"] = val_list[3]

    attrs = {**attrs, **temp}
    return attrs


def generate_dims(attrs):
    """Generate dims from specman acquisition parameters

    Args:
        attrs (dict): Dictionary of specman acqusition parameters

    Returns:
        dims (list): a new dims

    """
    kw = ["sweep_T", "sweep_X", "sweep_Y", "sweep_Z"]
    dims = [
        attrs[key + "_dim"] if key != "sweep_T" else "t2"
        for key in kw
        if key + "_dim" in attrs
    ]
    dims.append("x")

    return dims


def calculate_specman_coords(attrs, dims=None):
    """Generate coords from specman acquisition parameters

    Args:
        attrs (dict): Dictionary of specman acqusition parameters
        dims (list): (Optional) a list of dims

    Returns:
        coords (list): a calculated coords
    """

    kw = ["sweep_T", "sweep_X", "sweep_Y", "sweep_Z"]
    coords = []
    lengths = [attrs[key + "_length"] for key in kw if key + "_length" in attrs]
    lengths.append(2)

    if not dims:
        dims = generate_dims(attrs)
        print("Warning: the coords might not be correct")

    for index, dim in enumerate(dims):
        length = lengths[index]
        if dim in attrs and dim + "_step" in attrs:
            start = attrs[dim]
            step = attrs[dim + "_step"]
            stop = start + step * (length - 1)
            coord = _np.linspace(start, stop, length)
        elif dim in attrs and dim + "_stop" in attrs:
            start = attrs[dim]
            stop = attrs[dim + "_stop"]
            coord = _np.linspace(start, stop, length)
        elif dim in attrs and dim + "_step" not in attrs and dim + "_stop" not in attrs:
            val_string = attrs["params_" + dim].split(";")[0]
            coord = _np.array(
                [float(f) for f in val_string.split() if f.replace(".", "").isdigit()]
            )
        else:
            coord = _np.arange(0.0, length)

        coords.append(_np.array(coord))

    return coords


def _convert_unit(unit_string=None) -> float:
    if not unit_string:
        return 1.0

    if len(unit_string) != 1 and unit_string.lower() != "hz":
        if unit_string[0] in scale_dict:
            return scale_dict[unit_string[0]]

    return 1.0
