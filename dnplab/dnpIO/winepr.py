import numpy as np
import os
from .. import dnpdata
import warnings

__all__ = ["import_winepr", "load_par", "load_spc"]

rename_dict = {
    "MF": "frequency",
    "MP": "power",
    "MPD": "attenuation",
    "HCF": "center_field",
    "RCT": "conversion_time",
    "RTC": "time_constant",
    "GST": "sweep_start",
    "GSI": "sweep_extent",
    "RMA": "modulation_amplitude",
    "RRG": "receiver_gain",
    "JSD": "nscans",
    "TE": "temperature",
    "JUN": "x_unit",
    "SSX": "x_points",
    "HSW": "x_width",
    "XXUN": "x_unit",
    "XXLB": "x_min",
    "XYUN": "y_unit",
    "SSY": "y_points",
    "XYWI": "y_width",
    "XYLB": "y_min",
}

float_params = [
    "frequency",
    "power",
    "center_field",
    "conversion_time",
    "time_constant",
    "sweep_start",
    "sweep_extent",
    "modulation_amplitude",
    "receiver_gain",
    "temperature",
    "x_width",
    "x_min",
    "y_min",
    "y_width",
]

int_params = [
    "attenuation",
    "nscans",
    "x_points",
    "y_points",
]


def import_winepr(path):
    """
    Import Bruker par/spc data and return dnpdata object

    Args:
        path (str) : Path to either .par or .spc file

    Returns:
        parspc_data (object) : dnpdata object containing Bruker par/spc data
    """

    pathexten = os.path.splitext(path)[1]
    path = os.path.splitext(path)[0]
    if pathexten == ".par" or pathexten == ".spc":
        path_par = path + ".par"
        path_spc = path + ".spc"

    else:
        raise TypeError("data file must be .spc or .par")

    params = load_par(path_par)
    values, coords, dims, attrs = load_spc(path_spc, params)

    parspc_data = dnpdata(values, coords, dims, attrs)

    return parspc_data


def load_par(path):
    """
    Import contents of .par file

    Args:
        path (str) : Path to .par file

    Returns:
        params (dict) : dictionary of parameters
    """

    params = {}
    with open(path, "r") as f:

        for line in f:
            line = line.rstrip()
            split_line = line.split(" ", 1)

            if len(split_line) == 2:
                key = split_line[0].strip()
                value = split_line[1].strip()
                params[key] = value

    for key in rename_dict:
        if key in params:
            new_key = rename_dict[key]
            if new_key in float_params:
                params[new_key] = float(params[key])
            elif new_key in int_params:
                params[new_key] = int(float(params[key]))
            else:
                params[new_key] = params[key]

    if "DOS" in params:
        params["endian"] = "LIT"
        params["data_type"] = "float32"
    else:
        params["endian"] = "BIG"
        params["data_type"] = "int32"

    return params


def load_spc(path, params):
    """
    Import data and axes of .spc file

    Args:
        path (str) : Path to .spc file

    Returns:
        abscissa (ndarray) : coordinates for spectrum or spectra
        spec (ndarray) : data values
        params (dict) : updated dictionary of parameters
        dims (list) : dimension labels
    """

    data_format = np.dtype(params["data_type"]).newbyteorder(params["endian"])
    file_opened = open(path, "rb")
    file_bytes = file_opened.read()
    spec = np.frombuffer(file_bytes, dtype=data_format)

    params.pop("data_type", None)
    params.pop("endian", None)

    if "x_points" not in params.keys():
        if "y_points" not in params.keys():
            params["x_points"] = int(len(spec))
        else:
            params["x_points"] = int(len(spec) / params["y_points"])

    if "center_field" not in params.keys() or "x_width" not in params.keys():
        if "sweep_start" in params.keys() and "sweep_extent" in params.keys():
            abscissa = [
                np.linspace(
                    params["sweep_start"],
                    params["sweep_extent"],
                    params["x_points"],
                )
            ]
        else:
            warnings.warn("not axis information, axis is indexed only")
            abscissa = [range(params["x_points"])]
    elif "center_field" in params.keys() and "x_width" in params.keys():
        abscissa = [
            np.linspace(
                params["center_field"] - params["x_width"] / 2,
                params["center_field"] + params["x_width"] / 2,
                params["x_points"],
            )
        ]
    else:
        warnings.warn("unable to define axis, indexed only")
        abscissa = [range(params["x_points"])]

    if "x_unit" in params.keys() and params["x_unit"] in ["G", "T"]:
        if params["x_unit"] == "G":
            abscissa = [x / 10 for x in abscissa]
        elif params["x_unit"] == "T":
            abscissa = [x * 1000 for x in abscissa]
    dims = ["t2"]

    if "y_points" in params.keys() and params["y_points"] != 1:
        spec = np.reshape(spec, (params["x_points"], params["y_points"]), order="F")
        dims.append("t1")

        if "y_min" in params.keys() and "y_width" in params.keys():
            abscissa.append(
                np.linspace(
                    params["y_min"],
                    params["y_min"] + params["y_width"],
                    params["y_points"],
                )
            )
        else:
            warnings.warn("unable to define indirect axis")

    file_opened.close()

    return spec, abscissa, dims, params
