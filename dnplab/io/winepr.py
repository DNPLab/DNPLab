import numpy as _np
import os
from .. import DNPData
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
    """Import Bruker par/spc data and return DNPData object

    Args:
        path (str) : Path to either .par or .spc file

    Returns:
        parspc_data (object) : DNPData object containing Bruker par/spc data
    """

    pathexten = os.path.splitext(path)[1]
    path = os.path.splitext(path)[0]
    if pathexten == ".par" or pathexten == ".spc":
        path_par = path + ".par"
        path_spc = path + ".spc"

    else:
        raise TypeError("data file must be .spc or .par")

    attrs = load_par(path_par)
    values, dims, coords, attrs = load_spc(path_spc, attrs)

    # Assign data/spectrum type
    attrs["experiment_type"] = "epr_spectrum"

    parspc_data = DNPData(values, dims, coords, attrs)

    return parspc_data


def load_par(path):
    """Import contents of .par file

    Args:
        path (str) : Path to .par file

    Returns:
        attrs (dict) : dictionary of parameters
    """

    attrs = {}
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip()
            split_line = line.split(" ", 1)

            if len(split_line) == 2:
                key = split_line[0].strip()
                value = split_line[1].strip()
                attrs[key] = value

    for key in rename_dict:
        if key in attrs:
            new_key = rename_dict[key]
            if new_key in float_params:
                attrs[new_key] = float(attrs[key])
            elif new_key in int_params:
                attrs[new_key] = int(float(attrs[key]))
            else:
                attrs[new_key] = attrs[key]

    if "DOS" in attrs:
        attrs["endian"] = "LIT"
        attrs["data_type"] = "float32"
    else:
        attrs["endian"] = "BIG"
        attrs["data_type"] = "int32"

    return attrs


def load_spc(path, attrs):
    """Import data and axes of .spc file

    Args:
        path (str) : Path to .spc file

    Returns:
        coords (ndarray) : coordinates for spectrum or spectra
        values (ndarray) : data values
        attrs (dict) : updated dictionary of parameters
        dims (list) : dimension labels
    """

    data_format = _np.dtype(attrs["data_type"]).newbyteorder(attrs["endian"])
    file_opened = open(path, "rb")
    file_bytes = file_opened.read()
    values = _np.frombuffer(file_bytes, dtype=data_format)

    attrs.pop("data_type", None)
    attrs.pop("endian", None)

    if "x_points" not in attrs.keys():
        if "y_points" not in attrs.keys():
            attrs["x_points"] = int(len(values))
        else:
            attrs["x_points"] = int(len(values) / attrs["y_points"])

    if "center_field" not in attrs.keys() or "x_width" not in attrs.keys():
        if "sweep_start" in attrs.keys() and "sweep_extent" in attrs.keys():
            coords = [
                _np.linspace(
                    attrs["sweep_start"],
                    attrs["sweep_start"] + attrs["sweep_extent"],
                    attrs["x_points"],
                )
            ]
        else:
            warnings.warn("not axis information, axis is indexed only")
            coords = [range(attrs["x_points"])]
    elif "center_field" in attrs.keys() and "x_width" in attrs.keys():
        coords = [
            _np.linspace(
                attrs["center_field"] - attrs["x_width"] / 2,
                attrs["center_field"] + attrs["x_width"] / 2,
                attrs["x_points"],
            )
        ]
    else:
        warnings.warn("unable to define axis, indexed only")
        coords = [range(attrs["x_points"])]

    if "x_unit" in attrs.keys() and attrs["x_unit"] in ["G", "T"]:
        if attrs["x_unit"] == "G":
            coords = [x / 10 for x in coords]
        elif attrs["x_unit"] == "T":
            coords = [x * 1000 for x in coords]
        dims = ["B0"]
    else:
        dims = ["t2"]

    if "y_points" in attrs.keys() and attrs["y_points"] != 1:
        values = _np.reshape(values, (attrs["x_points"], attrs["y_points"]), order="F")
        dims.append("t1")

        if "y_min" in attrs.keys() and "y_width" in attrs.keys():
            coords.append(
                _np.linspace(
                    attrs["y_min"],
                    attrs["y_min"] + attrs["y_width"],
                    attrs["y_points"],
                )
            )
        else:
            warnings.warn("unable to define indirect axis")

    file_opened.close()

    return values, dims, coords, attrs
