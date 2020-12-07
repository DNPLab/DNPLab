import numpy as np
import os
from dnplab import dnpdata
import warnings


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

    file_opened = open(path, "r")
    parfile_contents = file_opened.readlines()

    params = {}
    for ix in range(len(parfile_contents)):
        par = parfile_contents[ix].rstrip("\n")
        if "MF" in par and "EMF" not in par:
            params["frequency"] = float(par.replace("MF", "").strip())
        elif "MP" in par and "MPD" not in par:
            params["power"] = float(par.replace("MP", "").strip())
        elif "MPD" in par:
            params["attenuation"] = int(float(par.replace("MPD", "").strip()))
        elif "HCF" in par:
            params["center_field"] = float(par.replace("HCF", "").strip())

        elif "RCT" in par:
            params["conversion_time"] = float(par.replace("RCT", "").strip())
        elif "RTC" in par:
            params["time_constant"] = float(par.replace("RTC", "").strip())
        elif "GST" in par:
            params["sweep_start"] = float(par.replace("GST", "").strip())
        elif "GSI" in par:
            params["sweep_extent"] = float(par.replace("GSI", "").strip())
        elif "RMA" in par:
            params["modulation_amplitude"] = float(par.replace("RMA", "").strip())
        elif "RRG" in par:
            params["receiver_gain"] = float(par.replace("RRG", "").strip())
        elif "JSD" in par:
            params["nscans"] = int(par.replace("JSD", "").strip())
        elif "TE" in par and "TE1" not in par:
            params["temperature"] = float(par.replace("TE", "").strip())
        elif "JUN" in par and "JDA" not in par:
            params["x_unit"] = par.replace("JUN", "").strip()
        elif "SSX" in par:
            params["x_points"] = int(par.replace("SSX", "").strip())
        elif "HSW" in par:
            params["x_width"] = float(par.replace("HSW", "").strip())

        elif "XXUN" in par:
            params["x_unit"] = par.replace("XXUN", "").strip()
        elif "XXLB" in par:
            params["x_min"] = float(par.replace("XXLB", "").strip())

        elif "XYUN" in par:
            params["y_unit"] = par.replace("XYUN", "").strip()
        elif "SSY" in par:
            params["y_points"] = int(par.replace("SSY", "").strip())
        elif "XYLB" in par:
            params["y_min"] = float(par.replace("XYLB", "").strip())
        elif "XYWI" in par:
            params["y_width"] = float(par.replace("XYWI", "").strip())

        elif "DOS" in par:
            params["endian"] = "LIT"
            params["data_type"] = "float32"

    if "data_type" not in params.keys():
        params["endian"] = "BIG"
        params["data_type"] = "int32"

    file_opened.close()

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

    if "x_unit" not in params.keys():
        dims = ["t2"]
    else:
        dims = [params["x_unit"]]

    if "y_points" in params.keys() and params["y_points"] != 1:
        spec = np.reshape(spec, (params["x_points"], params["y_points"]), order="F")

        if "y_unit" not in params.keys():
            dims.append("t1")
        else:
            dims.append(params["y_unit"])
        if dims[0] == dims[1]:
            dims = ["t2", "t1"]

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
