import numpy as np
import os
from .. import dnpdata
import warnings


def import_bes3t(path):
    """
    Import Bruker BES3T data and return dnpdata object

    Args:
        path (str) : Path to either .DSC or .DTA file

    Returns:
        bes3t_data (object) : dnpdata object containing Bruker BES3T data
    """

    pathexten = os.path.splitext(path)[1]
    path = os.path.splitext(path)[0]
    path_xgf = "none"
    path_ygf = "none"
    path_zgf = "none"
    if pathexten == ".DSC" or pathexten == ".DTA":
        path_dsc = path + ".DSC"
        path_dta = path + ".DTA"
        if os.path.isfile(path + ".XGF"):
            path_xgf = path + ".XGF"
        if os.path.isfile(path + ".YGF"):
            path_ygf = path + ".YGF"
        if os.path.isfile(path + ".ZGF"):
            path_zgf = path + ".ZGF"

    elif pathexten == ".XGF":
        path_xgf = path + ".XGF"
        path_dsc = path + ".DSC"
        path_dta = path + ".DTA"
        if os.path.isfile(path + ".YGF"):
            path_ygf = path + ".YGF"
        if os.path.isfile(path + ".ZGF"):
            path_zgf = path + ".ZGF"

    elif pathexten == ".YGF":
        path_ygf = path + ".YGF"
        path_dsc = path + ".DSC"
        path_dta = path + ".DTA"
        if os.path.isfile(path + ".XGF"):
            path_xgf = path + ".XGF"
        if os.path.isfile(path + ".ZGF"):
            path_zgf = path + ".ZGF"

    elif pathexten == ".ZGF":
        path_zgf = path + ".ZGF"
        path_dsc = path + ".DSC"
        path_dta = path + ".DTA"
        if os.path.isfile(path + ".XGF"):
            path_xgf = path + ".XGF"
        if os.path.isfile(path + ".YGF"):
            path_ygf = path + ".YGF"

    else:
        raise TypeError("data file must be .DTA, .DSC, .YGF, or .ZGF")

    params = load_dsc(path_dsc)
    values, coords, dims, attrs = load_dta(
        path_dta, path_xgf, path_ygf, path_zgf, params
    )

    bes3t_data = dnpdata(values, coords, dims, attrs)

    return bes3t_data


def load_dsc(path):
    """
    Import contents of .DSC file

    Args:
        path (str) : Path to .DSC file

    Returns:
        params (dict) : dictionary of parameters
    """

    file_opened = open(path, "r")
    dscfile_contents = file_opened.readlines()
    file_opened.close()

    params = {}
    sweep_domain = None
    for ix in range(len(dscfile_contents)):
        try:
            par = dscfile_contents[ix].rstrip("\t").rstrip("\n")
            if "MWFQ" in par:
                params["frequency"] = float(par.replace("MWFQ", "").strip()) / 1e9
            elif "Power" in par and "Atten" not in par:
                params["power"] = float(
                    par.replace("Power", "").replace("mW", "").strip()
                )
            elif "PowerAtten" in par:
                params["attenuation"] = int(
                    float(par.replace("PowerAtten", "").replace("dB", "").strip())
                )
            elif "Attenuation" in par:
                params["pulse_attenuation"] = int(
                    float(par.replace("Attenuation", "").replace("dB", "").strip())
                )
            elif "CenterField" in par:
                params["center_field"] = float(
                    par.replace("CenterField", "").replace("G", "").strip()
                )
            elif "ConvTime" in par:
                params["conversion_time"] = float(
                    par.replace("ConvTime", "").replace("ms", "").strip()
                )
            elif "TimeConst" in par:
                params["time_constant"] = float(
                    par.replace("TimeConst", "").replace("ms", "").strip()
                )
            elif "ModAmp" in par:
                params["modulation_amplitude"] = float(
                    par.replace("ModAmp", "").replace("G", "").strip()
                )
            elif "ModFreq" in par:
                params["modulation_frequency"] = float(
                    par.replace("ModFreq", "").replace("kHz", "").strip()
                )
            elif "NbScansDone" in par:
                params["nscans"] = int(par.replace("NbScansDone", "").strip())
            elif "Temperature" in par:
                params["temperature"] = float(
                    par.replace("Temperature", "").replace("K", "").strip()
                )
            elif "XNAM" in par:
                sweep_domain = par.replace("XNAM", "").replace("'", "").strip()
            elif "XUNI" in par:
                params["x_unit"] = par.replace("XUNI", "").replace("'", "").strip()
            elif "XPTS" in par:
                params["x_points"] = int(par.replace("XPTS", "").strip())
            elif "XMIN" in par:
                params["x_min"] = float(par.replace("XMIN", "").strip())
            elif "XWID" in par:
                params["x_width"] = float(par.replace("XWID", "").strip())
            elif "XTYP" in par:
                xtyp = par.replace("XTYP", "").strip()
                if xtyp == "IGD":
                    params["x_type"] = "nonlinear"
                elif xtyp == "IDX":
                    params["x_type"] = "linear"

            elif "XFMT" in par:
                dfmt = par.replace("XFMT", "").strip()
                if dfmt == "D":
                    params["x_format"] = "float64"
                elif dfmt == "F":
                    params["x_format"] = "float32"
                elif dfmt == "C":
                    params["x_format"] = "int8"
                elif dfmt == "S":
                    params["x_format"] = "int16"
                elif dfmt == "I":
                    params["x_format"] = "int32"

            elif "YNAM" in par:
                y_domain = par.replace("YNAM", "").replace("'", "").strip()
            elif "YUNI" in par:
                params["y_unit"] = par.replace("YUNI", "").replace("'", "").strip()
            elif "YPTS" in par:
                params["y_points"] = int(par.replace("YPTS", "").strip())
            elif "YMIN" in par:
                params["y_min"] = float(par.replace("YMIN", "").strip())
            elif "YWID" in par:
                params["y_max"] = float(par.replace("YWID", "").strip())
            elif "YTYP" in par:
                ytyp = par.replace("YTYP", "").strip()
                if ytyp == "IGD":
                    params["y_type"] = "nonlinear"
                elif ytyp == "IDX":
                    params["y_type"] = "linear"

            elif "YFMT" in par:
                dfmt = par.replace("YFMT", "").strip()
                if dfmt == "D":
                    params["y_format"] = "float64"
                elif dfmt == "F":
                    params["y_format"] = "float32"
                elif dfmt == "C":
                    params["y_format"] = "int8"
                elif dfmt == "S":
                    params["y_format"] = "int16"
                elif dfmt == "I":
                    params["y_format"] = "int32"

            elif "ZNAM" in par:
                z_domain = par.replace("ZNAM", "").replace("'", "").strip()
            elif "ZUNI" in par:
                params["z_unit"] = par.replace("ZUNI", "").replace("'", "").strip()
            elif "ZPTS" in par:
                params["z_points"] = int(par.replace("ZPTS", "").strip())
            elif "ZMIN" in par:
                params["z_min"] = float(par.replace("ZMIN", "").strip())
            elif "ZWID" in par:
                params["z_max"] = float(par.replace("ZWID", "").strip())
            elif "ZTYP" in par:
                ztyp = par.replace("ZTYP", "").strip()
                if ztyp == "IGD":
                    params["z_type"] = "nonlinear"
                elif ztyp == "IDX":
                    params["z_type"] = "linear"

            elif "ZFMT" in par:
                dfmt = par.replace("ZFMT", "").strip()
                if dfmt == "D":
                    params["z_format"] = "float64"
                elif dfmt == "F":
                    params["z_format"] = "float32"
                elif dfmt == "C":
                    params["z_format"] = "int8"
                elif dfmt == "S":
                    params["z_format"] = "int16"
                elif dfmt == "I":
                    params["z_format"] = "int32"

            elif "IRFMT" in par:
                dfmt = par.replace("IRFMT", "").strip()
                if dfmt == "D":
                    params["real_format"] = "float64"
                elif dfmt == "F":
                    params["real_format"] = "float32"
                elif dfmt == "C":
                    params["real_format"] = "int8"
                elif dfmt == "S":
                    params["real_format"] = "int16"
                elif dfmt == "I":
                    params["real_format"] = "int32"

            elif "IIFMT" in par:
                ifmt = par.replace("IIFMT", "").strip()
                if ifmt == "D":
                    params["imag_format"] = "float64"
                elif ifmt == "F":
                    params["imag_format"] = "float32"
                elif ifmt == "C":
                    params["imag_format"] = "int8"
                elif ifmt == "S":
                    params["imag_format"] = "int16"
                elif ifmt == "I":
                    params["imag_format"] = "int32"

            elif "IKKF" in par:
                params["data_type"] = par.replace("IKKF", "").strip()
            elif "BSEQ" in par:
                params["endian"] = par.replace("BSEQ", "").strip()

        except ValueError:
            continue

    if sweep_domain == "Time" and int(params["attenuation"]) == 60:
        params.pop("attenuation", None)
        params.pop("power", None)
    elif sweep_domain == "Time" and int(params["pulse_attenuation"]) == 60:
        params["pulse_attenuation"] = params["attenuation"]
        params.pop("attenuation", None)

    if all(
        [
            y not in params.keys()
            for y in ["real_format", "imag_format", "x_format", "y_format", "z_format"]
        ]
    ):
        raise TypeError(
            "data format unknown. None of IRFMT, IIFMT, XFMT, YFMT, ZFMT found in DSC"
        )

    for x in ["x", "y", "z"]:
        if x + "_format" not in params.keys():
            if params["data_type"] == "REAL":
                params[x + "_format"] = params["real_format"]
            elif params["data_type"] == "CPLX":
                params[x + "_format"] = params["imag_format"]

    return params


def load_dta(path_dta, path_xgf=None, path_ygf=None, path_zgf=None, params={}):
    """
    Import data from .DTA file. Uses .DSC and .XGF, .YGF, or .ZGF files if they exists

    Args:
        path_dta (str) : Path to .DTA file
        path_xgf (str) : path to .XGF file for 1D data with nonlinear axis, "none" otherwise
        path_ygf (str) : path to .YGF file for 2D data, "none" if 1D or linear y axis
        path_zgf (str) : path to .ZGF file for 3D data, "none" if 1D/2D or linear z axis
        params (dict) : dictionary of parameters

    Returns:
        abscissa (ndarray) : coordinates for spectrum or spectra
        spec (ndarray) : spectrum for 1D or spectra for 2D
        params (dict) : updated dictionary of parameters
        dims (list) : dimensions
    """

    dta_dtype = np.dtype(params["x_format"]).newbyteorder(params["endian"])
    file_opened = open(path_dta, "rb")
    file_bytes = file_opened.read()
    file_opened.close()
    spec = np.frombuffer(file_bytes, dtype=dta_dtype)
    if params["x_type"] == "nonlinear":
        abscissa = [
            load_gf_files(
                path_xgf,
                axis_type=params["x_type"],
                axis_format=params["x_format"],
                axis_points=params["x_points"],
                axis_min=params["x_min"],
                axis_max=params["x_max"],
                endian=params["endian"],
            )
        ]
    elif params["x_type"] == "linear":
        abscissa = [
            np.linspace(
                params["x_min"],
                params["x_min"] + params["x_width"],
                params["x_points"],
            )
        ]

    if "x_unit" not in params.keys() or params["x_unit"] in [
        "s",
        "ms",
        "ns",
        "ps",
        "Time",
        "time",
    ]:
        dims = ["t2"]
    elif params["x_unit"] in ["G", "mT", "T", "Field", "field"]:
        if params["x_unit"] == "G":
            abscissa = [x / 10 for x in abscissa]
        elif params["x_unit"] == "T":
            abscissa = [x * 1000 for x in abscissa]
        dims = ["B0"]
    else:
        dims = ["t2"]

    if params["data_type"] == "CPLX":
        spec = spec.astype(dtype=params["imag_format"]).view(dtype=np.dtype("complex"))
    elif params["data_type"] == "REAL":
        spec = spec.astype(dtype=params["real_format"]).view()

    if (
        "z_points" not in params.keys()
        or ("z_points" in params.keys() and params["z_points"] == 1)
    ) and ("y_points" in params.keys() and params["y_points"] != 1):
        spec = np.reshape(spec, (params["x_points"], params["y_points"]), order="F")

        if "y_unit" in params.keys():
            dims.append(params["y_unit"])
        else:
            dims.append("t1")

        abscissa.append(
            load_gf_files(
                path_ygf,
                axis_type=params["y_type"],
                axis_format=params["y_format"],
                axis_points=params["y_points"],
                axis_min=params["y_min"],
                axis_max=params["y_max"],
                endian=params["endian"],
            )
        )

    elif "z_points" in params.keys() and params["z_points"] != 1:
        spec = np.reshape(
            spec,
            (params["x_points"], params["y_points"]),
            params["z_points"],
            order="F",
        )

        if "y_unit" in params.keys():
            dims.append(params["y_unit"])
        else:
            dims.append("t1")

        if "z_unit" in params.keys():
            dims.append(params["z_unit"])
        else:
            dims.append("t0")

        abscissa.append(
            load_gf_files(
                path_ygf,
                axis_type=params["y_type"],
                axis_format=params["y_format"],
                axis_points=params["y_points"],
                axis_min=params["y_min"],
                axis_max=params["y_max"],
                endian=params["endian"],
            )
        )
        abscissa.append(
            load_gf_files(
                path_zgf,
                axis_type=params["z_type"],
                axis_format=params["z_format"],
                axis_points=params["z_points"],
                axis_min=params["z_min"],
                axis_max=params["z_max"],
                endian=params["endian"],
            )
        )

    params.pop("endian", None)
    params.pop("x_format", None)
    params.pop("y_format", None)
    params.pop("z_format", None)
    params.pop("real_format", None)
    params.pop("imag_format", None)
    params.pop("data_type", None)

    return spec, abscissa, dims, params


def load_gf_files(
    path, axis_type="", axis_format="", axis_points=1, axis_min=1, axis_max=1, endian=""
):
    """
    Import data from .XGF, .YGF, or .ZGF files

    Args:
        path (str) : Path to ._GF file
        axis_type (str) : linear or nonlinear
        axis_format (str) : format of file data
        axis_points (int) : number of points in axis
        axis_min (float) : minimum value of axis
        axis_max (float) : maximum value of axis
        endian (float) : endian of data

    Returns:
        abscissa (ndarray) : axis coordinates
    """

    if path != "none":
        if axis_type == "linear":
            warnings.warn("axis format is linear, confirm that the axis is correct")
        gf_type = np.dtype(axis_format).newbyteorder(endian)
        file_opened = open(path, "rb")
        file_bytes = file_opened.read()
        file_opened.close()
        abscissa = np.frombuffer(file_bytes, dtype=axis_format)
    elif path == "none":
        if axis_type == "nonlinear":
            warnings.warn(
                "axis is nonlinear, confirm that file is not needed and the axis is correct"
            )
        abscissa = np.linspace(axis_min, axis_max, axis_points)
    else:
        warnings.warn("axis format not supported, axis is only indexed")
        abscissa = [range(axis_points)]

    return abscissa
