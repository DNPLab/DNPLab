"""Functions to import Bruker EPR data"""

import numpy as _np
import os
from .. import DNPData
import warnings


def import_bes3t(path):
    """
    Import Bruker BES3T data and return dnpdata object

    Args:
        path (str) : Path to either .DSC or .DTA file

    Returns:
        bes3t_data (object) : DNPData object containing Bruker BES3T data
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

    attrs = load_dsc(path_dsc)

    values, dims, coords, attrs = load_dta(
        path_dta, path_xgf, path_ygf, path_zgf, attrs
    )
    attrs["spectrometer_format"] = "xepr"
    attrs["experiment_type"] = "epr_spectrum"
    attrs["nrScans"] = attrs["nscans"]

    bes3t_data = DNPData(values, dims, coords, attrs)

    return bes3t_data


def load_dsc(path):
    """
    Import contents of .DSC file

    Args:
        path (str) : Path to .DSC file

    Returns:
        attrs (dict) : dictionary of parameters
    """

    file_opened = open(path, "r")
    dscfile_contents = file_opened.readlines()
    file_opened.close()

    attrs = {}
    sweep_domain = None
    for ix in range(len(dscfile_contents)):
        try:
            par = dscfile_contents[ix].rstrip("\t").rstrip("\n")
            if "MWFQ" in par:
                attrs["frequency"] = float(par.replace("MWFQ", "").strip()) / 1e9
            elif "Power" in par and "Atten" not in par:
                attrs["power"] = float(
                    par.replace("Power", "").replace("mW", "").strip()
                )
            elif "PowerAtten" in par:
                attrs["attenuation"] = int(
                    float(par.replace("PowerAtten", "").replace("dB", "").strip())
                )
            elif "Attenuation" in par:
                attrs["pulse_attenuation"] = int(
                    float(par.replace("Attenuation", "").replace("dB", "").strip())
                )
            elif "CenterField" in par:
                attrs["center_field"] = float(
                    par.replace("CenterField", "").replace("G", "").strip()
                )
            elif "ConvTime" in par:
                attrs["conversion_time"] = float(
                    par.replace("ConvTime", "").replace("ms", "").strip()
                )
            elif "TimeConst" in par:
                attrs["time_constant"] = float(
                    par.replace("TimeConst", "").replace("ms", "").strip()
                )
            elif "ModAmp" in par:
                attrs["modulation_amplitude"] = float(
                    par.replace("ModAmp", "").replace("G", "").strip()
                )
            elif "ModFreq" in par:
                attrs["modulation_frequency"] = float(
                    par.replace("ModFreq", "").replace("kHz", "").strip()
                )
            elif "NbScansDone" in par:
                attrs["nscans"] = int(par.replace("NbScansDone", "").strip())
            elif "Temperature" in par:
                attrs["temperature"] = float(
                    par.replace("Temperature", "").replace("K", "").strip()
                )
            elif "XNAM" in par:
                sweep_domain = par.replace("XNAM", "").replace("'", "").strip()
            elif "XUNI" in par:
                attrs["x_unit"] = par.replace("XUNI", "").replace("'", "").strip()
            elif "XPTS" in par:
                attrs["x_points"] = int(par.replace("XPTS", "").strip())
            elif "XMIN" in par:
                attrs["x_min"] = float(par.replace("XMIN", "").strip())
            elif "XWID" in par:
                attrs["x_width"] = float(par.replace("XWID", "").strip())
            elif "XTYP" in par:
                xtyp = par.replace("XTYP", "").strip()
                if xtyp == "IGD":
                    attrs["x_type"] = "nonlinear"
                elif xtyp == "IDX":
                    attrs["x_type"] = "linear"

            elif "XFMT" in par:
                attrs["x_format"] = _return_data_type(par, "XFMT")

            elif "YNAM" in par:
                y_domain = par.replace("YNAM", "").replace("'", "").strip()
            elif "YUNI" in par:
                attrs["y_unit"] = par.replace("YUNI", "").replace("'", "").strip()
            elif "YPTS" in par:
                attrs["y_points"] = int(par.replace("YPTS", "").strip())
            elif "YMIN" in par:
                attrs["y_min"] = float(par.replace("YMIN", "").strip())
            elif "YWID" in par:
                attrs["y_width"] = float(par.replace("YWID", "").strip())
            elif "YTYP" in par:
                ytyp = par.replace("YTYP", "").strip()
                if ytyp == "IGD":
                    attrs["y_type"] = "nonlinear"
                elif ytyp == "IDX":
                    attrs["y_type"] = "linear"

            elif "YFMT" in par:
                attrs["y_format"] = _return_data_type(par, "YFMT")

            elif "ZNAM" in par:
                z_domain = par.replace("ZNAM", "").replace("'", "").strip()
            elif "ZUNI" in par:
                attrs["z_unit"] = par.replace("ZUNI", "").replace("'", "").strip()
            elif "ZPTS" in par:
                attrs["z_points"] = int(par.replace("ZPTS", "").strip())
            elif "ZMIN" in par:
                attrs["z_min"] = float(par.replace("ZMIN", "").strip())
            elif "ZWID" in par:
                attrs["z_width"] = float(par.replace("ZWID", "").strip())
            elif "ZTYP" in par:
                ztyp = par.replace("ZTYP", "").strip()
                if ztyp == "IGD":
                    attrs["z_type"] = "nonlinear"
                elif ztyp == "IDX":
                    attrs["z_type"] = "linear"

            elif "ZFMT" in par:
                attrs["z_format"] = _return_data_type(par, "ZFMT")

            elif "IRFMT" in par:
                attrs["real_format"] = _return_data_type(par, "IRFMT")

            elif "IIFMT" in par:
                attrs["imag_format"] = _return_data_type(par, "IIFMT")

            elif "IKKF" in par:
                attrs["data_type"] = par.replace("IKKF", "").strip()
            elif "BSEQ" in par:
                attrs["endian"] = par.replace("BSEQ", "").strip()

        except ValueError:
            continue

    if sweep_domain == "Time" and int(attrs["attenuation"]) == 60:
        attrs.pop("attenuation", None)
        attrs.pop("power", None)
    elif sweep_domain == "Time" and int(attrs["pulse_attenuation"]) == 60:
        attrs["pulse_attenuation"] = attrs["attenuation"]
        attrs.pop("attenuation", None)

    if all(
        [
            y not in attrs.keys()
            for y in ["real_format", "imag_format", "x_format", "y_format", "z_format"]
        ]
    ):
        raise TypeError(
            "data format unknown. None of IRFMT, IIFMT, XFMT, YFMT, ZFMT found in DSC"
        )

    for x in ["x", "y", "z"]:
        if x + "_format" not in attrs.keys():
            if attrs["data_type"] == "REAL":
                attrs[x + "_format"] = attrs["real_format"]
            elif attrs["data_type"] == "CPLX":
                attrs[x + "_format"] = attrs["imag_format"]

    return attrs


def load_dta(path_dta, path_xgf=None, path_ygf=None, path_zgf=None, attrs={}):
    """
    Import data from .DTA file. Uses .DSC and .XGF, .YGF, or .ZGF files if they exists

    Args:
        path_dta (str) : Path to .DTA file
        path_xgf (str) : path to .XGF file for 1D data with nonlinear axis, "none" otherwise
        path_ygf (str) : path to .YGF file for 2D data, "none" if 1D or linear y axis
        path_zgf (str) : path to .ZGF file for 3D data, "none" if 1D/2D or linear z axis
        attrs (dict) : dictionary of parameters

    Returns:
        values (ndarray) : Spectrum for 1D or spectra for 2D
        dims (list) : dimensions
        coords (ndarray) : coordinates for spectrum or spectra
        attrs (dict) : updated dictionary of parameters
    """

    dta_dtype = _np.dtype(attrs["x_format"]).newbyteorder(attrs["endian"])
    file_opened = open(path_dta, "rb")
    file_bytes = file_opened.read()
    file_opened.close()
    values = _np.frombuffer(file_bytes, dtype=dta_dtype)
    if attrs["x_type"] == "nonlinear":
        coords = [
            load_gf_files(
                path_xgf,
                axis_type=attrs["x_type"],
                axis_format=attrs["x_format"],
                axis_points=attrs["x_points"],
                axis_min=attrs["x_min"],
                axis_width=attrs["x_width"],
                endian=attrs["endian"],
            )
        ]
    elif attrs["x_type"] == "linear":
        coords = [
            _np.linspace(
                attrs["x_min"],
                attrs["x_min"] + attrs["x_width"],
                attrs["x_points"],
            )
        ]

    if "x_unit" in attrs.keys() and attrs["x_unit"] in ["G", "T"]:
        if attrs["x_unit"] == "G":
            coords = [x / 10 for x in coords]
        elif attrs["x_unit"] == "T":
            coords = [x * 1000 for x in coords]
        dims = ["B0"]
    else:
        dims = ["t2"]

    if attrs["data_type"] == "CPLX":
        values = values.astype(dtype=attrs["imag_format"]).view(
            dtype=_np.dtype("complex")
        )
    elif attrs["data_type"] == "REAL":
        values = values.astype(dtype=attrs["real_format"]).view()

    if (
        "z_points" not in attrs.keys()
        or ("z_points" in attrs.keys() and attrs["z_points"] == 1)
    ) and ("y_points" in attrs.keys() and attrs["y_points"] != 1):
        values = _np.reshape(values, (attrs["x_points"], attrs["y_points"]), order="F")

        dims.append("t1")

        coords.append(
            load_gf_files(
                path_ygf,
                axis_type=attrs["y_type"],
                axis_format=attrs["y_format"],
                axis_points=attrs["y_points"],
                axis_min=attrs["y_min"],
                axis_width=attrs["y_width"],
                endian=attrs["endian"],
            )
        )

    elif "z_points" in attrs.keys() and attrs["z_points"] != 1:
        values = _np.reshape(
            values,
            (attrs["x_points"], attrs["y_points"]),
            attrs["z_points"],
            order="F",
        )

        coords.append(
            load_gf_files(
                path_ygf,
                axis_type=attrs["y_type"],
                axis_format=attrs["y_format"],
                axis_points=attrs["y_points"],
                axis_min=attrs["y_min"],
                axis_width=attrs["y_width"],
                endian=attrs["endian"],
            )
        )
        dims.append("t1")

        coords.append(
            load_gf_files(
                path_zgf,
                axis_type=attrs["z_type"],
                axis_format=attrs["z_format"],
                axis_points=attrs["z_points"],
                axis_min=attrs["z_min"],
                axis_width=attrs["z_width"],
                endian=attrs["endian"],
            )
        )
        dims.append("t0")

    attrs = {
        x: attrs[x]
        for x in attrs.keys()
        if x
        not in [
            "endian",
            "x_format",
            "x_type",
            "x_points",
            "x_min",
            "x_width",
            "y_format",
            "y_type",
            "y_points",
            "y_min",
            "y_width",
            "z_format",
            "z_type",
            "z_points",
            "z_min",
            "z_width",
            "real_format",
            "imag_format",
            "data_type",
        ]
    }

    return values, dims, coords, attrs


def load_gf_files(
    path,
    axis_type="",
    axis_format="",
    axis_points=1,
    axis_min=1,
    axis_width=1,
    endian="",
):
    """
    Import data from .XGF, .YGF, or .ZGF files

    Args:
        path (str) : Path to ._GF file
        axis_type (str) : linear or nonlinear
        axis_format (str) : format of file data
        axis_points (int) : number of points in axis
        axis_min (float) : minimum value of axis
        axis_width (float) : total width of axis
        endian (float) : endian of data

    Returns:
        coords (ndarray) : axis coordinates
    """

    if path != "none":
        if axis_type == "linear":
            warnings.warn("axis format is linear, confirm that the axis is correct")
        gf_type = _np.dtype(axis_format).newbyteorder(endian)
        file_opened = open(path, "rb")
        file_bytes = file_opened.read()
        file_opened.close()
        coords = _np.frombuffer(file_bytes, dtype=axis_format)
    elif path == "none":
        if axis_type == "nonlinear":
            warnings.warn(
                "axis is nonlinear, confirm that file is not needed and the axis is correct"
            )
        coords = _np.linspace(axis_min, axis_min + axis_width, axis_points)
    else:
        warnings.warn("axis format not supported, axis is only indexed")
        coords = _np.array([range(axis_points)])

    return coords


def _return_data_type(par, key):
    fmt = par.replace(key, "").strip()
    if fmt == "D":
        return "float64"
    elif fmt == "F":
        return "float32"
    elif fmt == "C":
        return "int8"
    elif fmt == "S":
        return "int16"
    elif fmt == "I":
        return "int32"
