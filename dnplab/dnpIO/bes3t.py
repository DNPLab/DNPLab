import numpy as np
import os
from dnplab import dnpdata
import warnings


def import_bes3t(path):
    """
    Import Bruker BES3T data and return dnpdata object

    Args:
        path (str) : Path to either .DSC or .DTA file

    Returns:
        bes3t_data (object) : dnpdata object containing Bruker BES3T data
    """

    pathexten = path[-3:]
    path = path[:-3]
    path_ygf = "none"
    if pathexten == "DSC" or pathexten == "DTA":
        path_dsc = path + "DSC"
        path_dta = path + "DTA"
        if os.path.isfile(path + "YGF"):
            path_ygf = path + "YGF"

    elif pathexten == "YGF":
        path_ygf = path + "YGF"
        path_dsc = path + "DSC"
        path_dta = path + "DTA"

    else:
        raise TypeError("data file must be .DTA, .DSC, or .YGF")

    params = load_bes3t_dsc(path_dsc)
    coords, values, attrs, dims = load_bes3t_dta(path_dta, path_ygf, params)

    bes3t_data = dnpdata(values, coords, dims, attrs)

    return bes3t_data


def load_bes3t_dsc(path):
    """
    Import contents of .DSC file

    Args:
        path (str) : Path to .DSC file

    Returns:
        params (dict) : dictionary of parameters
    """

    file_opened = open(path, "r")
    dscfile_contents = file_opened.readlines()

    params = {}
    for ix in range(len(dscfile_contents)):
        par = dscfile_contents[ix].rstrip("\t").rstrip("\n")
        if "MWFQ" in par:
            params["frequency"] = float(par.replace("MWFQ", "").strip()) / 1e9
        elif "Power" in par and "Atten" not in par:
            params["power"] = float(par.replace("Power", "").replace("mW", "").strip())
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
                params["sweep_type"] = "nonlinear"
            elif xtyp == "IDX":
                params["sweep_type"] = "linear"

        elif "IRFMT" in par:
            dfmt = par.replace("IRFMT", "").strip()
            if dfmt == "D":
                params["sweep_format"] = "float64"
            elif dfmt == "F":
                params["sweep_format"] = "float32"
            elif dfmt == "C":
                params["sweep_format"] = "int8"
            elif dfmt == "S":
                params["sweep_format"] = "int16"
            elif dfmt == "I":
                params["sweep_format"] = "int32"

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
                params["slice_type"] = "nonlinear"
            elif ytyp == "IDX":
                params["slice_type"] = "linear"

        elif "YFMT" in par:
            yfmt = par.replace("YFMT", "").strip()
            if yfmt == "D":
                params["slice_format"] = "float64"
            elif yfmt == "F":
                params["slice_format"] = "float32"
            elif yfmt == "C":
                params["slice_format"] = "int8"
            elif yfmt == "S":
                params["slice_format"] = "int16"
            elif yfmt == "I":
                params["slice_format"] = "int32"

        elif "IKKF" in par:
            params["data_type"] = par.replace("IKKF", "").strip()
        elif "BSEQ" in par:
            params["endian"] = par.replace("BSEQ", "").strip()

    file_opened.close()

    if sweep_domain == "Time" and int(params["attenuation"]) == 60:
        params.pop("attenuation", None)
        params.pop("power", None)
    elif sweep_domain == "Time" and int(params["pulse_attenuation"]) == 60:
        params["pulse_attenuation"] = params["attenuation"]
        params.pop("attenuation", None)

    return params


def load_bes3t_dta(path_dta, path_ygf, params):
    """
    Import data from .DTA file, and from .YGF file if one exists

    Args:
        path_dta (str) : Path to .DTA file
        path_ygf (str) : Path to .YGF file if 2D, otherwise is 'none'

    Returns:
        abscissa (ndarray) : coordinates for spectrum or spectra
        spec (ndarray) : spectrum for 1D or spectra for 2D
        params (dict) : updated dictionary of parameters
        dims (list) : dimensions
    """

    dta_dtype = np.dtype(params["sweep_format"]).newbyteorder(params["endian"])
    file_opened = open(path_dta, "rb")
    file_bytes = file_opened.read()
    spec = np.frombuffer(file_bytes, dtype=dta_dtype)
    abscissa = [
        np.linspace(
            params["x_min"],
            params["x_min"] + params["x_width"],
            params["x_points"],
        )
    ]
    dims = [params["x_unit"]]

    if params["data_type"] == "CPLX":
        spec = spec.astype(dtype=params["sweep_format"]).view(dtype=np.dtype("complex"))

    if "y_points" in params.keys() and params["y_points"] != 1:
        spec = np.reshape(spec, (params["x_points"], params["y_points"]), order="F")
        dims.append(params["y_unit"])
        if path_ygf != "none":
            if params["slice_type"] == "linear":
                warnings.warn("axis is linear, confirm that indirect axis is correct")
            ygf_type = np.dtype(params["slice_format"]).newbyteorder(params["endian"])
            file_opened = open(path_ygf, "rb")
            file_bytes = file_opened.read()
            abscissa.append(np.frombuffer(file_bytes, dtype=ygf_type))
        elif path_ygf == "none":
            if params["slice_type"] == "nonlinear":
                warnings.warn(
                    "axis is nonlinear, confirm that .YGF file is not needed and indirect axis is correct"
                )
            abscissa.append(
                np.linspace(params["y_min"], params["y_max"], params["y_points"])
            )
        else:
            warnings.warn("indirect axis format not supported, axis is only indexed")
            abscissa.append([range(params["y_points"])])

    params.pop("endian", None)
    params.pop("sweep_format", None)
    params.pop("slice_format", None)
    params.pop("data_type", None)

    file_opened.close()

    return abscissa, spec, params, dims
