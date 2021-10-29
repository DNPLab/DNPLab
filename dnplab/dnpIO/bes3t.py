import numpy as np
import os
from .. import dnpdata
import warnings

import re
import io

__all__ = ['import_bes3t', 'load_dsc', 'load_dta', 'load_gf_files',]

rename_dict = {
        'MWFQ': 'frequency',
        'mW': "power",
        'PowerAtten':'attenuation',
        'Attenuation':'pulse_attenuation',
        'CenterField':'center_field',
        'ConvTime':'conversion_time',
        'TimeConst':'time_constant',
        'ModAmp':'modulation_amplitude',
        'ModFreq':'modulation_frequency',
        'NbScansDone':'nscans',
        'Temperature':'temperature',
        'XNAM':'x_name',
        'XUNI':'x_unit',
        'XPTS':'x_points',
        'XMIN':'x_min',
        'XWID':'x_width',
        'XTYP':'x_type',
        'XFMT':'x_format',
        'YUNI':'y_unit',
        'YPTS':'y_points',
        'YMIN':'y_min',
        'YWID':'y_width',
        'YTYP':'y_type',
        'YFMT':'y_format',
        'ZNAM':'z_name',
        'ZUNI':'z_unit',
        'ZPTS':'z_points',
        'ZMIN':'z_min',
        'ZWID':'z_width',
        'ZTYP':'z_type',
        'ZFMT':'z_format',
        'IRFMT':'real_format',
        'IKKF':'data_type',
        'BSEQ':'endian',
        }

float_params = [
        'frequency',
        'power',
        'center_field',
        'conversion_time',
        'time_constant',
        'modulation_amplitude',
        'modulation_frequency',
        'temperature',
        'x_min',
        'y_min',
        'z_min',
        'x_width',
        'y_width',
        'z_width',
        ]

int_params = [
        'attenuation',
        'pulse_attenuation',
        'nscans',
        'x_points',
        'y_points',
        'z_points',
        ]

find_string = r'\n\\n\\'
replace_string = r'\\n\\n\\'

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


    desc_params = {}
    spl_params = {}
    dsl_params = {}
    dvc_params = {}
    with open(path, 'r') as f:

        # Section 1: DESC Descriptor Information
        line = f.readline().strip()
        for line in f:
            if line[0] == '#':
                break
            if line[0] == '*':
                continue

            split_line = line.split('\t', 1)

            if len(split_line) == 2:
                key = split_line[0].strip()
                value = split_line[1].strip()
                desc_params[key] = value

        # Section 2: SPL Standard Parameter Layer
        for line in f:
            if line[0] == '#':
                break
            if line[0] == '*':
                continue

            split_line = line.split(' ', 1)

            if len(split_line) == 2:
                key = split_line[0].strip()
                value = split_line[1].strip()
                spl_params[key] = value

        dsl_lines = f.read()

        # Fixes ^M character sometimes present in pulse program
        lines = re.sub(find_string, replace_string, dsl_lines)

    # Section 3: DSL Device Specific Layer
    with io.StringIO(lines) as f:
        for line in f:
            if line[0] == '#':
                break
            if line[0] == '*':
                continue
            if line[0:4] == '.DVC':
                split_line = line[4:].split(',', 1)
                key = split_line[0].strip()
                value = split_line[1].strip()
                dvc_params[key] = value
                continue

            split_line = line.split(' ', 1)

            if len(split_line) == 2:
                key = split_line[0].strip()
                value = split_line[1].strip()

                if len(value) == 0:
                    continue

                ix = 0
                while value[-1] == '\\':
                     
                    ix+=1
                    line = f.readline().rstrip()
                    value += line


                dsl_params[key] = value

    if 'PlsSPELPrgTxt' in dsl_params:
        dsl_params['PlsSPELPrgTxt'] = dsl_params['PlsSPELPrgTxt'].replace('\\n','\n').replace('\\','')

    if 'PlsSPELGlbTxt' in dsl_params:
        dsl_params['PlsSPELGlbTxt'] = dsl_params['PlsSPELGlbTxt'].replace('\\n','\n').replace('\\','')

    params = {}
    params.update(desc_params)
    params.update(spl_params)
    params.update(dsl_params)
    params.update(dvc_params)

    for key in rename_dict:
        if key in params:
            new_key = rename_dict[key]
            value = params[key].strip()
            split_value = value.split(' ', 1)
            if len(split_value) == 1:
                value = split_value[0]
            elif len(split_value) == 2:
                params[key + '_units'] = split_value[1]
                value = split_value[0]

            if new_key in float_params:
                value = float(value)
            elif new_key in int_params:
                value = int(float(value))

            params[new_key] = value

    params['frequency'] = params['frequency'] / 1e9

    if 'XFMT' in params:
        params["x_format"] = _return_data_type(params['XFMT'], "XFMT")
    if 'YFMT' in params:
        params["y_format"] = _return_data_type(params['YFMT'], "YFMT")
    if 'ZFMT' in params:
        params["z_format"] = _return_data_type(params['ZFMT'], "ZFMT")

    if 'IRFMT' in params:
        params["real_format"] = _return_data_type(params['IRFMT'], "IRFMT")
    if 'IIFMT' in params:
        params["imag_format"] = _return_data_type(params['IIFMT'], "IIFMT")


    if params['XTYP'] == "IGD":
        params["x_type"] = "nonlinear"
    else:
        params["x_type"] = "linear"

    if params['YTYP'] == "IGD":
        params["y_type"] = "nonlinear"
    else:
        params["y_type"] = "linear"

    if params['ZTYP'] == "IGD":
        params["z_type"] = "nonlinear"
    else:
        params["z_type"] = "linear"


    for x in ["x", "y", "z"]:
        if x + "_format" not in params.keys():
            if params["data_type"] == "REAL":
                params[x + "_format"] = params["real_format"]
            elif params["data_type"] == "CPLX":
                params[x + "_format"] = params["imag_format"]

#    print(params)
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

    print(params)
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
                axis_width=params["x_width"],
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

    if "x_unit" in params.keys() and params["x_unit"] in ["G", "T"]:
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

        dims.append("t1")

        abscissa.append(
            load_gf_files(
                path_ygf,
                axis_type=params["y_type"],
                axis_format=params["y_format"],
                axis_points=params["y_points"],
                axis_min=params["y_min"],
                axis_width=params["y_width"],
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

        abscissa.append(
            load_gf_files(
                path_ygf,
                axis_type=params["y_type"],
                axis_format=params["y_format"],
                axis_points=params["y_points"],
                axis_min=params["y_min"],
                axis_width=params["y_width"],
                endian=params["endian"],
            )
        )
        dims.append("t1")

        abscissa.append(
            load_gf_files(
                path_zgf,
                axis_type=params["z_type"],
                axis_format=params["z_format"],
                axis_points=params["z_points"],
                axis_min=params["z_min"],
                axis_width=params["z_width"],
                endian=params["endian"],
            )
        )
        dims.append("t0")

    params = {
        x: params[x]
        for x in params.keys()
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

    return spec, abscissa, dims, params


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
        abscissa = np.linspace(axis_min, axis_min + axis_width, axis_points)
    else:
        warnings.warn("axis format not supported, axis is only indexed")
        abscissa = np.array([range(axis_points)])

    return abscissa


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
