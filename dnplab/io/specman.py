import numpy as np
import os
from .. import DNPData


def import_specman(path):
    """
    Import specman data and return DNPData object

    Args:
        path (str) : Path to either .d01 or .exp file

    Returns:
        specman_data (object) : DNPData object containing specman data
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

    params = load_specman_exp(file_name_exp)
    values, dims, coords, attrs = load_specman_d01(file_name_d01, params)

    specman_data = DNPData(values, dims, coords, attrs)

    return specman_data


def load_specman_exp(path):
    """
    Import parameter fields of specman data

    Args:
        path (str) : Path to either .d01 or .exp file

    Returns:
        params (dict) : dictionary of parameter fields and values
    """
    exp_file_opened = open(path, encoding="utf8", errors="ignore")
    file_contents = exp_file_opened.read().splitlines()
    exp_file_opened.close()
    params = {}
    c = ""
    for i in range(0, len(file_contents)):
        exp_content = str(file_contents[i])
        splt_exp_content = exp_content.split(" = ")
        if "[" in exp_content and "]" in exp_content and "=" not in exp_content:
            c = splt_exp_content[0].replace("[", "").replace("]", "")
        elif exp_content == "":
            c = "param"
        elif len(splt_exp_content) > 1:
            params[c + "_" + splt_exp_content[0]] = splt_exp_content[1]
        elif len(splt_exp_content) == 1 and exp_content != "":
            params[c + "_" + str(i)] = splt_exp_content
        else:
            pass

    return params


def load_specman_d01(path, params):
    """
    Import spectrum or spectra of specman data

    Args:
        path (str) : Path to either .d01 or .exp file
        params (dict) : dictionary of parameters from exp file

    Returns:
        abscissa (list) : coordinates of axes
        y_data (ndarray) : spectrum or spectra if >1D
        dims (list) : axes names
        params (dict) : updated parameters dictionary
    """

    if not params:
        params = {}

    file_opened = open(path, "rb")
    uint_read = np.fromfile(file_opened, dtype=np.uint32)
    file_opened.close()

    file_opened = open(path, "rb")
    float_read = np.fromfile(file_opened, dtype="<f4")
    file_opened.close()
    float_data_real = float_read[14 : uint_read[7] + 14]
    float_data_complex = float_read[uint_read[7] + 14 : len(float_read)]
    float_data_folded = float_data_real + 1j * float_data_complex

    if uint_read[2] == 1:
        y_data = np.reshape(float_data_folded, (uint_read[3]))
    elif uint_read[2] == 2:
        y_data = np.reshape(float_data_folded, (uint_read[4], uint_read[3]))
    elif uint_read[2] == 3:
        y_data = np.reshape(
            float_data_folded, (uint_read[5], uint_read[4], uint_read[3])
        )
    elif uint_read[2] == 4:
        y_data = np.reshape(
            float_data_folded, (uint_read[6], uint_read[5], uint_read[4], uint_read[3])
        )
    else:
        raise TypeError("DNPLab currently only supports up to 4D data")

    y_data = np.transpose(y_data)

    dims_full = ["t2", "t1", "t0", "t"]
    dims = dims_full[0 : uint_read[2]]
    axes_lengths = uint_read[9:13]

    abscissa = []
    for k in range(0, len(dims)):
        if dims[k] in params.keys():
            abscissa.append(params[dims[k]])
        else:
            abscissa.append(np.array(range(0, axes_lengths[k])))

    return y_data, dims, abscissa, params
