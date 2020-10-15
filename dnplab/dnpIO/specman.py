import numpy as np
import os
from dnplab import dnpdata


def import_specman(path):
    """
    Import specman data and return dnpdata object

    Args:
        path (str) : Path to either .d01 or .exp file

    Returns:
        specman_data (object) : dnpdata object containing specman data
    """
    if path[-1] == os.sep:
        path = path[:-1]
    if path[-1] == "p":
        file_name_d01 = path.replace("exp", "d01")
        file_name_exp = path
    elif path[-1] == "1":
        file_name_exp = path.replace("d01", "exp")
        file_name_d01 = path

    attrs = load_specman_exp(file_name_exp)
    values, dims, data_lengths = load_specman_values(file_name_d01)
    coords = load_specman_coords(attrs, data_lengths, dims)

    specman_data = dnpdata(values, coords, dims, attrs)

    return specman_data


def load_specman_coords(attrs, data_lengths, dims):
    """
    Import axes coordinates of specman data

    Args:
        attrs (dict) : dictionary of parameter fields and values
        data_lengths (ndarray) : axes lengths
        dims (list) : axes names

    Returns:
        abscissa (list) : coordinates of axes
    """
    abscissa = []
    for k in range(0, len(dims) - 1):
        if dims[k] in attrs.keys():
            abscissa.append(attrs[dims[k]])
        else:
            abscissa.append(np.array(range(0, data_lengths[k])))

    return abscissa


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

    return params


def load_specman_values(path):
    """
    Import spectrum or spectra of specman data

    Args:
        path (str) : Path to either .d01 or .exp file

    Returns:
        y_data (ndarray) : spectrum or spectra if >1D
        dims (list) : axes names
        axes_lengths (ndarray) : axes lengths
    """
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
    y_data = np.transpose(y_data)

    dims_full = ["t", "x", "y", "z"]
    dims = dims_full[0 : uint_read[2]]
    axes_lengths = uint_read[9:13]

    return y_data, dims, axes_lengths
