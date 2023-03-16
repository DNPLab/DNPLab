import numpy as _np
import _os
import dnplab as _dnp


def import_specman(path):
    """Import SpecMan data and return DNPData object

    DNPLab function to import SpecMan4EPR data (https://specman4epr.com/). The function returns a DNPdata object with the spectra data.
    The structure of the DNPdata object can be complex and the variables saved by SpecMan depend on the individual spectrometer configuration. Therefore, the import function returns a numpy array with the dimension "x0", "x1", "x2", "x3", "x4". In any case, the dimension "x0" corresponds to the variables stored in the data file. The spectroscopic data is stored in "x1" to "x4", depending on how many dimensions were recorded.
    The import function will require a parser script to properly assign the spectroscopic data and proper coordinates.

    Args:
        path (str) : Path to either .d01 or .exp file

    Returns:
        data (DNPData) : DNPData object containing SpecMan EPR data
    """

    if path[-1] == _os.sep:
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

    # Assign data/spectrum type
    attrs["experiment_type"] = "epr_spectrum"

    specman_data = _dnp.DNPData(data, dims, coords, attrs)

    return specman_data


def load_specman_exp(path):
    """Import SpecMan parameter fields

    Args:
        path (str) : Path to either .d01 or .exp file

    Returns:
        attrs (dict) : dictionary of parameter fields and values (DNPLab attributes)

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

    Args:
        path (str) : Path to either .d01 or .exp file

    Returns:
        data (ndarray) :    SpecMan data as numpy array
        params (dict) :     Dictionary with import updated parameters dictionary
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

    # Number of dimensions of stored data
    attrs["dims"] = uint_read[2]

    dataShape = uint_read[2 : 2 + uint_read[2] * 6]

    attrs["dataStreamShape"] = dataShape

    attrs["dataStartIndex"] = 2 + attrs["numberOfVariables"] * 6 - 1

    if verbose == True:
        print("Data paramters:")
        print(attrs)

    data = float_read[attrs["dataStartIndex"] : -1]

    if attrs["dims"] == 1:
        data = _np.reshape(data, (attrs["numberOfVariables"], uint_read[3]))

    elif attrs["dims"] == 2:
        data = _np.reshape(data, (uint_read[2], uint_read[3], uint_read[4]))

    elif attrs["dims"] == 3:
        data = _np.reshape(
            data, (uint_read[2], uint_read[3], uint_read[4], uint_read[5])
        )

    elif attrs["dims"] == 4:
        data = _np.reshape(
            data, (uint_read[2], uint_read[3], uint_read[4], uint_read[5], uint_read[6])
        )

    elif attrs["dims"] >= 4:
        print("Maximum dimensionality for SpecMan data is 4D")
        return None

    # Swap first axis with last
    data = _np.swapaxes(data, 0, -1)

    # SpecMan data can have a maximum of four dimensions
    dims_full = ["x0", "x1", "x2", "x3", "x4"]
    dims = dims_full[0 : dataShape[0] + 1]

    coords = []

    shape = _np.shape(data)

    for index in range(data.ndim):
        coords.append(_np.arange(0, shape[index]))

    return data, dims, coords, attrs
