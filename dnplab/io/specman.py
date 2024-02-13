import numpy as _np
import os
import dnplab as _dnp
import re

def import_specman(path):
    """Import SpecMan data and return DNPData object

    DNPLab function to import SpecMan4EPR data (https://specman4epr.com/). The function returns a DNPdata object with the spectral data.

    The structure of the DNPdata object can be complex and the variables saved by SpecMan depend on the individual spectrometer configuration. Therefore, the import function returns a numpy array with the dimension "x0", "x1", "x2", "x3", "x4". In any case, the dimension "x0" corresponds to the variables stored in the data file. The spectroscopic data is stored in "x1" to "x4", depending on how many dimensions were recorded.
    The import function will require a parser script to properly assign the spectroscopic data and proper coordinates.

    Args:
        path (str):         Path to either .exp file

    Returns:
        data (DNPData):     DNPData object containing SpecMan EPR data
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

    attrs = load_specman_exp(file_name_exp)
    data, dims, coords, attrs = load_specman_d01(file_name_d01, attrs)

    # Add import path
    attrs["import_path"] = path

    # Assign data/spectrum type
    attrs["experiment_type"] = "epr_spectrum"

    specman_data = _dnp.DNPData(data, dims, coords, attrs)

    return specman_data


def load_specman_exp(path):
    """Import SpecMan parameters

    DNPLab function to read and import the SpecMan exp file. The .exp file is a text file that stores the experimental data, the pulse program, and other spectrometer configuration files.

    Args:
        path (str):     Path to either .d01 or .exp file

    Returns:
        attrs (dict):   Dictionary of parameter fields and values (DNPLab attributes)

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
    
    attrs = analyze_attrs(attrs)
    return attrs


def load_specman_d01(path, attrs, verbose=False):
    """Import SpecMan d01 data file

    DNPLab function to import the SpecMan d01 data file. The format of the SpecMan data file is described here:


    Args:
        path (str):         Path to either .d01 or .exp file

    Returns:
        data (ndarray):     SpecMan data as numpy array
        params (dict):      Dictionary with import updated parameters dictionary
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

    # Number of dimensions of stored data. Note: Not clear why this is repeated for every dimension. This seems to be the same number for all dimensions, but is repeated.
    attrs["dims"] = uint_read[2]

    dataShape = uint_read[2 : 2 + uint_read[0] * 6]

    attrs["dataStreamShape"] = dataShape

    attrs["dataStartIndex"] = 2 + attrs["numberOfVariables"] * 6

    if verbose == True:
        print("** Data paramters **")
        print("numberOfVariables : ", attrs["numberOfVariables"])
        print("dataFormat        : ", attrs["dataFormat"])
        print("dims              : ", attrs["dims"])
        print("dataStreamShape   : ", attrs["dataStreamShape"])
        print("dataStartIndex    : ", attrs["dataStartIndex"])

    data = float_read[attrs["dataStartIndex"] :]

    if attrs["dims"] == 1:
        data = _np.reshape(data, (uint_read[0], uint_read[3]), order="C")

    elif attrs["dims"] == 2:
        data = _np.reshape(data, (uint_read[0], uint_read[4], uint_read[3]), order="C")

    elif attrs["dims"] == 3:
        data = _np.reshape(
            data, (uint_read[0], uint_read[5], uint_read[4], uint_read[3]), order="C"
        )

    elif attrs["dims"] == 4:
        data = _np.reshape(
            data,
            (uint_read[0], uint_read[4], uint_read[5], uint_read[6], uint_read[3]),
            order="C",
        )

    elif attrs["dims"] >= 4:
        print("Maximum dimensionality for SpecMan data is 4D")
        return None

    # Swap first axis with last
    data = _np.swapaxes(data, 0, -1)

    # SpecMan data can have a maximum of four dimensions
    dims_full = ["x0", "x1", "x2", "x3", "x4"]
    dims = dims_full[0 : dataShape[0] + 1]

    dims, coords = specman_coords(attrs, data, dims)
    
    return data, dims, coords, attrs

def analyze_attrs(attrs):
    """
    Analyze the attrs and add some important attrs to existing dictionary

    Args:
        attrs (dict): Dictionary of specman acqusition parameters
    
    Returns:
        attrs (dict): The dictionary of specman acqusition parameters and added parameters

    """
    
    temp = {}
    for key, val in attrs.items():
        if 'params_' in key:
            new_key = key.split('params_')[1] # get key value for temp dictionary
            val = val.split(';')[0] # remove non value related information
            val_list = val.split(' ') # split value string for further analyze
            val = val_list[0]
            temp[new_key] = int(val) if '.' not in val else float(val)
            if 'step' in val_list: # when it indicate the step
                step_index = val_list.index('step') + 1 # the index of the value of 'step' is equal to the index of string 'index' + 1
                step = float(val_list[step_index])
                temp[new_key + '_step'] = step

        if 'sweep_' in key:
            val_list = val.split(',')
            val = val_list[1] # get value
            new_key = 'sweep_' + val_list[0]
            temp[new_key + '_length'] = int(val)
            # new_key += '_dim' # last item is the key to the parameters, such as t, p...
            temp[new_key + '_dim'] = val_list[-1]
        

    attrs = {**attrs, **temp}
    return attrs

def specman_coords(attrs, data: _dnp.DNPData, dims):
    """Generate coords from specman acquisition parameters

    Args:
        attrs (dict): Dictionary of specman acqusition parameters
        data (DNPData Object): The specman data
        dims (list): the list of current dims

    Returns:
        tuple: dims and coords
    """
    shape = _np.shape(data)
    length = attrs['sweep_X_length']
    new_dims = [attrs['sweep_X_dim']] # new dim name list
    coords = []
    for index in range(data.ndim):
        if index < len(new_dims):
            dim = new_dims[index]
            dims[index] = dim
            start = attrs[dim]
            step = attrs[dim + '_step']
            coord = []
            for i in range(length):
                coord.append(start + step * i)
        else:
            coord = _np.arange(0, shape[index])

        coords.append(_np.array(coord))
    return dims, coords


