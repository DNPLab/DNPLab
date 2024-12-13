import numpy as _np
from struct import unpack
from .. import DNPData
from matplotlib.pyplot import *

DELTA_DATA_FORMAT_DICT = {
    1: ["One_D", 8, 8**1],
    2: ["Two_D", 32, 32**2],
    3: ["Three_D", 8, 8**3],
    4: ["Four_D", 8, 8**4],
    5: ["Five_D", 4, 4**5],
    6: ["Six_D", 4, 4**6],
    7: ["Seven_D", 2, 2**7],
    8: ["Eight_D", 2, 2**8],
    12: ["Small_Two_D", 4, 4**2],
    13: ["Small_Three_D", 4, 4**3],
    14: ["Small_Four_D", 4, 4**4],
}

DELTA_DATA_FIELD_DICT = {
    # attrs : [dtype, offset, size, step, index]
    "File_Identifier": [str, 0, 8, 1, None],
    "Endian": [">B", 8, 1, 1, 0],
    "Major_Version": [">B", 9, 1, 1, 0],
    "Minor_Version": [">B", 10, 2, 1, 0],
    "Data_Dimension_Number": [">B", 12, 1, 1, 0],
    "Data_Dimension_Exist": [">B", 13, 1, 1, 0],
    "Data_Format": [">B", 14, 1, 1, 0],
    "Instrument": [">B", 15, 1, 1, 0],
    "Translate": [">B", 16, 8, 1, 0],
    "Data_Axis_Type": [">B", 24, 8, 1, None],
    "Data_Units": [">B", 32, 16, 2, None],
    "Title": [str, 48, 124, 1, None],
    "Data_Axis_Ranged": [">B", 172, 4, 1, None],
    "Data_Points": [">I", 176, 32, 4, None],
    "Data_Offset_Start": [">I", 208, 32, 4, None],
    "Data_Offset_Stop": [">I", 240, 32, 4, None],
    "Data_Axis_Start": [">d", 272, 64, 8, None],
    "Data_Axis_Stop": [">d", 336, 64, 8, None],
    "Node_Name": [str, 408, 16, 1, None],
    "Site": [str, 424, 128, 1, None],
    "Author": [str, 552, 128, 1, None],
    "Comment": [str, 608, 128, 1, None],
    "Data_Axis_Title": [str, 808, 256, 32, None],
    "Base_Freq": [">d", 1064, 64, 8, None],
    "Zero_Points": [">d", 1128, 64, 8, None],
    "Reversed": [">B", 1192, 8, 1, None],
    "History_Used": [">I", 1204, 4, 4, 0],
    "History_Length": [">I", 1208, 4, 4, 0],
    "Param_Start": [">I", 1212, 4, 4, 0],
    "Param_Length": [">I", 1216, 4, 4, 0],
    "List_Start": [">I", 1220, 32, 4, None],
    "List_Length": [">I", 1252, 32, 4, None],
    "Data_Start": [">I", 1284, 4, 4, 0],
    "Data_Length": [">I", 1288, 8, 4, 1],
    "Context_Start": [">I", 1296, 8, 4, 1],
    "Context_Length": [">I", 1304, 4, 4, 0],
    "Annote_Start": [">I", 1308, 8, 4, 1],
    "Annote_Length": [">I", 1316, 4, 4, 0],
    "Total_Size": [">I", 1320, 8, 4, 1],
    "Unit_Location": [">I", 1328, 8, 4, None],
}


# Data_Start = _np.array([
#     unpack(">I", file_contents[1284 + ix : 1288 + ix])[0] for ix in range(0, 4, 4)
# ][0])
# Data_Length = int.from_bytes(file_contents[1288:1296], byteorder="big")
def import_delta(path, verbose=False):
    """Import Delta data and return DNPData object

    Currently only 1D and 2D data sets are supported.

    Args:
        path (str)          : Path to .jdf file

    Returns:
        dnpdata (DNPData)   : DNPData object containing Delta data
    """

    # params = import_delta_pars(path)
    values, dims, coords, attrs = import_delta_data(path, verbose=verbose)

    out = DNPData(values, dims, coords, attrs)

    return out


def import_delta_pars(path, context_start):
    """Import parameter fields of Delta data

    Args:
        path (str) : Path to .jdf file
        context_start (int): the index where the context starts

    Returns:
        params (dict) : dictionary of parameter fields and values
    """
    file_opened = open(path, "rb")
    file_opened.seek(context_start)
    lines = file_opened.readlines()
    params = {}
    in_when_condition = False  # detecting when condition
    in_if_condition = False  # detecting if condition
    if_params_key = None
    for line in lines:
        try:
            line = (
                str(line.decode("utf-8"))
                .replace("\x00", "")
                .replace(" ", "")
                .replace(";", "")
                .replace("\n", "")
                .replace('"', "")
            )
            if not line:
                continue  # ignore empty line

            if "endwhen" in line:  # detect the end of when loop and disable detecting
                in_when_condition = False
                when_condition_dict = {}

            if (
                in_if_condition and if_params_key and "if" not in line
            ):  # when it is not the if condition in if condition
                if "else" in line:  # end of if condition
                    in_if_condition = False
                line = (
                    line.replace("then", "")
                    .replace("else", "")
                    .replace("[", "")
                    .replace("]", "")
                )
                for condition_key, acceptance in if_condition_dict.items():
                    if condition_key in params and params[condition_key] in acceptance:
                        in_if_condition = False  # end if condition
                        break

                params[if_params_key] = line
                if_params_key = None
            if "=" in line:
                if "when" in line:  # is a when condition
                    when_condition_dict = {}
                    conditions = line.split("or")
                    for condition in conditions:
                        condition_key, condition_val = (
                            condition.replace("when", "").replace("do", "").split("=")
                        )
                        if condition_val == "TRUE":
                            condition_val = True
                        elif condition_val == "False":
                            condition_val = False
                        if condition_key not in when_condition_dict:
                            when_condition_dict[condition_key] = []
                        when_condition_dict[condition_key].append(condition_val)
                    in_when_condition = True

                elif "if" in line:  # is a if condition
                    if "else" not in line:
                        if_condition_dict = {}
                        elements = line.replace("if", "").split("=")
                        if len(elements) == 3:
                            if_params_key, condition_key, condition_val = elements
                        elif len(elements) == 2:
                            if_params_key, conditions = elements
                            if "not" in conditions:
                                condition_val = False
                            else:
                                condition_val = True
                            condition_key = conditions.replace("not", "")
                    else:
                        condition_key, condition_val = line.replace("if", "").split("=")
                    if condition_val == "TRUE":
                        condition_val = True
                    elif condition_val == "False":
                        condition_val = False
                    if condition_key not in if_condition_dict:
                        if_condition_dict[condition_key] = []
                    if_condition_dict[condition_key].append(condition_val)
                    in_if_condition = True

                else:
                    if "help" in line:
                        line = line[
                            : line.find("help") - 1
                        ]  # remove help information with the comma
                    info = line.split("=")
                    key = info[0]
                    val = info[1]
                    if val[0] == ">":
                        val = val[1:]
                    if val == "TRUE":  # Boolean True
                        val = True
                    elif val == "FALSE":  # Boolean True
                        val = False
                    elif val.isdigit():  # just number
                        val = int(val)
                    elif (
                        "[" in val
                        and "]" in val
                        and "[" not in val[val.find("[") + 1 :]
                        and "?" not in val
                    ):  # only a number and a unit
                        val, unit = (
                            val.replace("[", ",").replace("]", "").split(",")
                        )  # separate value and unit
                        params[key + "_unit"] = unit
                        val = float(val)
                    elif val[0] == "?" or "round" in val:  # math operation
                        val = val.replace("?", "")
                        try:
                            if val[0] not in ["#", "*"]:
                                # get all variables
                                vars = (
                                    val.replace("round", "")
                                    .replace("(", " ")
                                    .replace(")", " ")
                                    .replace("+", " ")
                                    .replace("-", " ")
                                    .replace("*", " ")
                                    .replace("/", " ")
                                    .split()
                                )
                                temp_val = val
                                for var in vars:
                                    if (
                                        not var.isdigit()
                                    ):  # variable is not a number, then it is a key in dictionary
                                        temp_val = temp_val.replace(
                                            var, 'params["%s"]' % var
                                        )
                                temp_val = temp_val.replace("/", "*1/").replace(
                                    "-", "+ -1 * "
                                )  # for safty, remove '-' and '/' from string
                                val = eval("%s" % temp_val)
                        except:
                            pass

                    elif "{" in val and "}" in val:  # is a list
                        val = val.replace("{", "[").replace("}", "]")
                        # list_index = val.find("[")

                        # if (
                        #     list_index == 0
                        # ):  # when there are some information before this val
                        #     try:
                        #         if "(" not in val and ")" not in val:
                        #             val = eval("%s" % val)

                        #     except:
                        #         pass

                    if in_when_condition:
                        for condition_key, acceptance in when_condition_dict.items():
                            if (
                                condition_key in params
                                and params[condition_key].split(",")[0] in acceptance
                            ):
                                params[key] = val
                                break
                    else:
                        params[key] = val

        except UnicodeDecodeError:
            pass
    return params


def import_delta_data(path, params={}, verbose=False):
    """Import spectrum or spectra of Delta data

    Currently only 1D and 2D data sets are supported.

    Args:
        path (str) : Path to .jdf file
        params (dict) : dictionary of parameters

    Returns:
        y_data (ndarray) : spectrum or spectra if >1D
        abscissa (list) : coordinates of axes
        dims (list) : axes names
        params (dict) : updated dictionary of parameters
    """

    file_opened = open(path, "rb")
    file_contents = file_opened.read()
    file_opened.close()
    params = params

    for key, val in DELTA_DATA_FIELD_DICT.items():
        dtype, offset, size, step, index = val
        if dtype == str:
            if step == 1:
                params[key] = str(
                    file_contents[offset : offset + size].decode("utf-8")
                ).replace("\x00", "")
            else:
                res = [
                    str(
                        file_contents[offset + ix : offset + ix + step].decode("utf-8")
                    ).replace("\x00", "")
                    for ix in range(0, size, step)
                ]

        else:
            if "B" in dtype:
                bits = 1
            elif "I" in dtype:
                bits = 4
            elif "d" in dtype:
                bits = 8

            res = _np.array(
                [
                    unpack(val[0], file_contents[offset + ix : offset + bits + ix])[0]
                    for ix in range(0, size, step)
                ]
            )
            if index != None:
                res = res[index]
            params[key] = res

    Endian = "<d" if params["Endian"] else ">d"
    Data_Dimension_Number = params["Data_Dimension_Number"]
    Data_Format = params["Data_Format"]

    # Data_Axis_Type
    # Array of 8 enumerations. Each element indicates the type of data for that axis.
    # Together with Data_Format this values determines the data layout in the data
    # section. Values are: 0 = None, axis not used, 1 = Real, Axis has real data only
    # 2 = TPPI, 3 = Complex, Axis has complex data, 4 = Real_Complex, Axis should be
    # accessed as complex when it is the major axis, accessed as real otherwise. This
    # is only valid when all axes in use have this setting, 5 = Envelop.
    Data_Axis_Type = params["Data_Axis_Type"]
    Data_Units = params["Data_Units"]

    params["units"] = []
    for ix in range(Data_Dimension_Number):
        if Data_Units[ix] == 1:
            params["units"].append("abundance")
        elif Data_Units[ix] == 13:
            params["units"].append("Hz")
        elif Data_Units[ix] == 26:
            params["units"].append("ppm")
        elif Data_Units[ix] == 27:
            params["units"].append("rad")
        elif Data_Units[ix] == 28:
            params["units"].append("s")
        else:
            params["units"].append("indexed")

    Data_Points = params["Data_Points"]
    Data_Offset_Start = params["Data_Offset_Start"]
    Data_Offset_Stop = params["Data_Offset_Stop"]
    Valid_pts = Data_Offset_Stop - Data_Offset_Start + 1
    Data_Axis_Start = params["Data_Axis_Start"]
    Data_Axis_Stop = params["Data_Axis_Stop"]
    Base_Freq = params["Base_Freq"]
    params["nmr_frequency"] = _np.array(Base_Freq[0]) * 1e6
    Data_Start = params["Data_Start"]
    Data_Length = params["Data_Length"]
    Total_Size = params["Total_Size"]
    context = import_delta_pars(path, params["Context_Start"])
    params = {**params, **context}
    abscissa = []

    # analyze the coords
    pre_defined_axis = ["x", "y", "z", "t"]
    for k in range(Data_Dimension_Number):
        axis = pre_defined_axis[k]
        interval = params.get("interval")
        if interval and axis in interval:
            val = params["interval"].replace("%s_acq" % axis, "")
            val = val.replace("[ks]", "e3")
            val = val.replace("[s]", "")
            val = val.replace("[ms]", "e-3")
            val = val.replace("[us]", "e-6")
            val = val.replace("[ns]", "e-9")
            val = val.replace("[ps]", "e-12")
            abscissa.append(np.array(eval(val)))
        else:
            abscissa.append(
                _np.linspace(Data_Axis_Start[k], Data_Axis_Stop[k], Valid_pts[k])
            )

    file_opened = open(path, "rb")
    file_opened.seek(Data_Start)

    # Read data from file
    if all(
        Data_Axis_Type[:Data_Dimension_Number] == 4
    ):  # Only valid if all axis types are 4
        sections = (
            2  # Axis is accessed as complex when it is major axis and others are real.
        )
    else:
        sections = 2 ** sum(
            Data_Axis_Type > 2
        )  # if an axis type is complex (type 3 or 4), then there are two sections for this axis
    read_pts = (
        _np.prod(Data_Points) * sections
    )  # if Data_Points is in complex, the data is saved in the Data_Points * number of sections

    if sections > 2:  # then the data is in hyper complex format
        raise TypeError("Hyper Complex Format is detected and not supported.")

    data = _np.fromfile(file_opened, Endian, read_pts)
    file_opened.close()

    data_format, submatrix_edge, submatrix_points = DELTA_DATA_FORMAT_DICT[Data_Format]
    # 1D data reshaping
    if Data_Dimension_Number == 1:
        if Data_Axis_Type[0] == 1:
            temp = data

        elif Data_Axis_Type[0] == 3 or Data_Axis_Type[0] == 4:
            temp = _np.split(data, 2)[0] - 1j * _np.split(data, 2)[1]

        else:
            raise TypeError("Data format not recognized")

        # Extract the actual data by using Data_Offset_Start and Data_Offset_Stop indexes
        out = temp[Data_Offset_Start[0] : Data_Offset_Stop[0] + 1]
        dims = ["t2"]

    # 2D data reshaping
    elif Data_Dimension_Number == 2:
        """
        Data is saved as the order of submatrices.
        E.g, assume a 2D dataset has 4 spectra, each spectrum has 4 data points (4*4), then:

            16 total submatrices laid out 2*2, each submatrix is 2*2

            m1 = |1 2|, m2 = |5 6|, m3 = |9  10 |, m4 = |13 14|
                 |3 4|       |7 8|       |11 12 |       |15 16|

            M = | 1   2   5   6  |
                | 3   4   7   8  |
                | 9   10  11  12 |
                | 13  14  15  16 |

        The data are save in the file:
            data = [1, 2, 3, 4, 5, 6, ....]

        But the spectrum data is:
            s1 = [1, 2, 5, 6]
            s2 = [3, 4, 7, 8]
            s3 = [9 ,10 ,11, 12]
            s4 = [13, 14, 15, 16]

        Each M is called a 'section' in JEOL dataset. If there are two sections, e.g. 2D complex data,
        the whole dataset must be separated into 2 sections, with the 1st section is for real, and 2nd is for image.

        The information can be found in JEOL documentation.

        ** Please be aware of the dataset layout in matrix. For m1 in 2d numpy.array, it is [[1,3],
                                                                                             [2,4]]

        """
        if all(Data_Axis_Type != 2):  # TPPI axis format is not supported
            # Step 1, separate real and image sections: first section is real and second section is image
            data_folded = (
                _np.split(data, sections)[0] - 1j * _np.split(data, sections)[1]
            )

            # Step 2: reshape to the layout of submatrices, shape = (matrix_x, matrix_y, submatrix_edge, submatrix_edge)
            # maxtrix_x is the number of submatrice in row and matrix_y is the number of submatrice in column
            # at this point, first two and second two axes are swapped
            temp = _np.reshape(
                data_folded,
                (
                    Data_Points[1] // submatrix_edge,
                    Data_Points[0] // submatrix_edge,
                    submatrix_edge,
                    submatrix_edge,
                ),
            )

            # Step 3: swap axes
            ndims = temp.ndim
            for dim in range(ndims - 1, 0, -2):
                temp = _np.swapaxes(temp, dim, dim - 1)

            # Step 4: stack data horizontally twice to get full matrix
            temp = _np.hstack(_np.hstack(temp))

            # Step 5: select data
            temp1 = temp[
                Data_Offset_Start[0] : Data_Offset_Stop[0] + 1,
                Data_Offset_Start[1] : Data_Offset_Stop[1] + 1,
            ]

            out = temp1

        else:
            raise ValueError("Data format not recognized")

        dims = ["t2", "t1"]

    else:
        raise TypeError("Only 1D or 2D are supported")

    if verbose == True:
        print("Endian: ", Endian)
        print("Data Dimension Number: ", Data_Dimension_Number)
        print("Data Format: ", data_format)
        print("Data Axis Type: ", Data_Axis_Type)
        print("Data Units: ", Data_Units)
        print("Data Points: ", Data_Points)
        print("Data Offset Start: ", Data_Offset_Start)
        print("Data Offset Stop: ", Data_Offset_Stop)
        print("Valid points: ", Valid_pts)
        print("Data Axis Start", Data_Axis_Start)
        print("data Axis Stop: ", Data_Axis_Stop)
        print("Base Freq: ", Base_Freq)
        print("Data Start: ", Data_Start)
        print("Data Length: ", Data_Length)
        print("Submatrix Points: ", submatrix_points)
        print("Submatrix Edge: ", submatrix_edge)
        print("Total Size: ", Total_Size)
        print("Output data shape: ", _np.shape(out))

    return out, dims, abscissa, params
