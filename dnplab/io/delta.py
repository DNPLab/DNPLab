import numpy as _np
import re
from struct import unpack
from .. import DNPData


def import_delta(path):
    """Import Delta data and return DNPData object

    Currently only 1D and 2D data sets are supported.

    Args:
        path (str)          : Path to .jdf file

    Returns:
        dnpdata (DNPData)   : DNPData object containing Delta data
    """

    params = import_delta_pars(path)
    values, dims, coords, attrs = import_delta_data(path, params)

    out = DNPData(values, dims, coords, attrs)

    return out


def import_delta_pars(path):
    """Import parameter fields of Delta data

    Args:
        path (str) : Path to .jdf file

    Returns:
        params (dict) : dictionary of parameter fields and values
    """

    file_opened = open(path, "rb")
    file_contents = file_opened.readlines()
    file_opened.close()
    params = {}
    for ix in range(len(file_contents)):
        try:
            file_contents[ix] = str(file_contents[ix], "utf8")
            if "=" in file_contents[ix]:
                new_line = file_contents[ix].split("=")
                new_name = new_line[0].replace("/", "").strip()
                if "@" in new_name or new_name == "":
                    pass
                else:
                    params[new_name] = (
                        new_line[1]
                        .replace(">", "")
                        .replace(";", "")
                        .replace("?", "")
                        .strip()
                    )
        except UnicodeDecodeError:
            pass

    return params


def import_delta_data(path, params, verbose=True):
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

    if not params:
        params = {}

    file_opened = open(path, "rb")
    file_contents = file_opened.read(1360)
    file_opened.close()

    ##### Reading field header section #####
    Endian = [unpack(">B", file_contents[8 + ix : 9 + ix])[0] for ix in range(0, 1, 1)][
        0
    ]
    if Endian == 0:
        Endian = ">d"
    elif Endian == 1:
        Endian = "<d"
    else:
        raise UnicodeTranslateError("Failed to determine endianness")

    Data_Dimension_Number = [
        unpack(">B", file_contents[12 + ix : 13 + ix])[0] for ix in range(0, 1, 1)
    ][0]

    # Data_Type =  int.from_bytes(file_contents[14:16], byteorder="big")

    Data_Format = [
        unpack(">B", file_contents[14 + ix : 15 + ix])[0] for ix in range(0, 1, 1)
    ][0]

    # Data_Axis_Type = [
    #     unpack(">B", file_contents[24 + ix : 25 + ix])[0] for ix in range(0, 8, 1)
    # ][:Data_Dimension_Number]

    # Data_Axis_Type
    # Array of 8 enumerations. Each element indicates the type of data for that axis.
    # Together with Data_Format this values determines the data layout in the data
    # section. Values are: 0 = None, axis not used, 1 = Real, Axis has real data only
    # 2 = TPPI, 3 = Complex, Axis has complex data, 4 = Real_Complex, Axis should be
    # accessed as complex when it is the major axis, accessed as real otherwise. This
    # is only valid when all axes in use have this setting, 5 = Envelop.
    Data_Axis_Type = [
        unpack(">B", file_contents[24 + ix : 25 + ix])[0] for ix in range(0, 8, 1)
    ]

    # Data_Units = [
    #     unpack(">B", file_contents[32 + ix : 33 + ix])[0] for ix in range(1, 16, 2)
    # ][:Data_Dimension_Number]

    Data_Units = [
        unpack(">B", file_contents[32 + ix : 33 + ix])[0] for ix in range(1, 16, 2)
    ]

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



    #### We need to decide whether we want these parameters as a list or an array. I would
    # suggest array to be able to do math.


    Data_Points = _np.array([
        unpack(">I", file_contents[176 + ix : 180 + ix])[0] for ix in range(0, 32, 4)
    ])

    Data_Offset_Start = _np.array([
        unpack(">I", file_contents[208 + ix : 212 + ix])[0] for ix in range(0, 32, 4)
    ])

    Data_Offset_Stop = _np.array([
        unpack(">I", file_contents[240 + ix : 244 + ix])[0] for ix in range(0, 32, 4)
    ])

    Valid_pts = Data_Offset_Stop - Data_Offset_Start + 1





    Data_Axis_Start = [
        unpack(">d", file_contents[272 + ix : 280 + ix])[0] for ix in range(0, 64, 8)
    ]

    Data_Axis_Stop = [
        unpack(">d", file_contents[336 + ix : 344 + ix])[0] for ix in range(0, 64, 8)
    ]

    Base_Freq = [
        unpack(">d", file_contents[1064 + ix : 1072 + ix])[0] for ix in range(0, 64, 8)
    ]

    params["nmr_frequency"] = [
        unpack(">d", file_contents[1064 + ix : 1072 + ix])[0] for ix in range(0, 64, 8)
    ][
        0
    ] * 1e6  # convert from MHz to Hz, is this always in Hz?

    Data_Start = [
        unpack(">I", file_contents[1284 + ix : 1288 + ix])[0] for ix in range(0, 4, 4)
    ][0]

    Data_Length = int.from_bytes(file_contents[1288:1296], byteorder="big")

    Total_Size = int.from_bytes(file_contents[1320:1328], byteorder="big")

    # Number_sections = 2**

    abscissa = []

    for k in range(Data_Dimension_Number):
        abscissa.append(
            # _np.linspace(Data_Axis_Start[k], Data_Axis_Stop[k], Data_Points[k])
            _np.linspace(Data_Axis_Start[k], Data_Axis_Stop[k], Valid_pts[k])
        )

    file_opened = open(path, "rb")
    file_opened.seek(Data_Start)




    # Read data from file
    if Data_Dimension_Number == 1:
        read_pts = _np.prod(Data_Points) * 2

    elif Data_Dimension_Number == 2:
        if Data_Axis_Type[0] == 3 and Data_Axis_Type[1] == 1:
            read_pts = _np.prod(Data_Points) * 2

        elif Data_Axis_Type[0] == 3 and Data_Axis_Type[1] == 3:
            read_pts = _np.prod(Data_Points) * 4


    # if Data_Dimension_Number == 2 and Data_Axis_Type[0] == 3 and Data_Axis_Type[1] == 3:
    #     read_pts = _np.prod(Data_Points) * 4
    # else:
    #     read_pts = _np.prod(Data_Points) * 2

    data = _np.fromfile(file_opened, Endian, read_pts)
    file_opened.close()




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
        if Data_Axis_Type[0] == 3 and Data_Axis_Type[1] == 1:


        # if Data_Axis_Type[0] == 4 or (
        #     Data_Axis_Type[0] == 3 and Data_Axis_Type[1] == 1
        # ):

            data_folded = _np.split(data, 2)[0] - 1j * _np.split(data, 2)[1]

            print("data: ", _np.shape(data))
            print("data_folded: ", _np.shape(data_folded))


            # data_shaped = _np.reshape(
            #     data_folded,
            #     [int(Data_Points[0] / 4), int(Data_Points[1] / 4), 4, 4],
            #     order="C",
            # )

            
            temp = _np.reshape(data_folded, (Data_Points[0], Data_Points[1]))
            temp1 = temp[Data_Offset_Start[0]:Data_Offset_Stop[0]+1, Data_Offset_Start[1]:Data_Offset_Stop[1]+1]

            out = temp1

            # print("Shape temp:", _np.shape(temp))
            # out = _np.concatenate(_np.concatenate(data_shaped, 1), 1)




        # elif Data_Axis_Type[0] == 3 and Data_Axis_Type[1] == 3:
        #     print("Complex")
        #     data_folded = [
        #         _np.split(data, 4)[0] - 1j * _np.split(data, 4)[1],
        #         _np.split(data, 4)[2] - 1j * _np.split(data, 4)[3],
        #     ]
        #     for idx in enumerate(data_folded):
        #         data_shaped[idx] = _np.reshape(
        #             data_folded[idx],
        #             [int(Data_Points[0] / 32), int(Data_Points[1] / 32), 32, 32],
        #             order="F",
        #         )
        #         out[idx] = _np.concatenate(_np.concatenate(data_shaped[idx], 1), 1)

        else:
            raise ValueError("Data format not recognized")

        dims = ["t2", "t1"]

    else:
        raise TypeError("Only 1D or 2D are supported")
    
    if verbose == True:
        print("Endian: ", Endian)
        print("Data Dimension Number: ", Data_Dimension_Number)
        print("Data Format: ", Data_Format)
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
        print("Total Size: ", Total_Size)
        print("Output data shape: ", _np.shape(out))

    return out, dims, abscissa, params
