from ..core.data import DNPData
import numpy as _np
from struct import unpack
import warnings
import os
import glob


def import_prospa(path, parameters_filename=None, experiment=None, verbose=False):
    """Import Kea data

    Args:
        path (str): Path to data
        parameters_filename (str):
        experiment (str): Prospa experiment, used when calculating coords from parameters
        verbose (bool): If true, prints additional information for troubleshooting

    Returns:
        dnpdata object with Kea data
    """

    if parameters_filename == None:
        parameters_filename = "acqu.par"
    extension = ""

    if os.path.isfile(path):
        path, filename = os.path.split(path)
        filename, extension = os.path.splitext(filename)
    elif os.path.isdir(path):
        filesList = glob.glob(os.path.join(path, "*.[1-4]d"))
        if len(filesList) == 0:
            raise ValueError("No binary data file in directory:", path)
        elif len(filesList) > 1:
            raise ValueError("More than one binary data file in directory:", filesList)
        else:
            data_filename = filesList[0]
            filename, extension = os.path.splitext(os.path.split(data_filename)[-1])
    else:
        raise OSError("%s not found" % path)

    try:
        attrs = import_par(os.path.join(path, parameters_filename))
    except:
        warnings.warn("No parameters file in directory")
        attrs = {}

    if extension == ".csv":
        # Import csv data
        x, data = import_csv(os.path.join(path, filename + extension))
    else:
        # Import Binary data
        x, data = import_nd(os.path.join(path, filename + extension))

    if "b1Freq" in attrs:
        nmr_frequency = attrs["b1Freq"]
        try:
            nmr_frequency = float(nmr_frequency)
        except:
            nmr_frequency = float(nmr_frequency.replace("d", ""))

        attrs["nmr_frequency"] = nmr_frequency * 1e6
        attrs["spectrometer_format"] = "prospa"

    # Assume direct dimension is 1st dimension
    data_shape = _np.shape(_np.squeeze(data))

    dims, coords = prospa_coords(attrs, data_shape, experiment=experiment)

    kea_data = DNPData(data, dims, coords, attrs)

    return kea_data


# Depreciated
def import_prospa_dir(path, exp_list=None):
    """Import directory of prospa experiments

    Args:
        path (str): Directory of all data
        exp_list (list): List of directorys to include

    Returns:
        dict: Dictionary of DNPData objects
    """

    dirs = [x for x in os.listdir(path) if os.path.isdir(os.path.join(path, x))]

    if exp_list is not None:
        dirs = [dir_ for dir_ in dirs if dir_ in exp_list]

    ws = {}

    for ix, dir_ in enumerate(dirs):
        try:
            tmp = import_prospa(os.path.join(path, dir_))
            ws[dir] = tmp
        except ValueError as e:
            print("Skipping folder: %s" % str(e))

    return ws


def import_nd(path):
    """Import Kea binary 1d, 2d, 3d, 4d files

    Args:
        path (str): Path to file

    Returns:
        tuple:
            x (None, numpy.array): Axes if included in binary file, None otherwise
            data (numpy.array): Numpy array of data
    """

    ascii_header_size = 12
    ascii_header_format = "<4s4s4s"

    with open(path, "rb") as f:
        ascii_header_bytes = f.read(ascii_header_size)

        ascii_header = unpack(ascii_header_format, ascii_header_bytes)

        owner = ascii_header[0][::-1].decode("utf-8")
        format_ = ascii_header[1][::-1].decode("utf-8")
        version = ascii_header[2][::-1].decode("utf-8")

        # header size depends on version
        if version == "V1.0":
            header_format = "4i"
            header_size = 16
        else:
            header_format = "5i"
            header_size = 20

        header_bytes = f.read(header_size)

        header = unpack(header_format, header_bytes)

        dataType = header[0]
        xDim = header[1]
        yDim = header[2]
        zDim = header[3]
        if version != "V1.0":
            qDim = header[4]
        else:
            qDim = 1

        raw = f.read()

        x = None
        if dataType == 500:  # float
            raw_data = unpack("<%if" % (xDim * yDim * zDim * qDim), raw)
            data = _np.array(raw_data)
        elif dataType == 501:  # complex
            raw_data = unpack("<%if" % (xDim * yDim * zDim * qDim * 2), raw)
            data = _np.array(raw_data)
            data = data[0::2] - 1j * data[1::2]
        elif dataType == 502:  # double
            raw_data = unpack("<%id" % (xDim * yDim * zDim * qDim), raw)
            data = _np.array(raw_data)
        elif dataType == 503:
            raw_data = unpack("<%if" % (xDim * yDim * zDim * qDim * 2), raw)
            raw_data = _np.array(raw_data)
            x = raw_data[0:xDim]
            data = raw_data[xDim:]
        elif dataType == 504:
            raw_data = unpack("<%if" % (xDim * yDim * zDim * qDim * 3), raw)  # 504
            raw_data = _np.array(raw_data)
            x = raw_data[0:xDim]
            data = raw_data[xDim:]
            data = data[0::2] - 1j * data[1::2]
        else:
            raise ValueError("Data %i type not recognized" % dataType)

        # re-shape using F-ordering (Fortran)
        data = data.reshape(xDim, yDim, zDim, qDim, order="F")

        data = data.squeeze()  # remove length 1 dimensions

    return x, data


def import_par(path):
    """Import Kea parameters .par file

    Args:
        path (str): Path to parameters file

    Returns:
        dict: Dictionary of Kea Parameters
    """

    attrs = {}

    with open(path, "r") as f:
        raw = f.read()
        lines = raw.rstrip().rsplit("\n")

        for line in lines:
            key, value = line.rstrip().rsplit(" = ")

            if value[0] == '"' and value[-1] == '"':
                value = value[1:-1]

            if "." in value:
                try:
                    value = float(value)
                except:
                    pass
            else:
                try:
                    value = int(value)
                except:
                    pass
            attrs[key] = value

    return attrs


def import_csv(path, return_raw=False, is_complex=True):
    """Import Kea csv file

    Args:
        path (str): Path to csv file

    Returns:
        tuple:
            x(numpy.array): axes if return_raw = False
            data(numpy.array): Data in csv file
    """

    raw = _np.loadtxt(path, delimiter=",")

    if not return_raw:
        x = raw[:, 0]
        if is_complex:
            data = raw[:, 1::2] + 1j * raw[:, 2::2]
        else:
            data = raw[:, 1:]
        data = _np.squeeze(data)
        return x, data
    else:
        return raw


def prospa_coords(attrs, data_shape, experiment):
    """Generate coords from prospa acquisition parameters

    Args:
        attrs (dict): Dictionary of prospa acqusition parameters
        data_shape (tuple): Shape of data

    Returns:
        tuple: dims and coords
    """

    if experiment is None:
        experiment = attrs["experiment"]
    dims = []
    coords = []

    if experiment == "1Pulse":
        pts = attrs["nrPnts"]
        dwell_time = attrs["dwellTime"]
        x = _np.arange(0.0, pts * dwell_time, dwell_time) / 1e6
        dims.append("t2")
        coords.append(x)

    elif experiment == "B12T_1Pulse" or experiment == "B12T_1Pulse_MPS":
        pts = attrs["nrPnts"]
        dwell_time = attrs["dwellTime"]
        x = _np.arange(0.0, pts * dwell_time, dwell_time) / 1e6
        dims.append("t2")
        coords.append(x)

        dims.append("Average")
        coords.append(_np.arange(attrs["nrScans"]))

    elif experiment == "T1-IR-FID":
        pts = attrs["nrPnts"]
        dwell_time = attrs["dwellTime"]
        x = _np.arange(0.0, pts * dwell_time, dwell_time) / 1e6
        dims.append("t2")
        coords.append(x)

        T1_steps = attrs["nrSteps"]
        T1_min_delay = attrs["minDelay"]
        T1_max_delay = attrs["maxDelay"]

        if attrs["delaySpacing"] == "lin":
            T1 = _np.linspace(T1_min_delay, T1_max_delay, T1_steps) / 1000.0
        elif attrs["delaySpacing"] == "log":
            T1 = (
                _np.logspace(_np.log10(T1_min_delay), _np.log10(T1_max_delay), T1_steps)
                / 1000.0
            )
        else:
            raise ValueError(
                f"Unable to determine delaySpacing {attrs['delaySpacing']}"
            )

        dims.append("t1")
        coords.append(T1)
    elif experiment == "B12T_T1-IR-FID" or experiment == "B12T_T1-IR-FID_MPS":
        pts = attrs["nrPnts"]
        dwell_time = attrs["dwellTime"]
        x = _np.arange(0.0, pts * dwell_time, dwell_time) / 1e6
        dims.append("t2")
        coords.append(x)

        T1_steps = attrs["nrSteps"]
        T1_min_delay = attrs["minDelay"]
        T1_max_delay = attrs["maxDelay"]

        if attrs["delaySpacing"] == "lin":
            T1 = _np.linspace(T1_min_delay, T1_max_delay, T1_steps) / 1000.0
        else:
            T1 = (
                _np.logspace(_np.log10(T1_min_delay), _np.log10(T1_max_delay), T1_steps)
                / 1000.0
            )

        dims.append("t1")
        coords.append(T1)

        dims.append("Average")
        coords.append(_np.arange(attrs["nrScans"]))
    elif experiment == "B12T_jres2D":
        pts = attrs["nrPnts"]
        dwell_time = attrs["dwellTime"]
        x = _np.arange(0.0, pts * dwell_time, dwell_time) / 1e6
        dims.append("t2")
        coords.append(x)

        inter_pulse_delay = attrs["interPulseDelay"]
        increment = attrs["increment"]
        steps = attrs["nrSteps"]

        t1 = (
            _np.r_[
                inter_pulse_delay : inter_pulse_delay
                + increment * (steps - 1) : 1j * steps
            ]
            / 1e6
        )

        dims.append("t1")
        coords.append(t1)

        dims.append("Average")
        coords.append(_np.arange(attrs["nrScans"]))

    else:
        dims_list = ["x", "y", "z", "q"]
        for ix in range(len(data_shape)):
            dims.append(dims_list[ix])  # call dimensions in order: x, y, z, q
            coords.append(_np.arange(data_shape[ix]))  # set coords to index for now
        if ("nrPnts" in attrs) and ("dwellTime" in attrs):
            pts = attrs["nrPnts"]
            dwell_time = attrs["dwellTime"]

            x = _np.arange(0.0, pts * dwell_time, dwell_time) / 1e6

            dims[0] = "t2"
            coords[0] = x

    return dims, coords
