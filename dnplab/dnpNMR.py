from . import dnpdata, dnpdata_collection
import numpy as _np


def return_data(all_data):

    is_workspace = False
    if isinstance(all_data, dnpdata):
        data = all_data.copy()
    elif isinstance(all_data, dict):
        raise ValueError("Type dict is not supported")
    elif isinstance(all_data, dnpdata_collection):
        is_workspace = True
        if all_data.processing_buffer in all_data.keys():
            data = all_data[all_data.processing_buffer]
        else:
            raise ValueError("No data in processing buffer")
    else:
        raise ValueError("Data type not supported")

    return data, is_workspace


def update_parameters(proc_parameters, requiredList, default_parameters):
    """Add default parameter to processing parameters if a processing parameter is missing

    Args:
        proc_parameters (dict): Dictionary of initial processing parameters
        requiredList (list): List of requrired processing parameters
        default_parameters (dict): Dictionary of default processing parameters

    Returns:
        dict: Updated processing parameters dictionary
    """
    updatedProc_parameters = proc_parameters
    for requiredParameter in requiredList:
        if not requiredParameter in updatedProc_parameters:
            updatedProc_parameters[requiredParameter] = default_parameters[
                requiredParameter
            ]

    return updatedProc_parameters


def remove_offset(all_data, dim="t2", offset_points=10):
    """Remove DC offset from FID by averaging the last few data points and subtracting the average

    Args:
        all_data (dnpdata, dict): Data container for data
        proc_parameters (dict,procParam): Processing _parameters

    +---------------+------+---------+----------------------------------------------------------+
    | parameter     | type | default | description                                              |
    +---------------+------+---------+----------------------------------------------------------+
    | dim           | str  | 't2'    | Dimension to calculate DC offset                         |
    +---------------+------+---------+----------------------------------------------------------+
    | offset_points | int  | 10      | Number of points at end of data to average for DC offset |
    +---------------+------+---------+----------------------------------------------------------+

    Returns:
        dnpdata_collection: If workspace is given returns dnpdata_collection with data in processing buffer updated
        dnpdata: If dnpdata object is given, return dnpdata object.

    Example::

       workspace = dnplab.dnpNMR.remove_offset(workspace)
    """

    # Determine if data is dictionary or dnpdata object
    data, isDict = return_data(all_data)

    proc_parameters = {
        "dim": dim,
        "offset_points": offset_points,
    }

    dim = proc_parameters["dim"]
    offset_points = int(proc_parameters["offset_points"])

    offsetData = data["t2", -1 * offset_points :].values
    offsetData = offsetData.reshape(-1)
    offset = _np.mean(offsetData)

    data -= offset

    proc_attr_name = "remove_offset"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data


def fourier_transform(
    all_data, dim="t2", zero_fill_factor=2, shift=True, convert_to_ppm=True
):
    """Perform Fourier Transform down dim dimension given in proc_parameters

    .. Note::
        Assumes dt = t[1] - t[0]

    Args:
        all_data (dnpdata, dict): Data container

    +------------------+------+---------+--------------------------------------------------+
    | parameter        | type | default | description                                      |
    +------------------+------+---------+--------------------------------------------------+
    | dim              | str  | 't2'    | dimension to Fourier transform                   |
    +------------------+------+---------+--------------------------------------------------+
    | zero_fill_factor | int  | 2       | factor to increase dim with zeros                |
    +------------------+------+---------+--------------------------------------------------+
    | shift            | bool | True    | Perform fftshift to set zero frequency to center |
    +------------------+------+---------+--------------------------------------------------+
    | convert_to_ppm   | bool | True    | Convert dim from Hz to ppm                       |
    +------------------+------+---------+--------------------------------------------------+

    Returns:
        all_data (dnpdata, dict): Processed data in container

    Example:

    .. code-block:: python

        all_data = dnplab.dnpNMR.fourier_transform(all_data, proc_parameters)
    """

    # Determine if data is dictionary or dnpdata object
    data, isDict = return_data(all_data)

    # handle zero_fill_factor
    zero_fill_factor = int(zero_fill_factor)
    if zero_fill_factor <= 0:
        zero_fill_factor = 1

    proc_parameters = {
        "dim": dim,
        "zero_fill_factor": zero_fill_factor,
        "shift": shift,
        "convert_to_ppm": convert_to_ppm,
    }

    index = data.dims.index(dim)
    dt = data.coords[dim][1] - data.coords[dim][0]
    n_pts = zero_fill_factor * len(data.coords[dim])
    f = (1.0 / (n_pts * dt)) * _np.r_[0:n_pts]
    if shift == True:
        f -= 1.0 / (2 * dt)

    if convert_to_ppm:
        nmr_frequency = data.attrs["nmr_frequency"]
        f /= -1 * nmr_frequency / 1.0e6

    data.values = _np.fft.fft(data.values, n=n_pts, axis=index)
    if shift:
        data.values = _np.fft.fftshift(data.values, axes=index)
    data.coords[dim] = f

    proc_attr_name = "fourier_transform"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data


def window(all_data, dim="t2", linewidth=10):
    """Apply Apodization to data down given dimension

    Args:
        all_data (dnpdata, dict): data container
        proc_parameters (dict, procParam): parameter values

    .. note::
        Axis units assumed to be seconds

    +-----------+-------+---------+--------------------------------------------+
    | parameter | type  | default | description                                |
    +-----------+-------+---------+--------------------------------------------+
    | dim       | str   | 't2'    | Dimension to apply exponential apodization |
    +-----------+-------+---------+--------------------------------------------+
    | linewidth | float | 10      | Linewidth of broadening to apply in Hz     |
    +-----------+-------+---------+--------------------------------------------+

    Returns:
        dnpdata_collection or dnpdata: data object with window function applied

    Example:

    .. code-block:: python

        proc_parameters = {
                'linewidth' : 10,
                'dim' : 't2',
                }
        all_data = dnplab.dnpNMR.window(all_data,proc_parameters)

    """

    data, isDict = return_data(all_data)
    proc_parameters = {"dim": dim, "linewidth": linewidth}

    index = data.dims.index(dim)

    reshape_size = [1 for k in data.dims]
    reshape_size[index] = len(data.coords[dim])

    # Must include factor of 2 in exponential to get correct linewidth ->
    window_array = _np.exp(-1.0 * data.coords[dim] * 2.0 * linewidth).reshape(
        reshape_size
    )
    window_array = _np.ones_like(data.values) * window_array
    data.values *= window_array

    proc_attr_name = "window"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data


def integrate(all_data, dim="t2", integrate_center=0, integrate_width=100):
    """Integrate data down given dimension

    Args:
        all_data (dnpdata,dict): Data container
        proc_parameters (dict, procParam): Processing Parameters

    +------------------+-------+---------+------------------------------+
    | parameter        | type  | default | description                  |
    +------------------+-------+---------+------------------------------+
    | dim              | str   | 't2'    | dimension to integrate       |
    +------------------+-------+---------+------------------------------+
    | integrate_center | float | 0       | center of integration window |
    +------------------+-------+---------+------------------------------+
    | integrate_width  | float | 100     | width of integration window  |
    +------------------+-------+---------+------------------------------+

    Returns:
        all_data (dnpdata,dict): Processed data

    Example::

        dnplab.dnpNMR.integrate(all_data)

    """

    data, isDict = return_data(all_data)

    proc_parameters = {
        "dim": dim,
        "integrate_center": integrate_center,
        "integrate_width": integrate_width,
    }

    integrateMin = integrate_center - _np.abs(integrate_width) / 2.0
    integrateMax = integrate_center + _np.abs(integrate_width) / 2.0

    data = data[dim, (integrateMin, integrateMax)]

    index = data.index(dim)
    data.values = _np.trapz(data.values, x=data.coords[dim], axis=index)

    data.coords.pop(dim)

    proc_attr_name = "integrate"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data


def align(all_data, dim="t2"):
    """Alignment of NMR spectra down given dim dimension

    Example::

        data = dnp.dnpNMR.align(data, {})
    """

    data, isDict = return_data(all_data)

    if len(_np.shape(data.values)) != 2:
        raise ValueError("Only 2-dimensional data is currently supported")

    proc_parameters = {"dim": dim}
    originalAxesOrder = data.dims
    data.reorder([dim])
    dimIter = data.dims[-1]

    refData = data[dimIter, 0].values.reshape(-1)

    for ix in range(len(data.coords[dimIter])):
        tempData = data[dimIter, ix].values.reshape(-1)

        corrData = _np.correlate(_np.abs(tempData), _np.abs(refData), mode="same")
        shiftIx = _np.argmax(corrData) - (
            len(corrData) / 2
        )  # subtract half length so spectrum is shifted relative to center, not edge
        shiftData = _np.roll(tempData, -1 * int(_np.round(shiftIx, 0)))
        data.values[:, ix] = shiftData
    data.reorder(originalAxesOrder)

    proc_attr_name = "align"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data


def autophase(workspace, method="arctan"):
    """Automatically phase data

    Args:
        workspace (dnpdata_collection, dnpdata): Data object to autophase
        parameters (dict):

    Returns:
        dnpdata_collection, dnpdata: Autophased data

    Example::

        ws = dnp.dnpNMR.autophase(ws, {})

    """

    proc_parameters = {"method": "arctan"}

    data, is_workspace = return_data(workspace)
    phase = _np.arctan(_np.sum(_np.imag(data.values)) / _np.sum(_np.real(data.values)))

    data.values *= _np.exp(-1j * phase)
    if _np.sum(_np.real(data.values)) < 0:
        data.values *= -1.0

    proc_attr_name = "autophase"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if is_workspace:
        workspace[workspace.processing_buffer] = data
        return workspace
    else:
        return data
