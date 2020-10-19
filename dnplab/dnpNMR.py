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


def exp_fit_func_1(x_axis, c1, c2, c3):
    return c1 + c2 * _np.exp(-1.0 * c3 * x_axis)


def exp_fit_func_2(x_axis, c1, c2, c3, c4, c5):
    return c1 + c2 * _np.exp(-1.0 * c3 * x_axis) + c4 * _np.exp(-1.0 * c5 * x_axis)


def baseline_fit(temp_coords, temp_data, type, order):

    if type == "poly":
        base_line = _np.polyval(_np.polyfit(temp_coords, temp_data, order), temp_coords)
    elif type == "expon":
        temp_data = temp_data.real
        if order == 1:
            x0 = [temp_data[-1], temp_data[0], 1]
            out, cov = curve_fit(
                exp_fit_func_1, temp_coords, temp_data, x0, method="lm"
            )
            base_line = exp_fit_func_1(temp_coords, out[0], out[1], out[2])
        elif order == 2:
            x0 = [temp_data[-1], temp_data[0], 1, temp_data[0], 1]
            out, cov = curve_fit(
                exp_fit_func_2, temp_coords, temp_data, x0, method="lm"
            )
            base_line = exp_fit_func_2(
                temp_coords, out[0], out[1], out[2], out[3], out[4]
            )
        else:
            raise ValueError(
                "Use order=1 for mono-exponential, order=2 for bi-exponential"
            )

    else:
        raise TypeError

    return base_line


def baseline(all_data, dim="t1", type="poly", order=1, reference_slice=None):
    """Baseline correction of NMR spectra down given dim dimension

    Args:
        all_data (object) : dnpdata object
        dim (str) : dimension to correct along
        type (str) : type of baseline fit, either "poly" for polynomial or "expon" for exponential
        order (int) : polynomial order, or 1=mono and 2=bi for exponential


    returns:
        all_data (dnpdata, dict): Baseline corrected data in container
    """

    data, isDict = return_data(all_data)

    if reference_slice is not None:
        if len(_np.shape(data.values)) == 1:
            reference_slice = None
            print("ignoring reference_slice, this is 1D data")
        else:
            reference_slice -= 1

    if len(_np.shape(data.values)) == 2:
        temp_coords = data[dim, :].coords["t2"]
        if reference_slice is not None:
            bline = baseline_fit(
                temp_coords, data[dim, :].values[:, reference_slice], type, order
            )
            for ix in range(len(data.coords[dim]) - 1):
                data.values[:, ix] -= bline
        elif reference_slice is None:
            for ix in range(len(data.coords[dim]) - 1):
                bline = baseline_fit(temp_coords, data.values[:, ix], type, order)
                data.values[:, ix] -= bline
        else:
            raise TypeError("invalid reference_slice")

    elif len(_np.shape(data.values)) == 1:
        dim = "t2"
        bline = baseline_fit(data[dim, :].coords[dim], data[dim, :].values, type, order)
        data.values -= bline

    else:
        raise ValueError("1D or 2D only")

    proc_parameters = {
        "dim": dim,
        "type": type,
        "order": order,
        "reference_slice": reference_slice,
    }
    proc_attr_name = "baseline"
    data.add_proc_attrs(proc_attr_name, proc_parameters)
    data.attrs["baseline"] = bline

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


def autophase(
    all_data,
    method="arctan",
    reference_slice=None,
    force_positive=False,
    order="zero",
    pivot=0,
    rise=0,
):
    """Automatically phase data

    Args:
        all_data (dnpdata_collection, dnpdata): Data object to autophase
        method (str): "arctan" finds the arctan of the ratio of the sum of the imaginary to the sum of the real, or "search" finds the maximum of the sum of the real to the sum of the imaginary after phase correction by an array of test angles
        reference_slice (int, or None): slice of 2D data that is used to draw the baseline, none means each slice is corrected individually
        force_positive (boolean): If true, any phase correction orients data above the baseline
        order (str) : "zero" or "first" order phase corrections
        pivot (int) : number between 0 and the length of the data at which the first order phase correction is equal to the zeroth order phase correction
        rise (float) : change in phase angle from data[0] to data[end]

    Returns:
        all_data (dnpdata, dict): Autophased data in container

    """

    data, isDict = return_data(all_data)
    shape_data = _np.shape(data.values)

    if reference_slice is not None:
        if len(shape_data) == 1:
            reference_slice = None
            temp_data = data.values
            print("ignoring reference_slice, this is 1D data")
        else:
            reference_slice -= 1
            temp_data = data.values[:, reference_slice]
    else:
        temp_data = data.values

    if method == "arctan":
        phase_0 = _np.arctan(
            _np.sum(_np.imag(temp_data.reshape(-1, 1)))
            / _np.sum(_np.real(temp_data.reshape(-1, 1)))
        )
    elif method == "search":
        phases_0 = _np.linspace(-_np.pi / 2, _np.pi / 2, 180).reshape(-1)
        rotated_data = (temp_data.reshape(-1, 1)) * _np.exp(-1j * phases_0)
        real_imag_ratio = (_np.real(rotated_data) ** 2).sum(axis=0) / (
            (_np.imag(rotated_data) ** 2).sum(axis=0)
        )
        phase_0 = phases_0[_np.argmax(real_imag_ratio)]
    else:
        raise TypeError("Invalid method")

    if order == "zero":
        data.values *= _np.exp(-1j * phase_0)
    elif order == "first":
        pivot_ratio = pivot / len(data.values)
        phase_1 = _np.linspace(
            phase_0 - rise * pivot_ratio,
            phase_0 + rise * (1 - pivot_ratio),
            len(data.values),
        )
        data.attrs["phase_1"] = phase_1
        if len(shape_data) == 2:
            for ix in range(shape_data[1] - 1):
                data.values[:, ix] *= _np.exp(-1j * phase_1)
        else:
            data.values *= _np.exp(-1j * phase_1)
    else:
        raise TypeError("Invalid order")

    if force_positive:
        if len(shape_data) == 2:
            for ix in range(shape_data[1] - 1):
                if _np.sum(_np.real(data.values[:, ix])) < 0:
                    data.values[:, ix] *= -1.0
                else:
                    pass
        elif len(_np.shape(data.values)) == 1:
            if _np.sum(_np.real(data.values)) < 0:
                data.values *= -1.0
        else:
            raise ValueError("only 1D or 2D data are currently supported")

    proc_parameters = {
        "method": method,
        "reference_slice": reference_slice,
        "force_positive": force_positive,
        "order": order,
        "pivot": pivot,
        "rise": rise,
    }
    proc_attr_name = "autophase"
    data.add_proc_attrs(proc_attr_name, proc_parameters)
    data.attrs["phase_0"] = phase_0

    if isDict:
        all_data[all_data.processing_buffer] = data
        return all_data
    else:
        return data
