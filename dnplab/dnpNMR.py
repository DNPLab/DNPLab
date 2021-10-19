""" Functions specific to processing NMR data
"""
from warnings import warn

import numpy as _np

from . import return_data, dnpdata, dnpdata_collection
from . import dnpTools, dnpMath

import re
import copy


def ndalign(all_data, dim="f2", reference=None, center=None, width=None, average=None):
    """Alignment of NMR spectra using FFT Cross Correlation

    Args:
        all_data (object) : dnpdata object
        dim (str) : dimension to align along
        reference (numpy) : second dimension to align along
        center (float) : range center
        width (float) : range width

    returns:
        dnpdata: Aligned data in container
    """

    data, isDict = return_data(all_data)

    proc_parameters = {"dim": dim}

    original_order = data.dims  # Original order of dims

    data.reorder([dim])  # Move dim to first dimension

    all_values = data.values  # Extract Data Values for alignment

    if center != None and width != None:
        start = center - 0.5 * width
        stop = center + 0.5 * width
    elif center == None and width == None:
        start = data.coords[dim][-1]
        stop = data.coords[dim][0]
    else:
        raise ValueError("selected rangfe is not accpetable")

    values = data[dim, (start, stop)].values

    all_original_shape = all_values.shape
    original_shape = values.shape  # Preserve original shape

    all_align_dim_length = all_original_shape[0]
    align_dim_length = original_shape[0]  # length of dimension to align down

    all_values = all_values.reshape(all_align_dim_length, -1)
    values = values.reshape(align_dim_length, -1)  # Reshape to 2d

    new_shape = _np.shape(values)

    dim2 = new_shape[1]

    abs_values = _np.abs(values)

    if reference is None:
        reference = _np.abs(values[:, -1])
    elif isinstance(reference, dnpdata):
        reference = _np.abs(reference.values)
        if average != None:
            reference = _np.convolve(reference, _np.ones(average), "same") / average

    ref_max_ix = _np.argmax(reference)

    all_aligned_values = _np.zeros_like(all_values)

    for ix in range(dim2):
        if average != None:
            abs_values[:, ix] = (
                _np.convolve(abs_values[:, ix], _np.ones(average), "same") / average
            )
        cor = _np.correlate(
            abs_values[:, ix], reference, mode="same"
        )  # calculate cross-correlation
        max_ix = _np.argmax(cor)  # Maximum of cross correlation
        delta_max_ix = max_ix - ref_max_ix  # Calculate how many points to shift
        all_aligned_values[:, ix] = _np.roll(
            all_values[:, ix], -1 * delta_max_ix
        )  # shift values

    all_aligned_values = all_aligned_values.reshape(
        all_original_shape
    )  # reshape to original values shape

    data.values = all_aligned_values  # Add aligned values back to data object

    data.reorder(original_order)  # Back to original order

    proc_attr_name = "ndalign"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def align(all_data, dim="f2", dim2=None, center=None, width=None):
    """Alignment of NMR spectra down given dimension or dimensions

    Args:
        all_data (object) : dnpdata object
        dim (str) : dimension to align along
        dim2 (str) : second dimension to align along
        center (float) : range center
        width (float) : range width

    returns:
        dnpdata: Aligned data in container
    """

    data, isDict = return_data(all_data)

    if len(_np.shape(data.values)) > 3:
        raise ValueError("Greater than 3-dimensional data is currently not supported")

    proc_parameters = {"dim": dim, "dim2": dim2}
    originalAxesOrder = data.dims

    if (dim2 is None) or (len(data.shape) == 2):
        if len(data.shape) > 2:
            raise ValueError(
                "data has more than 2 dimensions, 2nd dimension is ambiguous"
            )
        dim2 = None
        data.reorder([dim])
        dimIter = data.dims[-1]

    else:
        data.reorder([dim, dim2])
        dimIter = data.dims[-1]
    if center != None and width != None:
        start = center - 0.5 * width
        stop = center + 0.5 * width
    elif center == None and width == None:
        start = None
        stop = None
    else:
        raise ValueError("selected range is not accpetale")
    if dim2 == None:
        if start != None and stop != None:
            refData = data[dimIter, 0, dim, (start, stop)].values.reshape(-1)
        elif start == None and stop == None:
            refData = data[dimIter, 0].values.reshape(-1)
        else:
            raise ValueError("selected range is not accpetale")

        for ix in range(len(data.coords[dimIter])):
            tempData = data[dimIter, ix].values.reshape(-1)
            if start != None and stop != None:
                rangeData = data[dimIter, ix, dim, (start, stop)].values.reshape(-1)
            elif start == None and stop == None:
                rangeData = tempData
            else:
                raise ValueError("selected range is not accpetale")

            corrData = _np.correlate(_np.abs(rangeData), _np.abs(refData), mode="same")
            shiftIx = _np.argmax(corrData) - (
                len(corrData) / 2
            )  # subtract half length so spectrum is shifted relative to center, not edge
            shiftData = _np.roll(tempData, -1 * int(_np.round(shiftIx, 0)))
            data.values[:, ix] = shiftData
    else:

        for ix1 in range(len(data.coords[-1])):
            if start != None and stop != None:
                refData = data[dim, (start, stop)].values[:, 0, 0]
            elif start == None and stop == None:
                refData = data.values[:, 0, 0]
            else:
                raise ValueError("selected range is not accpetale")

            for ix2 in range(len(data.coords[dim2])):
                tempData = data.values[:, ix2, ix1]
                if start != None and stop != None:
                    rangeData = data[dim, (start, stop)].values[:, ix2, ix1]
                elif start == None and stop == None:
                    rangeData = tempData
                else:
                    raise ValueError("selected range is not accpetale")

                corrData = _np.correlate(
                    _np.abs(rangeData), _np.abs(refData), mode="same"
                )
                shiftIx = _np.argmax(corrData) - (
                    len(corrData) / 2
                )  # subtract half length so spectrum is shifted relative to center, not edge
                shiftData = _np.roll(tempData, -1 * int(_np.round(shiftIx, 0)))
                data.values[:, ix2, ix1] = shiftData

    data.reorder(originalAxesOrder)

    proc_attr_name = "align"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def autophase(
    all_data,
    dim="f2",
    method="search",
    reference_range=None,
    pts_lim=None,
    order="zero",
    pivot=0,
    delta=0,
    phase=None,
    reference_slice=None,
    force_positive=False,
):
    """Automatically phase correct data, or apply manual phase correction

    .. math::

        \mathrm{data}         &= \exp(-1j * \mathrm{phase}) &

        \mathrm{phase(arctan)} &= \mathrm{arctan}(\mathrm{sum}(\mathrm{data.imag}) / \mathrm{sum}(\mathrm{data.real})) &

        \mathrm{phase(search)} &= \mathrm{argmax}(\mathrm{sum}(phased\_real^{2}) / \mathrm{sum}(phased\_imag^{2})) &

        phased\_real          &= \mathrm{data.real} * \exp(-1j * \mathrm{phase}) &

        phased\_imag          &= \mathrm{data.imag} * \exp(-1j * \mathrm{phase}) &

    Args:
        all_data (dnpdata_collection, dnpdata): Data object to autophase

    +-----------------+---------------+---------------+---------------------------------------------------+
    | parameter       | type          | default       | description                                       |
    +=================+===============+===============+===================================================+
    | method          | str           | 'search'      | method of searching for the best phase            |
    +-----------------+---------------+---------------+---------------------------------------------------+
    | reference_range | list or tuple | None          | data window to use for phase calculation          |
    +-----------------+---------------+---------------+---------------------------------------------------+
    | pts_lim         | int or None   | None          | specify the max points used in phase search       |
    +-----------------+---------------+---------------+---------------------------------------------------+
    | order           | str           | 'zero'        | order of phase correction                         |
    +-----------------+---------------+---------------+---------------------------------------------------+
    | pivot           | int           | 0             | pivot point for first order correction            |
    +-----------------+---------------+---------------+---------------------------------------------------+
    | delta           | float or int  | 0             | total change in phase magnitude for first order   |
    +-----------------+---------------+---------------+---------------------------------------------------+
    | phase           | float or int  | 0             | manual phase correction in radians                |
    +-----------------+---------------+---------------+---------------------------------------------------+
    | reference_slice | int, or None  | None          | slice of 2D data used to define the phase         |
    +-----------------+---------------+---------------+---------------------------------------------------+
    | force_positive  | boolean       | False         | force the entire spectrum to positive magnitude   |
    +-----------------+---------------+---------------+---------------------------------------------------+

    Returns:
        dnpdata: Autophased data, including attrs "phase0" for order="zero", and "phase1" if order="first"

    """
    if reference_slice == 0:
        raise ValueError(
            "please use indices from 1 to number of slices, i.e. use 1 instead of 0"
        )

    data, isDict = return_data(all_data)
    shape_data = _np.shape(data.values)
    index = data.dims.index(dim)

    if phase is not None:
        method = "manual"

    if method == "manual":
        if order == "zero" and isinstance(phase, (int, float)):
            data.attrs["phase0"] = phase
        elif order == "zero" and not isinstance(phase, (int, float)):
            raise ValueError(
                "for a zero order phase correction you must supply a single phase"
            )
        elif order == "first" and isinstance(phase, (int, float)):
            data.attrs["phase0"] = phase
            order = "zero"
            warn(
                "method=manual and order=first but only a single phase was given, switching to order=zero"
            )
        elif (
            order == "first"
            and isinstance(phase, (list, _np.ndarray))
            and len(phase) == shape_data[index]
        ):
            data.attrs["phase1"] = _np.array(phase)
        else:
            raise ValueError(
                "Invalid combination of phase order and phase value(s). Supply float for zero order, array or list for first order"
            )
    else:

        if isinstance(reference_range, (list, tuple)) and len(reference_range) == 2:
            check_data = data[dim, (reference_range[0], reference_range[1])]
        else:
            check_data = data.copy()
            if reference_range is not None:
                warn("reference_range must be None or list/tuple length=2")

        if reference_slice is not None:
            if len(shape_data) == 1:
                reference_slice = None
                temp_data = check_data.values
                warn("ignoring reference_slice, this is 1D data")
            else:
                temp_data = check_data.values[:, reference_slice - 1]
        else:
            temp_data = check_data.values

        if method == "arctan":
            data.attrs["phase0"] = _np.arctan(
                _np.sum(_np.imag(temp_data.reshape(-1, 1)))
                / _np.sum(_np.real(temp_data.reshape(-1, 1)))
            )
        elif method == "search":
            if pts_lim is not None:
                if len(check_data.coords[dim]) > pts_lim:
                    phasing_x = _np.linspace(
                        min(check_data.coords[dim]),
                        max(check_data.coords[dim]),
                        int(pts_lim),
                    ).reshape(-1)
                    if len(check_data.dims) == 1:
                        temp_data = _np.interp(
                            phasing_x, check_data.coords[dim], check_data.values
                        ).reshape(-1)
                    else:
                        ind_dim = list(set(data.dims) - set([dim]))[0]
                        ind_shape = data.shape[data.index(ind_dim)]
                        temp_data = _np.array(
                            [
                                _np.interp(
                                    phasing_x,
                                    check_data.coords[dim],
                                    check_data[dim, :].values[:, indx],
                                ).reshape(-1)
                                for indx in range(ind_shape)
                            ]
                        ).reshape(pts_lim, ind_shape)
            phases_0 = _np.linspace(-_np.pi / 2, _np.pi / 2, 180).reshape(-1)
            rotated_data = (temp_data.reshape(-1, 1)) * _np.exp(-1j * phases_0)
            real_imag_ratio = (_np.real(rotated_data) ** 2).sum(axis=0) / (
                (_np.imag(rotated_data) ** 2).sum(axis=0)
            )
            data.attrs["phase0"] = phases_0[_np.argmax(real_imag_ratio)]
        else:
            raise TypeError("Invalid autophase method")

    if order == "zero":
        data.values *= _np.exp(-1j * data.attrs["phase0"])
    elif order == "first":
        if method == "manual":
            data.attrs["phase1"] = phase
        else:
            pivot_ratio = pivot / len(data.values)
            data.attrs["phase1"] = _np.linspace(
                data.attrs["phase0"] - delta * pivot_ratio,
                data.attrs["phase0"] + delta * (1 - pivot_ratio),
                len(data.values),
            )
        data.values.T.dot(_np.exp(-1j * data.attrs["phase1"]))

    else:
        raise TypeError("Invalid order or order & phase pair")

    if force_positive:
        data.values = _np.absolute(data.values)

    proc_parameters = {
        "method": method,
        "reference_range": reference_range,
        "reference_slice": reference_slice,
        "force_positive": force_positive,
        "order": order,
        "pivot": pivot,
        "delta": delta,
    }
    proc_attr_name = "autophase"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def calculate_enhancement(
    all_data,
    off_spectrum=1,
    on_spectra="all",
    integrate_center=0,
    integrate_width="full",
    method="integrate",
    dim="f2",
):
    """Calculate enhancement from DNP data

    Args:
        all_data (dnpdata, dict): data container

    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | parameter        | type                       | default     | description                                                          |
    +==================+============================+=============+======================================================================+
    | off_spectrum     | int or dnpdata             | 1           | slice of 2D data to be used as p = 0 spectrum, or dnpdata            |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | on_spectra       | str or dnpdata             | "all"       | "all"  unless dnpdata given                                          |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | integrate_center | str, int, list, or ndarray | 0           | "max", center of integration window, or index used to find amplitude |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | integrate_width  | str, int, list, or ndarray | "full"      | "full" or width of integration window                                |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | method           | str                        | "integrate" | either "integrate" or "ampltiude"                                    |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+
    | dim              | str                        | "f2"        | dimension to integrate down or search down for max                   |
    +------------------+----------------------------+-------------+----------------------------------------------------------------------+

    Returns:
        dnpdata: data object with "enhancements" key added

    """
    if off_spectrum == 0:
        raise ValueError(
            "please use indices from 1 to number of slices, i.e. use 1 instead of 0"
        )

    _, isDict = return_data(all_data)

    if isDict and "integrals" in all_data.keys() and isinstance(off_spectrum, int):
        enh = _np.array(
            all_data["integrals"].values.real
            / all_data["integrals"].values.real[off_spectrum - 1]
        )
        enhancement_data = dnpdata(
            enh,
            [all_data["integrals"].coords[x] for x in all_data["integrals"].dims],
            all_data["integrals"].dims,
        )

    elif isinstance(off_spectrum, dnpdata) and isinstance(on_spectra, dnpdata):

        data_off, _ = return_data(off_spectrum)
        if len(data_off.shape) != 1:
            raise TypeError("off_spectrum should be 1D")
        data_on, _ = return_data(on_spectra)
        index_on = data_on.dims.index(dim)

        if integrate_width == "full":
            int_width_off = max(data_off.coords[dim]) - min(data_off.coords[dim])
            int_width_on = max(data_on.coords[dim]) - min(data_on.coords[dim])
        elif (
            isinstance(integrate_width, (list, _np.ndarray))
            and len(integrate_width) == 2
        ):
            int_width_off = integrate_width[0]
            int_width_on = integrate_width[1]
        elif isinstance(integrate_width, int):
            int_width_off = integrate_width
            int_width_on = integrate_width
        else:
            raise ValueError(
                "integrate_width must be an integer, a list/array with len = 2, or 'full'"
            )

        if integrate_center == "max":
            on_maxs = _np.argmax(data_on.values.real, axis=index_on)
            int_center_off = data_off.coords[dim][_np.argmax(data_off.values.real)]
            if on_maxs.size == 1:
                int_center_on = data_on.coords[dim][on_maxs]
            else:
                int_center_on = [data_on.coords[dim][x] for x in on_maxs]
        elif (
            isinstance(integrate_center, (list, _np.ndarray))
            and len(integrate_center) == 2
        ):
            int_center_off = integrate_center[0]
            int_center_on = integrate_center[1]
        elif isinstance(integrate_center, int):
            int_center_off = integrate_center
            int_center_on = integrate_center
        else:
            raise ValueError(
                "integrate_center must be an integer, a list/array with len = 2, or 'max'"
            )

        if method == "integrate":
            off_data = dnpTools.integrate(
                data_off,
                dim=dim,
                integrate_center=int_center_off,
                integrate_width=int_width_off,
            )

            on_data = dnpTools.integrate(
                data_on,
                dim=dim,
                integrate_center=int_center_on,
                integrate_width=int_width_on,
            )

            enh = _np.array(on_data.values.real / off_data.values.real)
            enhancement_data = dnpdata(
                enh, [on_data.coords[x] for x in on_data.dims], on_data.dims
            )

        elif method == "amplitude":
            on_maxs = _np.argmax(abs(data_on.values.real), axis=index_on)
            if integrate_center == "max":
                off_data = data_off.values.real[_np.argmax(abs(data_off.values.real))]
                if on_maxs.size == 1:
                    on_data = data_on.values.real[on_maxs]
                else:
                    on_data = [
                        data_on.values.real[x, indx] for indx, x in enumerate(on_maxs)
                    ]
            else:
                off_data = data_off.values.real[int_center_off]
                if on_maxs.size == 1:
                    on_data = data_on.values.real[int_center_on]
                else:
                    on_data = [
                        data_on.values.real[int_center_on, indx]
                        for indx, _ in enumerate(on_maxs)
                    ]

            if (isinstance(on_data, list) and len(on_data) == 1) or (
                isinstance(on_data, float) and on_data.size == 1
            ):
                enh = _np.array([on_data / off_data])
            else:
                enh = _np.array(on_data / off_data)

            remaining_dims = [x for x in data_on.dims if x != dim]
            if len(remaining_dims) == 0:
                remaining_dims = ["index"]
                remaining_coords = [_np.array([0])]
            else:
                remaining_coords = [data_on.coords[x] for x in data_on.dims if x != dim]

            enhancement_data = dnpdata(enh, remaining_coords, remaining_dims)

    else:
        raise TypeError(
            "Either use the integrate function first and define the index of the off spectrum, or pass dnpata objects for off_spectrum and on_spectra"
        )

    if isDict:
        all_data["enhancements"] = enhancement_data
    else:
        return enhancement_data


def fourier_transform(
    all_data,
    dim="t2",
    zero_fill_factor=2,
    shift=True,
    convert_to_ppm=True,
    output="complex",
):
    """Perform Fourier Transform down dim dimension given in proc_parameters

    .. Note::
        Assumes dt = t[1] - t[0]

    .. math::

        \mathrm{power spectrum}       &=  data.real^{2} + data.imag^{2} &

        \mathrm{magnitude spectrum}   &=  sqrt(data.real^{2} + data.imag^{2}) &

    Args:
        all_data (dnpdata, dict): Data container

    +------------------+------+-----------+--------------------------------------------------+
    | parameter        | type | default   | description                                      |
    +==================+======+===========+==================================================+
    | dim              | str  | 't2'      | dimension to Fourier transform                   |
    +------------------+------+-----------+--------------------------------------------------+
    | zero_fill_factor | int  | 2         | factor to increase dim with zeros                |
    +------------------+------+-----------+--------------------------------------------------+
    | shift            | bool | True      | Perform fftshift to set zero frequency to center |
    +------------------+------+-----------+--------------------------------------------------+
    | convert_to_ppm   | bool | True      | Convert dim from Hz to ppm                       |
    +------------------+------+-----------+--------------------------------------------------+
    | output           | str  | 'complex' | output complex, magnitude, or power spectrum     |
    +------------------+------+-----------+--------------------------------------------------+

    Returns:
        dnpdata: data object after FT
    """
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
        if "nmr_frequency" not in data.attrs.keys():
            warn(
                "NMR frequency not found in the attrs dictionary, coversion to ppm requires the NMR frequency. See docs."
            )
        else:
            nmr_frequency = data.attrs["nmr_frequency"]
            f /= -1 * nmr_frequency / 1.0e6

    data.values = _np.fft.fft(data.values, n=n_pts, axis=index)

    if shift:
        data.values = _np.fft.fftshift(data.values, axes=index)

    data.coords[dim] = f

    if output == "mag":
        data.values = _np.sqrt(data.values.real ** 2 + data.values.imag ** 2)
    elif output == "pow":
        data.values = data.values.real ** 2 + data.values.imag ** 2
    elif output == "complex":
        pass
    else:
        raise ValueError("options for output are 'complex' (default), 'mag', or 'pow'")

    if re.fullmatch("t[0-9]*", dim) is not None:
        new_dim = dim.replace("t", "f")
        data.rename(dim, new_dim)

    proc_attr_name = "fourier_transform"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def inverse_fourier_transform(
    all_data,
    dim="f2",
    zero_fill_factor=1,
    shift=True,
    convert_from_ppm=True,
    output="complex",
):
    """Perform Fourier Transform down dim dimension given in proc_parameters

    .. Note::
        Assumes dt = t[1] - t[0]

    .. math::

        \mathrm{power spectrum}       &=  data.real^{2} + data.imag^{2} &

        \mathrm{magnitude spectrum}   &=  sqrt(data.real^{2} + data.imag^{2}) &

    Args:
        all_data (dnpdata, dict): Data container

    +------------------+------+-----------+--------------------------------------------------+
    | parameter        | type | default   | description                                      |
    +==================+======+===========+==================================================+
    | dim              | str  | 't2'      | dimension to Fourier transform                   |
    +------------------+------+-----------+--------------------------------------------------+
    | zero_fill_factor | int  | 2         | factor to increase dim with zeros                |
    +------------------+------+-----------+--------------------------------------------------+
    | shift            | bool | True      | Perform fftshift to set zero frequency to center |
    +------------------+------+-----------+--------------------------------------------------+
    | convert_from_ppm | bool | True      | Convert dim from Hz to ppm                       |
    +------------------+------+-----------+--------------------------------------------------+
    | output           | str  | 'complex' | output complex, magnitude, or power spectrum     |
    +------------------+------+-----------+--------------------------------------------------+

    Returns:
        dnpdata: data object after inverse FT
    """
    data, isDict = return_data(all_data)

    # handle zero_fill_factor
    zero_fill_factor = int(zero_fill_factor)
    if zero_fill_factor <= 0:
        zero_fill_factor = 1

    proc_parameters = {
        "dim": dim,
        "zero_fill_factor": zero_fill_factor,
        "shift": shift,
        "convert_from_ppm": convert_from_ppm,
    }

    index = data.dims.index(dim)

    df = data.coords[dim][1] - data.coords[dim][0]
    if convert_from_ppm:
        if "nmr_frequency" not in data.attrs.keys():
            warn(
                "NMR frequency not found in the attrs dictionary, coversion from ppm requires the NMR frequency. See docs."
            )
        else:
            nmr_frequency = data.attrs["nmr_frequency"]
            df /= -1 / (nmr_frequency / 1.0e6)

    n_pts = zero_fill_factor * len(data.coords[dim])
    t = (1.0 / (n_pts * df)) * _np.r_[0:n_pts]

    if shift:
        data.values = _np.fft.fftshift(data.values, axes=index)

    data.values = _np.fft.ifft(data.values, n=n_pts, axis=index)
    data.coords[dim] = t

    if output == "mag":
        data.values = _np.sqrt(data.values.real ** 2 + data.values.imag ** 2)
    elif output == "pow":
        data.values = data.values.real ** 2 + data.values.imag ** 2
    elif output == "complex":
        pass
    else:
        raise ValueError("options for output are 'complex' (default), 'mag', or 'pow'")

    if re.fullmatch("f[0-9]*", dim) is not None:
        new_dim = dim.replace("f", "t")
        data.rename(dim, new_dim)

    proc_attr_name = "inverse_fourier_transform"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def left_shift(all_data, dim="t2", shift_points=0):
    """Remove points from the left of data

    Args:
        all_data (dnpdata, dict): Data container for data

    +---------------+------+---------+--------------------------------------------------+
    | parameter     | type | default | description                                      |
    +===============+======+=========+==================================================+
    | dim           | str  | "t2"    | dimension to shift                               |
    +---------------+------+---------+--------------------------------------------------+
    | shift_points  | int  | 0       | Number of points to remove from left of data     |
    +---------------+------+---------+--------------------------------------------------+

    Returns:
        dnpdata: data object with left-shifted data
    """

    data, isDict = return_data(all_data)

    shifted_data = data[dim, shift_points:]

    proc_attr_name = "left_shift"
    proc_parameters = {
        "dim": dim,
        "points": shift_points,
    }
    shifted_data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = shifted_data
    else:
        return shifted_data


def remove_offset(all_data, dim="t2", offset_points=10):
    """Remove DC offset from FID by averaging the last few data points and subtracting the average

    Args:
        all_data (dnpdata, dict): Data container for data

    +---------------+------+---------+----------------------------------------------------------+
    | parameter     | type | default | description                                              |
    +===============+======+=========+==========================================================+
    | dim           | str  | 't2'    | Dimension to calculate DC offset                         |
    +---------------+------+---------+----------------------------------------------------------+
    | offset_points | int  | 10      | Number of points at end of data to average for DC offset |
    +---------------+------+---------+----------------------------------------------------------+

    Returns:
        dnpdata: data object with offset removed
    """

    # Determine if data is dictionary or dnpdata object
    data, isDict = return_data(all_data)

    proc_parameters = {
        "dim": dim,
        "offset_points": offset_points,
    }

    dim = proc_parameters["dim"]
    offset_points = int(proc_parameters["offset_points"])

    offsetData = data[dim, -1 * offset_points :].values
    offsetData = offsetData.reshape(-1)
    offset = _np.mean(offsetData)

    data -= offset

    proc_attr_name = "remove_offset"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def window(
    all_data,
    dim="t2",
    type="exponential",
    linewidth=1,
    gaussian_max=0,
    inverse=False,
):
    """Apply Apodization to data down given dimension

    .. math::

        \mathrm{exponential}    &=  \exp(-2t * \mathrm{linewidth}) &

        \mathrm{gaussian}       &=  \exp((\mathrm{linewidth[0]} * t) - (\mathrm{linewidth[1]} * t^{2})) &

        \mathrm{hamming}        &=  0.53836 + 0.46164\cos(\pi * n/(N-1)) &

        \mathrm{han}            &=  0.5 + 0.5\cos(\pi * n/(N-1)) &

        \mathrm{sin2}           &=  \cos((-0.5\pi * n/(N - 1)) + \pi)^{2} &

        \mathrm{lorentz\_gauss} &=  \exp(L -  G^{2}) &

               L(t)    &=  \pi * \mathrm{linewidth[0]} * t &

               G(t)    &=  0.6\pi * \mathrm{linewidth[1]} * (\mathrm{gaussian\_max} * (N - 1) - t) &

        \mathrm{traf}           &=  (f1 * (f1 + f2)) / (f1^{2} + f2^{2}) &

               f1(t)   &=  \exp(-t * \pi * \mathrm{linewidth[0]}) &

               f2(t)   &=  \exp((t - T) * \pi * \mathrm{linewidth[1]}) &


    Args:
        all_data (dnpdata, dict): data container

    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | parameter       | type                    | default       | description                                       |
    +=================+=========================+===============+===================================================+
    | dim             | str                     | 't2'          | Dimension to apply exponential apodization        |
    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | type            | str                     | 'exponential' | type of apodization                               |
    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | linewidth       | float, list, or ndarray | 1             | linewidths  in Hz                                 |
    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | gaussian_max    | float                   | 0             | Location of gaussian component maximum            |
    +-----------------+-------------------------+---------------+---------------------------------------------------+
    | inverse         | boolean                 | False         | invert the window function                        |
    +-----------------+-------------------------+---------------+---------------------------------------------------+

    Returns:
        dnpdata: data object with window function applied, including attr "window"
    """
    data, isDict = return_data(all_data)
    dim_size = data.coords[dim].size
    shape_data = _np.shape(data.values)
    index = data.index(dim)

    if (isinstance(linewidth, _np.ndarray) or isinstance(linewidth, list)) and len(
        linewidth
    ) == 2:
        exp_lw = linewidth[0]
        gauss_lw = linewidth[1]
    elif isinstance(linewidth, (int, float)):
        exp_lw = linewidth
        gauss_lw = linewidth
    else:
        raise ValueError("linewidth must be int/float, or list/ndarray with len==2")

    if type == "exponential":
        apwin = dnpMath.exponential_window(all_data, dim, linewidth)
    elif type == "gaussian":
        apwin = dnpMath.gaussian_window(all_data, dim, [exp_lw, gauss_lw])
    elif type == "hamming":
        apwin = dnpMath.hamming_window(dim_size)
    elif type == "hann":
        apwin = dnpMath.hann_window(dim_size)
    elif type == "lorentz_gauss":
        apwin = dnpMath.lorentz_gauss_window(
            all_data, dim, exp_lw, gauss_lw, gaussian_max=gaussian_max
        )
    elif type == "sin2":
        apwin = dnpMath.sin2_window(dim_size)
    elif type == "traf":
        apwin = dnpMath.traf_window(all_data, dim, exp_lw, gauss_lw)
    else:
        raise ValueError("Invalid window type")

    apwin.reshape(dim_size)

    if inverse:
        apwin = 1 / apwin

    new_shape = [1 if ix != index else shape_data[index] for ix in range(data.ndim)]
    apwin = _np.reshape(apwin, new_shape)

    data.values *= apwin

    proc_parameters = {
        "type": type,
        "linewidth": linewidth,
        "dim": dim,
        "gaussian_max": gaussian_max,
        "inverse": inverse,
    }
    proc_attr_name = "window"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    if isDict:
        all_data[all_data.processing_buffer] = data
    else:
        return data


def phasecycle(all_data, dim, receiver_phase):
    """Phase cycle

    Args:
        all_data (dnpdata_collection, dnpdata): data to process
        dim (str): dimension to perform phase cycle
        receiver_phase (numpy.array, list): Receiver Phase 0 (x), 1 (y), 2 (-x), 3 (-y)

    Returns:
        dnpdata: data object after phase cycle applied
    """

    data, isDict = return_data(all_data)

    data = data.copy()

    if dim not in data.dims:
        raise ValueError("dim not in dims")

    coord = data.coords[dim]
    receiver_phase = _np.array(receiver_phase).ravel()

    proc_parameters = {"dim": dim, "receiver_phase": receiver_phase}

    receiver_phase = _np.tile(receiver_phase, int(coord.size / receiver_phase.size))

    index = data.dims.index(dim)

    reshape_size = [1 for k in data.dims]
    reshape_size[index] = len(data.coords[dim])

    data *= _np.exp(1j * (_np.pi / 2.0) * receiver_phase.reshape(reshape_size))

    proc_attr_name = "phasecycle"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data
