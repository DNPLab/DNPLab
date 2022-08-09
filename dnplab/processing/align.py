import numpy as np
from ..core.data import DNPData


def ndalign(data, dim="f2", reference=None, center=None, width=None, average=None):
    """Alignment of NMR spectra using FFT Cross Correlation

    Args:
        all_data (object) : dnpdata object
        dim (str) : dimension to align along
        reference (numpy) : second dimension to align along
        center (float) : range center
        width (float) : range width

    Returns:
        dnpdata: Aligned data in container
    
    Examples:

        >>> import numpy as np
        >>> from matplotlib.pylab import *
        >>> import dnplab as dnp
        >>> x = np.r_[-10:10:100j]
        >>> y = np.r_[0,1,2]
        >>> values = dnp.math.lineshape.lorentzian(x.reshape(-1,1),y.reshape(1,-1),1)
        >>> data = dnp.DNPData(values, ['f2','sample'], [x, y])
        >>> figure('Before Alignment')
        >>> dnp.plot(data)
        >>> data_aligned = dnp.ndalign(data)
        >>> figure('After Alignment')
        >>> dnp.plot(data_aligned)
        >>> show()
    """

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
        raise ValueError("selected range is not accpetable")

    values = data[dim, (start, stop)].values

    all_original_shape = all_values.shape
    original_shape = values.shape  # Preserve original shape

    all_align_dim_length = all_original_shape[0]
    align_dim_length = original_shape[0]  # length of dimension to align down

    all_values = all_values.reshape(all_align_dim_length, -1)
    values = values.reshape(align_dim_length, -1)  # Reshape to 2d

    new_shape = np.shape(values)

    dim2 = new_shape[1]

    abs_values = np.abs(values)

    if reference is None:
        reference = np.abs(values[:, -1])
    elif isinstance(reference, DNPData):
        reference = np.abs(reference.values)
        if average != None:
            reference = np.convolve(reference, np.ones(average), "same") / average

    ref_max_ix = np.argmax(reference)

    all_aligned_values = np.zeros_like(all_values)

    for ix in range(dim2):
        if average != None:
            abs_values[:, ix] = (
                np.convolve(abs_values[:, ix], np.ones(average), "same") / average
            )
        cor = np.correlate(
            abs_values[:, ix], reference, mode="same"
        )  # calculate cross-correlation
        max_ix = np.argmax(cor)  # Maximum of cross correlation
        delta_max_ix = max_ix - ref_max_ix  # Calculate how many points to shift
        all_aligned_values[:, ix] = np.roll(
            all_values[:, ix], -1 * delta_max_ix
        )  # shift values

    all_aligned_values = all_aligned_values.reshape(
        all_original_shape
    )  # reshape to original values shape

    data.values = all_aligned_values  # Add aligned values back to data object

    data.reorder(original_order)  # Back to original order

    proc_attr_name = "ndalign"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data


def align(data, dim="f2", dim2=None, center=None, width=None):
    """DEPRECIATED. Please use ndalign instead. Alignment of NMR spectra down given dimension or dimensions

    Args:
        all_data (object) : dnpdata object
        dim (str) : dimension to align along
        dim2 (str) : second dimension to align along
        center (float) : range center
        width (float) : range width

    Returns:
        dnpdata: Aligned data in container
    """

    if len(np.shape(data.values)) > 3:
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

            corrData = np.correlate(np.abs(rangeData), np.abs(refData), mode="same")
            shiftIx = np.argmax(corrData) - (
                len(corrData) / 2
            )  # subtract half length so spectrum is shifted relative to center, not edge
            shiftData = np.roll(tempData, -1 * int(np.round(shiftIx, 0)))
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

                corrData = np.correlate(np.abs(rangeData), np.abs(refData), mode="same")
                shiftIx = np.argmax(corrData) - (
                    len(corrData) / 2
                )  # subtract half length so spectrum is shifted relative to center, not edge
                shiftData = np.roll(tempData, -1 * int(np.round(shiftIx, 0)))
                data.values[:, ix2, ix1] = shiftData

    data.reorder(originalAxesOrder)

    proc_attr_name = "align"
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data
