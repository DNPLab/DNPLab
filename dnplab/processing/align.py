import numpy as _np
from ..core.data import DNPData
import warnings


def ndalign(data, dim="f2", reference=None, center=None, width=None, average=None):
    """Alignment of NMR spectra using FFT Cross Correlation

    Args:
        all_data (object) : DNPData object
        dim (str) : Dimension to align along
        reference (numpy) : Reference spectra for alignment
        center (float) : Center of alignment range, by default entire range
        width (float) : Width of alignment range, by default entire range
        average (int) : Number of averages to use for reference

    Returns:
        DNPData: Aligned data

    Examples:

        >>> data_aligned = dnp.ndalign(data)
        >>> data_aligned = dnp.ndalign(data, center = 10, width = 20)
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
        raise ValueError("selected range is not acceptable")

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
    elif isinstance(reference, DNPData):
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
        if ix == 0:
            first_shift = delta_max_ix
        delta_max_ix -= first_shift
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

    return data
