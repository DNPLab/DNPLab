import numpy as _np
from ..core.data import DNPData
import warnings


def ndalign(data, dim="f2", reference=None, center=None, width=None):
    """Alignment of NMR spectra using FFT Cross Correlation

    Args:
        all_data (object) : DNPData object
        dim (str) : Dimension to align along
        reference (numpy) : Reference spectra for alignment
        center (float) : Center of alignment range, by default entire range
        width (float) : Width of alignment range, by default entire range

    Returns:
        DNPData: Aligned data

    Examples:

        >>> data_aligned = dnp.ndalign(data)
        >>> data_aligned = dnp.ndalign(data, center = 10, width = 20)
    """

    out = data.copy()

    proc_parameters = {"dim": dim}

    if center != None and width != None:
        start = center - 0.5 * width
        stop = center + 0.5 * width
    elif center == None and width == None:
        start = out.coords[dim][-1]
        stop = out.coords[dim][0]
    else:
        raise ValueError("selected range is not acceptable")

    temp_out = out[dim, (start, stop)].copy()
    temp_out.unfold(dim)
    temp_values = temp_out.values.T

    out.unfold(dim)
    all_values = out.values.T

    abs_temp_values = _np.abs(temp_values)

    if reference is None:
        reference = _np.abs(temp_values[-1])
    elif isinstance(reference, DNPData):
        reference = _np.abs(reference.values)

    ref_max_ix = _np.argmax(reference)

    aligned_all_values = _np.zeros_like(all_values)

    for ix in range(len(abs_temp_values)):
        cor = _np.correlate(
            abs_temp_values[ix], reference, mode="same"
        )  # calculate cross-correlation
        max_ix = _np.argmax(cor)  # Maximum of cross correlation
        delta_max_ix = max_ix - ref_max_ix  # Calculate how many points to shift
        if ix == 0:
            first_shift = delta_max_ix
        delta_max_ix -= first_shift
        aligned_all_values[ix] = _np.roll(
            all_values[ix], -1 * delta_max_ix
        )  # shift values

    out.values = aligned_all_values.T  # Add aligned values back to data object

    out.fold()  # Back to original order

    proc_attr_name = "ndalign"
    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out
