import numpy as _np
from scipy.signal import savgol_filter

from ..core.data import DNPData
from ..processing.integration import integrate
from ..processing.offset import remove_background as dnp_remove_background

import dnplab as dnp


def calculate_enhancement(data, off_spectrum_index=0, return_complex_values=False):
    """Calculate enhancement of a power series. Needs integrals as input

    Args:
        integrals (DNPData):
        off_spectrum_index (int):
        return_complex_values (bool):

    Returns:
        enhancements (DNPData): Enhancement values
    """

    enhancements = data.copy()

    if not "experiment_type" in data.attrs.keys():
        raise KeyError("Experiment type not defined")

    if data.attrs["experiment_type"] != "integrals":
        raise ValueError("dnpdata object does not contain integrals.")

    if data.dims[0] == "Power":
        enhancements.attrs["experiment_type"] = "enhancements_P"

        enhancements.values = (
            enhancements.values / enhancements.values[off_spectrum_index]
        )

    elif data.dims[0] == "B0":
        enhancements.attrs["experiment_type"] = "enhancements_B0"
        print("This is a DNP enhancement profile. Not implemented yet.")

    else:
        raise TypeError(
            "Integration axis not recognized. First dimension should be Power or B0."
        )

    if return_complex_values == True:
        return enhancements

    elif return_complex_values == False:
        return enhancements.real


def signal_to_noise(
    data: DNPData,
    signal_region: list = slice(0, None),
    noise_region: list = (None, None),
    dim: str = "f2",
    remove_background: list = None,
    **kwargs
):
    """Find signal-to-noise ratio

    Simplest implementation: select largest value in a signal_region and divide this value by the estimated std. deviation of another noise_region. If the noise_region list contains (None,None) (the default) then all points except the points +10% and -10% around the maximum are used for the noise_region.

    Args:
        data: Spectrum data
        signal_region (list): list with a single tuple (start,stop) of a region where a signal should be searched, default is [slice(0,None)] which is the whole spectrum
        noise_region (list): list with tuples (start,stop) of regions that should be taken as noise, default is (None,None)
        dim (str): dimension of data that is used for snr calculation, default is 'f2'
        remove_background (list): if this is not None (a list of tuples, or a single tuple) this will be forwarded to dnp.remove_background, together with any kwargs
        kwargs : parameters for dnp.remove_background

    Returns:
        SNR (float): Signal to noise ratio

    Examples:

        A note for the usage: regions can be provided as (min,max), slices use indices.
        To use the standard values just use

            >>> snr = dnp.signal_to_noise(data)

        If you want to select a region for the noise and the signal:

            >>> snr = dnp.signal_to_noise(data,[(-1.23,300.4)],noise_region=[(-400,-240.5),(123.4,213.5)])

        With background subtracted:

            >>> snr = dnp.signal_to_noise(data,[(-1.23,300.4)],noise_region=[(-400,-240.5),(123.4,213.5)],remove_background=[(123.4,213.5)])

        This function allows to use a single tuple instead of a list with a single tuple for signal_region, noise_region and remove_background. This is for convenience, slices are currently only supoprted for signal_region and noise_region.

            >>> snr = dnp.signal_to_noise(data,(-1.23,300.4),noise_region=[(-400,-240.5),(123.4,213.5],remove_background=(123.4,213.5))

    """
    import warnings
    import scipy.optimize as _scipy_optimize

    # convenience for signal and noise region
    def _convenience_tuple_to_list(possible_region: list):
        if possible_region is None:
            return possible_region
        # we assume its iterable
        try:
            l = len(possible_region)
            if l != 2:
                return possible_region  # it seems to be iterable?
        except TypeError:
            return [possible_region]  # return as is in a list, might be slice
        try:
            # check whether we can interpret it as value
            a = int(possible_region[0])
            return [
                (possible_region[0], possible_region[1])
            ]  # make a list that contains a tuple
        except TypeError:
            if type(possible_region) == list:
                return possible_region
            return [possible_region]

    signal_region = _convenience_tuple_to_list(signal_region)
    if len(signal_region) > 1:
        raise ValueError(
            "More than one signal region ({0}) provided in signal_to_noise. This is not supported.".format(
                signal_region
            )
        )
    noise_region = _convenience_tuple_to_list(noise_region)
    remove_background = _convenience_tuple_to_list(remove_background)

    if not (dim in data.dims):
        raise ValueError("dim {0} not in data.dims ({1})".format(dim, data.dims))

    # remove background
    if remove_background is not None:
        deg = kwargs.pop("deg", 1)
        data = dnp_remove_background(data, dim, deg, remove_background)

    # currently only one method avaiable -> absolute value
    signal = _np.max(_np.abs(data[dim, signal_region[0]]))

    if (None, None) in noise_region:
        signal_arg = _np.argmax(_np.abs(data[dim, signal_region[0]]))
        datasize = data[dim, :].size
        noise_region = [
            slice(0, int(_np.maximum(2, int(signal_arg * 0.9)))),
            slice(int(_np.minimum(datasize - 2, int(signal_arg * 1.1))), None),
        ]

    # concatenate noise_regions
    noise_0 = _np.abs(data[dim, noise_region[0]])

    for k in noise_region[1:]:
        noise_0.concatenate(_np.abs(data[dim, k]), dim)

    noise = _np.std(_np.abs(noise_0[dim, slice(0, None)]))

    return signal / noise


def smooth(data, dim="t2", window_length=11, polyorder=3):
    """Apply Savitzky-Golay Smoothing

    This function is a wrapper function for the savgol_filter from the SciPy python package (https://scipy.org/). For a more detailed description see the SciPy help for this function.

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform smoothing
        window_length (int): Length of window (number of coefficients)
        polyorder (int): Polynomial order to fit samples

    Returns:
        data (DNPData): Data with Savitzky-Golay smoothing applied
    """
    out = data.copy()

    out.unfold(dim)

    out.values = savgol_filter(out.values, window_length, polyorder, axis=0)

    out.fold()

    return out


def left_shift(data, dim="t2", shift_points=0):
    """Remove points from the left

    Args:
        data (DNPData): Data object
        dim (str): Name of dimension to left shift, default is "t2"
        shift_points (int): Number of points to left shift, default is 0.

    Returns:
        data (DNPDdata): Shifted data object
    """

    data = data[dim, shift_points:]

    proc_attr_name = "left_shift"
    proc_parameters = {
        "dim": dim,
        "points": shift_points,
    }
    data.add_proc_attrs(proc_attr_name, proc_parameters)

    return data


def normalize(data, amplitude=True):
    """Normalize spectrum

    Args:
        data (DNPData): Data object
        amplitude (boolean): True: normalize amplitude, false: normalize area. The default is True

    Returns:
        data (DNPDdata): Normalized data object
    """

    out = data.copy()

    if amplitude == True:
        out.values = out.values / _np.max(out.values)
    elif amplitude == False:
        out.values = out.values  # Normalize to area = 1, not implemented yet

    proc_attr_name = "normalized"
    proc_parameters = {
        "amplitude": amplitude,
    }

    out.add_proc_attrs(proc_attr_name, proc_parameters)

    return out


def reference(data, dim="f2", old_ref=0, new_ref=0):
    """Function for referencing NMR spectra

    Args:
        data (DNPData): Data for referencing
        dim (str): dimension to perform referencing down. By default this dimension is "f2".
        old_ref (float): Value of old reference
        new_ref (float): New reference value

    Returns:
        DNPData: referenced data
    """

    data = data.copy()

    data.coords[dim] -= old_ref - new_ref

    return data
