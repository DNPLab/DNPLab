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
    signal_region: list,
    noise_region: list,
    dim: str = "f2",
    remove_background: list = None,
    fullOutput: bool = False,
    **kwargs
):
    """Find signal-to-noise ratio

    Simplest implementation: select largest value in a signal_region and divide this value by the estimated std. of the baseline in another region

    The signature of this function is not final and is subject to changes - do use at your own risk!

    Args:
        data: Spectrum data
        signal_region : tuple (start,stop) of region where signal should be searched
        noise_region : tuple (start,stop) of region that should be taken as noise
        dim (default f2): dimension of data that is used for snr
        remove_background :listl=None, if this is not none (a list of tuples, or a single tuple) this will be forwarded to dnp.remove_background, tgether with any kwargs
        fullOutput:bool=False whether signal and noise should also be returned

    Returns:
        data: DNPData, snr:float
        or:
        data: DNPData, snr:float, signal:float, noise: float
    """
    import warnings
    import scipy.optimize as _scipy_optimize

    warnings.warn(
        "helpers.signal_to_noise: The signature and implementation of this function is not final and is subject to changes - do use at your own risk!"
    )

    # convenience for signal and noise region
    def _convenience_tuple_to_list(possible_region: list):
        if possible_region is None:
            return possible_region
        # we assume its iterable
        try:
            l = len(possible_region)
            if l != 2:
                return possible_region
        except TypeError:
            return possible_region  # return as is
        try:
            #check whether we can interpret it as value
            a=int(possible_region[0])
            return [(possible_region[0], possible_region[1])]  # make a list that contains a tuple
        except TypeError:
            return possible_region

    signal_region = _convenience_tuple_to_list(signal_region)
    noise_region = _convenience_tuple_to_list(noise_region)
    remove_background = _convenience_tuple_to_list(remove_background)

    print(remove_background)

    # remove background
    if remove_background is not None:
        deg = kwargs.pop("deg", 2)
        data = dnp_remove_background(data, dim, deg, remove_background)

    # currently only one method avaiable -> absolute value
    signal = _np.max(_np.abs(data[dim, signal_region[0]]))
    noise = _np.std(_np.abs(data[dim, noise_region[0]]))

    if fullOutput:
        return data, signal / noise, signal, noise
    return data, signal / noise


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
