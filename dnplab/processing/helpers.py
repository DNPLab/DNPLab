import numpy as _np
from scipy.signal import savgol_filter

from ..core.data import DNPData
from ..processing.integration import integrate

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


def signal_to_noise():
    """Find signal-to-noise ratio

    Returns:
        NotImplemented
    """

    return NotImplemented


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


def normalize():
    return NotImplemented


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
