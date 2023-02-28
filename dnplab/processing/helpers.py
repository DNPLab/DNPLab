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


def signal_to_noise(
    data: DNPData,
    signal_region: tuple,
    noise_region: tuple,
    dim: str = "f2",
    method="absval",
    detrend: bool = True,
    fullOutput: bool = False,
):
    """Find signal-to-noise ratio

    Simplest implementation: select largest value in a signal_region and divide this value by the estimated std. of the baseline in another region

    The signature of this function is not final and is subject to changes - do use at your own risk!

    Args:
        data: Spectrum data
        signal_region : tuple (start,stop) of region where signal should be searched
        noise_region : tuple (start,stop) of region that should be taken as noise
        dim (default f2): dimension of data that is used for snr
        detrend:bool=False is whether a quadratic detrending should be done on the
        fullOutput:bool=False whether singal and noise should also be returned

        method:str=absval currently not used
    Returns:
        NotImplemented
    """
    import warnings

    warnings.warn(
        "helpers.signal_to_noise: The signature and implementation of this function is not final and is subject to changes - do use at your own risk!"
    )

    def signal_estimation(data, signal_region: tuple, dim: str = "f2"):
        signal = _np.max(_np.abs(data[dim, signal_region[0] : signal_region[1]]))
        return signal

    def noise_estimation(
        data, noise_region: tuple, dim: str = "f2", detrend: bool = detrend
    ):
        noise_data = _np.abs(
            data[dim, noise_region[0] : noise_region[1]]
        )  # make it behave as numpy array?
        f_index = noise_data.index(dim)
        f_noise = data.coords[f_index][noise_region[0] : noise_region[1]]
        if detrend:
            # detrend by making quadratic fit
            ind_mean = int(noise_data.size / 2)
            if noise_data.size < 3:
                raise ValueError(
                    "Please chose larger region, noise data region is too small ({0}) for quadratic detrending fit".format(
                        noise_region
                    )
                )
            ind_last = noise_data.size - 1
            guess_data = _np.array(
                (noise_data[dim, 0], noise_data[dim, ind_mean], noise_data[dim, -1])
            )
            # work on indices not on frequencies
            ind = _np.arange(0, noise_data.size)
            guess_mat = _np.array(
                [[0, 0, 1], [ind_mean**2, ind_mean, 1], [ind_last**2, ind_last, 1]]
            )
            init_guess = _np.squeeze(_np.linalg.solve(guess_mat, guess_data))
            detrend_fun = (
                lambda prm, x, data: prm[0] * x**2 + prm[1] * x + prm[2] - data
            )
            result = _scipy_optimize.least_squares(
                detrend_fun, init_guess, args=(ind, noise_data.values)
            )
            if not result.success:
                warnigns.warn(
                    "No success with detrending quadratic fit, do not consider the result as meaningful! maybe try with detrend=False"
                )
            # detrended, only mean value remains, should work with DNPdata
            noise_data = (
                _np.mean(noise_data) + noise_data - detrend_fun(result.x, ind, 0)
            )
        noise = _np.std(noise_data)
        return noise

    signal = signal_estimation(data, signal_region, dim)
    noise = noise_estimation(data, noise_region, dim, detrend)

    if fullOutput:
        return signal / noise, signal, noise
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
