""" Modules for data processing when using CPMG acquisition"""

import numpy as _np
import matplotlib.pyplot as _plt
import dnplab as _dnp
from scipy.signal import savgol_filter as _savgol_filter
from numpy import random as _random

from .integration import integrate

def cpmg_detect_first_echo(
    data, region=[0, -1], graphical_output=True, use_real=True, verbose=False
):
    """Find time of first echo in CPMG sequence

    Args:
        data (dnpdata):                 CPMG tansient (1D)
        region (list):                  List with start and end time to look for first echo
        graphical_output (boolean):     Display graphical output
        use_real (boolean):             Omit imaginary part when using complex data. Default is true
        verbose (boolean):              Enable additional output for debugging

    Return:
        t_echo (float):                 Time of first echo.
    """
    if verbose == True:
        print("********** This is: cpmg_detect_first_echo **********")
        print("Input data file: ")
        print(data)

    if _np.squeeze(_np.shape(data.dims)) > 2:
        print("Dataset with dimensions > 2 currently not supported")
        return

    if _np.shape(data.values)[1] == 1:
        n_dim = 1
        plot_title_pre_text = "CPMG Transient"
    elif _np.shape(data.values)[1] > 1:
        n_dim = 2
        plot_title_pre_text = "CPMG Transient (Integrated)"

    if verbose == True:
        if n_dim == 1:
            print("Data set is 1D.")
        if n_dim == 2:
            print(
                "Data set is 2D. For processing the data set is reduced to 1D by integrating along the second axis."
            )
            print("The second axis is: ", data.dims[1])

    if n_dim == 2:
        data = integrate(data, dim=data.dims[1])

    t = data.coords["t2"]

    if use_real == True:
        signal = data.values.real
        if verbose == True:
            print("Use only real part of data set")
    else:
        signal = data.values

    signal = _np.squeeze(signal)
    dwell_time = _np.diff(t)[0]

    region_start_index = _np.squeeze(_np.where(_np.isclose(t, region[0])))
    region_stop_index = _np.squeeze(_np.where(_np.isclose(t, region[1])))

    echo_index = _np.argmax(signal[region_start_index:region_stop_index], axis=0)

    t_echo = _np.round(t[echo_index + region_start_index], 1)

    if verbose == True:
        print("First Echo found at: " + str(t_echo) + " ns")

    if graphical_output == True:

        # This is a bit hacky, but the standard font size is too large when using subplots
        plot_font_size = 12
        _plt.rcParams["font.size"] = plot_font_size

        _plt.subplot(2, 1, 1)
        _plt.plot(t, signal, color=_dnp.colors.secondary2, label="Transient")
        _plt.plot(
            t[region_start_index:region_stop_index],
            signal[region_start_index:region_stop_index],
            color=_dnp.colors.primary1,
            label="Search Region",
        )
        _plt.plot(
            [t_echo, t_echo],
            [_np.min(signal), _np.max(signal)],
            color=_dnp.colors.accent,
            linestyle="--",
        )
        _plt.title(plot_title_pre_text + ". 1st Echo at: " + str(t_echo) + " ns", fontsize = plot_font_size)
        _plt.xlabel("Time (ns)")
        _plt.ylabel("Signal (a.u.)")
        _plt.grid(True)
        # _plt.legend(loc="lower right")

        _plt.subplot(2, 2, 3)
        _plt.plot(
            t[region_start_index:region_stop_index],
            signal[region_start_index:region_stop_index],
            color=_dnp.colors.primary1,
        )
        _plt.plot(
            [t_echo, t_echo],
            [
                _np.min(signal[region_start_index:region_stop_index]),
                _np.max(signal[region_start_index:region_stop_index]),
            ],
            color=_dnp.colors.accent,
            linestyle="--",
        )
        _plt.grid(True)
        _plt.xlabel("Time (ns)")
        _plt.ylabel("Signal (a.u.)")
        _plt.title("1st Echo", fontsize = plot_font_size)


        _plt.tight_layout()
        _plt.show()

    return t_echo


def cpmg_show_integration_region(
    data,
    n_echo,
    t_start,
    t_period,
    t_width,
    noise_region=[None, None],
    alternate=True,
    use_real=True,
    graphical_output=True,
    verbose=False,
):
    """Process data using CPMG acquisition

    Display integration regions and calculate integration regions for CPMG trace.

    Known issues:
        * Data has to be 1D. Need to come up with a good strategy how to deal with higher dimensions
        * Noise region only accepts index, not float values.

    Args:
        data (DNPdata):             Data object containing transient of CPMG detected signals
        n_echo (int):               Number of echoes to use for processing
        t_start (float):            Position of first echo
        t_period (float):           Period of echoes
        t_width (float):            Integration width
        noise_region (int/float):   Region to use to calculate noise. Currently requires index. Default is [None, None]. Only one noise region will be used (only one tuple)
        alternate=True (boolean):   Alternate sign of echo amplitude. Default is true
        use_real (boolean):         Omit imaginary part when using complex data. Default is true
        graphical_output (boolean): Show graphical output. Default is true
        verbose=False (boolean):    Enable additional output for debugging.

    Return:
        data (DNPData):             DNPData object containing integration regions. These regions are then used in cpmg_integrate() to process data.

    """
    if verbose == True:
        print("********** This is: cpmg_show_integration_region **********")
        print("Input data file: ")
        print(data)

    if _np.squeeze(_np.shape(data.dims)) > 2:
        print("Dataset with dimensions > 2 currently not supported")
        return

    if _np.shape(data.values)[1] == 1:
        n_dim = 1
        plot_title_pre_text = "CPMG Transient"
    elif _np.shape(data.values)[1] > 1:
        n_dim = 2
        plot_title_pre_text = "CPMG Transient (Integrated)"

    if verbose == True:
        if n_dim == 1:
            print("Data set is 1D.")
        if n_dim == 2:
            print(
                "Data set is 2D. For processing the data set is reduced to 1D by integrating along the second axis."
            )
            print("The second axis is: ", data.dims[1])

    if n_dim == 2:
        data = integrate(data, dim=data.dims[1])

    # Need to pull this out of the class for now. It's hacky but we can clean this up later
    t = data.coords["t2"]

    if use_real == True:
        signal = data.values.real
        if verbose == True:
            print("Use only real part of data set")
    else:
        signal = data.values

    signal = _np.squeeze(signal)

    dwell_time = _np.diff(t)[0]

    # Make sure all values are multiple of the dwell time
    t_period = int(t_period / dwell_time) * dwell_time
    t_start = int(t_start / dwell_time) * dwell_time
    t_width = int(t_width / dwell_time) * dwell_time

    t_echo = _np.linspace(t_start + 0, t_start + (n_echo - 1) * t_period, n_echo)
    t_int_start = t_echo - t_width / 2
    t_int_end = t_echo + t_width / 2

    t_width_pts = int(_np.round(t_width / dwell_time))
    t_int_window = _np.linspace(
        -t_width_pts / 2 * dwell_time, t_width_pts / 2 * dwell_time, t_width_pts
    )

    echo_slice = _np.empty([t_width_pts, n_echo], dtype=complex)

    regions = []

    # Slice up transient and create sub-array for each echo slice
    for index, value in enumerate(t_echo):

        index_echo_max = _np.squeeze(
            _np.where(_np.isclose(t, value, atol=dwell_time / 2))[0]
        )
        index_int_start = int(index_echo_max - t_width_pts / 2)
        index_int_stop = int(index_echo_max + t_width_pts / 2)

        regions.append((index_int_start * dwell_time, index_int_stop * dwell_time))

        if alternate == True:
            echo_slice[:, index] = (-1) ** index * signal[
                index_int_start:index_int_stop
            ]
        else:
            echo_slice[:, index] = signal[index_int_start:index_int_stop]

    # Calculate signal to noise ratios
    if noise_region[0] == None:
        noise_region[0] = 0

    if noise_region[1] == None:
        noise_region[1] = -1

    noise = signal[noise_region[0] : noise_region[1]]
    noise_std = _np.std(noise)

    if verbose == True:
        print("Noise value (std): ", noise_std)

    signal = _np.sum(_np.cumsum(echo_slice, axis=0), axis=0) / len(echo_slice[:, 0])
    signal_cumsum = _np.cumsum(signal)

    snr = signal_cumsum.real / noise_std

    n = _np.linspace(1, n_echo, n_echo)
    snr = (1 / _np.sqrt(n)) * snr / snr[0]

    if verbose == True:
        print("Maximum SNR increase: " + str(_np.round(_np.max(snr), 2)))

    # Plot results
    if graphical_output == True:

        # This is a bit hacky, but the standard font size is too large when using subplots
        plot_font_size = 12
        _plt.rcParams["font.size"] = plot_font_size

        _plt.subplot(2, 1, 1)
        _plt.plot(data.coords["t2"].real, data.values.real, color=_dnp.colors.primary1)
        _plt.title(plot_title_pre_text, fontsize = plot_font_size)
        _plt.xlabel("Time (ns)")
        _plt.ylabel("Signal (a.u.)")
        _plt.grid(True)

        min_signal = _np.min(data.values.real)
        max_signal = _np.max(data.values.real)

        for index, value in enumerate(t_echo):

            _plt.plot(
                [value, value],
                [min_signal, max_signal],
                color=_dnp.colors.accent,
                linestyle="--",
            )
            _plt.plot(
                [t_int_start[index], t_int_start[index]],
                [min_signal, max_signal],
                color=_dnp.colors.secondary2,
                linestyle="--",
            )
            _plt.plot(
                [t_int_end[index], t_int_end[index]],
                [min_signal, max_signal],
                color=_dnp.colors.secondary2,
                linestyle="--",
            )

        _plt.plot(
            [data.coords["t2"][noise_region[0]], data.coords["t2"][noise_region[0]]],
            [0.1 * min_signal, 0.1 * max_signal],
            color=_dnp.colors.primary2,
            linestyle="--",
        )
        _plt.plot(
            [data.coords["t2"][noise_region[0]], data.coords["t2"][noise_region[1]]],
            [0.1 * max_signal, 0.1 * max_signal],
            color=_dnp.colors.primary2,
            linestyle="--",
        )
        _plt.plot(
            [data.coords["t2"][noise_region[1]], data.coords["t2"][noise_region[1]]],
            [0.1 * min_signal, 0.1 * max_signal],
            color=_dnp.colors.primary2,
            linestyle="--",
        )
        _plt.plot(
            [data.coords["t2"][noise_region[0]], data.coords["t2"][noise_region[1]]],
            [0.1 * min_signal, 0.1 * min_signal],
            color=_dnp.colors.primary2,
            linestyle="--",
        )

        _plt.autoscale(enable=True, axis="x", tight=True)

        _plt.subplot(2, 2, 3)
        _plt.plot(t_int_window, echo_slice.real)
        _plt.title("Echo Slices", fontsize = plot_font_size)
        _plt.xlabel("Time (ns)")
        _plt.ylabel("Signal (a.u.)")
        _plt.grid(True)

        # Plot Results
        _plt.subplot(2, 2, 4)
        _plt.plot(_np.linspace(1, n_echo, n_echo), snr.real, ".-")
        # _plt.plot(echo_slice_smoothed[0])
        _plt.title("Relative SNR Increase", fontsize = plot_font_size)
        _plt.xlabel("Number of Echoes")
        _plt.ylabel("SNR")
        _plt.grid(True)

        _plt.tight_layout()
        _plt.show()

    return regions


def cpmg_integrate(data, regions, dim="t2", alternate=True):
    """Process data using CPMG acquisition

    First use cpmg_detect_first_echo() to get position of first echo.
    Then use cpmg_show_integration_region() to determine the integration regions.

    Integrate data along given dimension. If no region is given, the integral is calculated over the entire range.

    The function is a wrapper for the DNPlab integrate() function

    Args:
        data (DNPData):         Data object
        regions (None, list):   List of tuples defining the region to integrate
        dim (str):              Dimension to perform integration. Default is "t2"
        alternate (boolean):    Alternate the sign of the echo traces. Default is true
        normalize (boolean):    Normalize trace of integration. Default is true

    Returns:
        data (DNPData):         Integrals of data. If multiple regions are given the first value corresponds to the first region, the second value corresponds to the second region, etc.

    """
    out = data.copy()

    out = integrate(out, regions=regions, dim=dim)

    if alternate == True:
        sign = (-1) ** out.coords[1]

        # This next operation needs to be checked. I don't even understand why
        # this is working? How do we know it is multiplying the correct dimension?
        out.values = out.values * sign

    return out
