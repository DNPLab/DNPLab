""" Modules for data processing when using CPMG acquisition"""

import numpy as _np
import matplotlib.pyplot as _plt
from scipy.signal import savgol_filter as _savgol_filter
from numpy import random as _random

from ..config.config import DNPLAB_CONFIG
from .integration import integrate


dark_green = DNPLAB_CONFIG.get("COLORS", "dark_green")
light_green = DNPLAB_CONFIG.get("COLORS", "light_green")
dark_grey = DNPLAB_CONFIG.get("COLORS", "dark_grey")
light_grey = DNPLAB_CONFIG.get("COLORS", "light_grey")
orange = DNPLAB_CONFIG.get("COLORS", "orange")


def cpmg_detect_first_echo(data, region=[0, -1], graphical_output=True):
    """Find time for first echo

    Args:
        data (dnpdata):                 CPMG tansient (1D)
        region (list):                  List with start and end time to look for first echo
        graphical_output (boolean):     Display graphical output

    Return:
        t_echo:                         Time of first echo.
    """

    t = data.coords["t2"]
    signal = data.values
    signal = _np.squeeze(signal)
    dwell_time = _np.diff(t)[0]

    region_start_index = _np.squeeze(_np.where(_np.isclose(t, region[0])))
    region_stop_index = _np.squeeze(_np.where(_np.isclose(t, region[1])))

    echo_index = _np.argmax(signal[region_start_index:region_stop_index], axis=0)

    t_echo = _np.round(t[echo_index + region_start_index], 1)

    if graphical_output == True:

        _plt.subplot(2, 1, 1)
        _plt.plot(t, signal, color=light_grey, label="Transient")
        _plt.plot(
            t[region_start_index:region_stop_index],
            signal[region_start_index:region_stop_index],
            color=dark_green,
            label="Search Region",
        )
        _plt.plot(
            [t_echo, t_echo],
            [_np.min(signal), _np.max(signal)],
            color=orange,
            linestyle="--",
        )
        _plt.title("Transient. First Echo at: " + str(t_echo) + " ns")
        _plt.xlabel("Time (ns)")
        _plt.ylabel("Signal (a.u.)")
        _plt.grid(True)
        _plt.legend()

        _plt.subplot(2, 2, 3)
        _plt.plot(
            t[region_start_index:region_stop_index],
            signal[region_start_index:region_stop_index],
            color=dark_green,
        )
        _plt.plot(
            [t_echo, t_echo],
            [
                _np.min(signal[region_start_index:region_stop_index]),
                _np.max(signal[region_start_index:region_stop_index]),
            ],
            color=orange,
            linestyle="--",
        )
        _plt.grid(True)

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
    graphical_output=True,
    verbose=False,
):
    """
    Display integration regions and calculate integration regions for CPMG trace

    data needs to be dnpdata object, but needs to be 1D

    Return a list of tuples with regions

    """

    # Need to pull this out of the class for now. It's hacky but we can clean this up later
    t = data.coords["t2"]
    signal = data.values
    signal = _np.squeeze(signal)

    filter_window_length = 20
    filter_polynomial_order = 5

    t_echo = _np.linspace(t_start + 0, t_start + (n_echo - 1) * t_period, n_echo)
    t_int_start = t_echo - t_width / 2
    t_int_end = t_echo + t_width / 2

    dwell_time = _np.diff(t)[0]

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

    if graphical_output == True:

        # Plot results
        _plt.subplot(2, 1, 1)
        _plt.plot(data.coords["t2"], data.values, color=dark_green)
        _plt.title("CPMG Transient")
        _plt.xlabel("Time (ns)")
        _plt.ylabel("Signal (a.u.)")
        _plt.grid(True)

        min_signal = _np.min(data.values.real)
        max_signal = _np.max(data.values.real)

        for index, value in enumerate(t_echo):

            _plt.plot(
                [value, value], [min_signal, max_signal], color=orange, linestyle="--"
            )
            _plt.plot(
                [t_int_start[index], t_int_start[index]],
                [min_signal, max_signal],
                color=light_grey,
                linestyle="--",
            )
            _plt.plot(
                [t_int_end[index], t_int_end[index]],
                [min_signal, max_signal],
                color=light_grey,
                linestyle="--",
            )

        _plt.plot([data.coords["t2"][noise_region[0]], data.coords["t2"][noise_region[0]]], [0.1 * min_signal, 0.1 * max_signal], color=light_green, linestyle="--")
        _plt.plot([data.coords["t2"][noise_region[0]], data.coords["t2"][noise_region[1]]], [0.1 * max_signal, 0.1 * max_signal], color=light_green, linestyle="--")
        _plt.plot([data.coords["t2"][noise_region[1]], data.coords["t2"][noise_region[1]]], [0.1 * min_signal, 0.1 * max_signal], color=light_green, linestyle="--")
        _plt.plot([data.coords["t2"][noise_region[0]], data.coords["t2"][noise_region[1]]], [0.1 * min_signal, 0.1 * min_signal], color=light_green, linestyle="--")

        _plt.autoscale(enable=True, axis="x", tight=True)

        _plt.subplot(2, 2, 3)
        _plt.plot(t_int_window, echo_slice.real)
        _plt.title("Echo Slices")
        _plt.xlabel("Time (ns)")
        _plt.ylabel("Signal (a.u.)")
        _plt.grid(True)

        # Plot Results
        _plt.subplot(2, 2, 4)
        _plt.plot(_np.linspace(1, n_echo, n_echo), snr.real, ".-")
        # _plt.plot(echo_slice_smoothed[0])
        _plt.title("Relative SNR Increase")
        _plt.xlabel("Number of Echoes")
        _plt.ylabel("SNR")
        _plt.grid(True)

        _plt.tight_layout()
        _plt.show()

    return regions


def cpmg_integrate(data, regions, dim="f2", alternate=True):
    """

    This first uses the DNPLab integrate function to integrate the data
    and then changes the sign depending on whether the echo sign needs
    to be alternated.






    Integrate data along given dimension. If no region is given, the integral is calculated over the entire range.

    Args:
        data (DNPData): Data object
        dim (str): Dimension to perform integration. Default is "f2"
        regions (None, list): List of tuples defining the region to integrate

    Returns:
        data (DNPData): Integrals of data. If multiple regions are given the first value corresponds to the first region, the second value corresponds to the second region, etc.

    Examples:
        Integrated entire data region:

            >>> data = dnp.integrate(data)

        Integrate single peak/region:

            >>> data = dnp.integrate(data, regions=[(4, 5)])

        Integrate two regions:

            >>> data = dnp.integrate(data, regions=[(1.1, 2.1), (4.5, 4.9)])

    """
    out = data.copy()

    out = integrate(out, regions=regions, dim=dim)

    if alternate == True:
        sign = (-1) ** out.coords[1]

        # This next operation needs to be checked. I don't even understand why
        # this is working

        out.values = out.values * sign

    return out
