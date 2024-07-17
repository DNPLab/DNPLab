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


def cpmg_show_integration_region(
    data, n_echo, t_start, t_period, t_width, alternate=True, graphical_output=True
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

    filter_window_length = 10
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

    # This signal to noise calculation does not seem to be correct. This seems to work
    # on the actual field sweep spectrum, but not on the individual echoes. Not sure
    # why this shouldn't work.

    # Calculate signal to noise ratios for 1st echo and summed up echo
    echo_slice_cumsum = _np.cumsum(echo_slice, axis=1)
    echo_max_amplitude = _np.max(echo_slice_cumsum, axis=0)

    # To calculate signal to noise ratio, filter (smooth) signal and subtract filter signal
    # from signal. Calculate std on remaining signal.

    echo_slice_cumsum_smoothed_real = _savgol_filter(
        echo_slice_cumsum.real,
        window_length=filter_window_length,
        polyorder=filter_polynomial_order,
        axis=0,
    )
    echo_slice_cumsum_smoothed_imag = _savgol_filter(
        echo_slice_cumsum.imag,
        window_length=filter_window_length,
        polyorder=filter_polynomial_order,
        axis=0,
    )
    echo_slice_cumsum_smoothed = (
        echo_slice_cumsum_smoothed_real + echo_slice_cumsum_smoothed_imag
    )

    signal_noise = _np.std(echo_slice_cumsum - echo_slice_cumsum_smoothed, axis=0)
    snr = echo_max_amplitude / signal_noise

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

        _plt.autoscale(enable=True, axis="x", tight=True)

        _plt.subplot(2, 2, 3)
        _plt.plot(t_int_window, echo_slice.real)
        _plt.title("Echo Slices")
        _plt.xlabel("Time (ns)")
        _plt.ylabel("Signal (a.u.)")
        _plt.grid(True)

        # Plot Results
        # _plt.subplot(2,2,4)
        # _plt.plot(_np.linspace(1, n_echo, n_echo), snr.real,".-")
        # _plt.title("Signal to Noise")
        # _plt.xlabel("Number of Echoes")
        # _plt.ylabel("SNR")
        # _plt.grid(True)

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
