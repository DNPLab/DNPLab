"""Functions to calculate AWG pulse shapes"""

import dnplab as _dnp
import numpy as _np


def chirp(tp, BW, resolution=1.0e-9):
    """Calculate complex chirp pulse shape

    .. math::
        e^{i 2 \pi (k/2) (t - t_p/2)^2}

    Args:
        tp (float): Pulse length in (s)
        BW (float): Pulse excitation bandwidth in (Hz)

    Returns:
        tuple: tuple containing:
            t (*numpy.ndarray*): Time axis
            pulse (*numpy.ndarray*): Pulse shape
    """

    k = BW / tp
    t = _np.r_[0.0:tp:resolution]
    pulse = _np.exp(1.0j * 2.0 * _np.pi * ((k / 2.0) * ((t - tp / 2.0) ** 2.0)))
    return t, pulse


def wurst(tp, N, BW, B1, resolution=1.0e-9):
    """Calculate complex WURST pulse shape

    The WURST pulse can be constructed from two parts: 1) a hyperbolic secant, and 2) a frequency chirp.

    .. math::
        1 - \\text{abs} \left( \cos \left( \\frac{\pi}{t_p} (t - \\frac{t_p}{2}) + \\frac{\pi}{2} \\right) \\right) ^N

    Args:
        tp (float): Pulse length
        N (float): Exponent
        BW (float): Bandwidth in (Hz)
        B1 (float): Pulse strength in (Hz)
        resolution: Resolution (time increment) in (s). Default is 1 ns

    Returns:
        tuple: tuple containing:

        t (*numpy.ndarray*): Time axis
        pulse (*numpy.ndarray*): Real pulse shape of WURST pulse
    """

    t = _np.r_[0.0:tp:resolution]
    pulse_real = (
        1.0 - _np.abs(_np.cos(_np.pi * (t - tp / 2.0) / tp + _np.pi / 2.0)) ** N
    ) + 0j

    pulse_imag = (
        1.0 - _np.abs(_np.cos(_np.pi * (t - tp / 2.0) / tp + _np.pi / 2.0)) ** N
    ) + 0j

    # pulse_wurst = pulse_real + 1j * pulse_imag
    pulse_wurst = pulse_real

    t, pulse_chirp = _dnp.chirp(tp, BW, resolution=resolution)

    pulse = B1 * pulse_wurst * pulse_chirp

    return t, pulse
