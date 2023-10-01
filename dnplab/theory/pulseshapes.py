"""Functions to calculate AWG pulse shapes"""

import dnplab as _dnp
import numpy as _np


def wurst(tp, N, resolution=1.0e-9):
    """Calculate real value WURST envelope pulse shape

    .. math::
        1 - \\text{abs} \left( \cos \left( \\frac{\pi}{t_p} (t - \\frac{t_p}{2}) + \\frac{\pi}{2} \\right) \\right) ^N

    Args:
        tp (float): Pulse length
        N (float): exponential
        resolution: Resolution (time increment). Default is 1 ns

    Returns:
        tuple: tuple containing:

        t (*numpy.ndarray*): Time axis
        pulse (*numpy.ndarray*): Real pulse shape of WURST pulse
    """

    t = _np.r_[0.0:tp:resolution]
    pulse = (
        1.0 - _np.abs(_np.cos(_np.pi * (t - tp / 2.0) / tp + _np.pi / 2.0)) ** N
    ) + 0j

    return t, pulse
