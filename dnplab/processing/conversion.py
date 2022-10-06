import numpy as _np
from ..constants.constants import *


def dBm2w(power_in_dBm):
    """Convert power in dBm to power in W

    Convert a microwave power given in dBm to W

    Args:
        power_in_dBm:   Power in (dBm)

    """

    power_in_W = 10.0 ** (power_in_dBm / 10.0) / 1000.0

    return power_in_W


def w2dBm(power_in_W):
    """Convert power in W to power in dbM

    Convert a microwave power given in W to dBm

    Args:
        power_in_W:   Power in (W)

    """

    power_in_dBm = 10.0 * _np.log10(1000 * power_in_W)

    return power_in_dBm


def B2Hz(B, g=g_e):
    """Convert magnetic field value to Hz

    Convert a magnetic field value in (T) to a frequency given in (Hz)

    Args:
        B:   Magnetic field value in (T)
        g:  g value. If no value is given g_free is used.

    Returns:
        single value or ndarray:     Magnetic field value in (Hz)

    """

    v = g * mu_B * B / (h)

    return v


def Hz2B(Hz, g=g_e):
    """Convert frequency to magnetic field value

    Convert a magnetic field value in (T) to a frequency given in (Hz)

    Args:
        B:   Magnetic field value in (T)
        g:  g value. If no value is given g_free is used.

    Returns:
        single value or ndarray:     Magnetic field value in (T)

    """

    B = h * Hz / (g * mu_B)

    return B
