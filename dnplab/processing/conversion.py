import numpy as _np


def dBm2w(power_in_dBm):
    """Convert power in dBm to power in W

    Convert a microwave power given in dBm to W. Note that for values <= -199 dBm the output power in W is set to 0.

    Args:
        power_in_dBm: Power in (dBm)

    Returns:
        float (array): Power in (W)

    """
    print('This is dBm2w')

    power_in_W = 10.0 ** (power_in_dBm / 10.0) / 1000.0

    # Set values below -199 dBm to 0 W
    index = _np.where(power_in_dBm <= -199)
    power_in_W[index] = 0.0

    return power_in_W


def w2dBm(power_in_W):
    """Convert power in W to power in dbM

    Convert a microwave power given in W to dBm

    Args:
        power_in_W:   Power in (W)

    """

    power_in_dBm = 10.0 * _np.log10(1000 * power_in_W)

    return power_in_dBm
