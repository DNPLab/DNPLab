import numpy as _np

def dBm2w(power_in_dBm):
    """Convert power in dBm to power in W

    Convert a microwave power given in dBm to W

    Args:
        power_in_dBm:   Power in (dBm)

    """

    power_in_W =  10.0**(power_in_dBm/10.0)/1000.0

    return power_in_W

def w2dBm(power_in_W):
    """Convert power in W to power in dbM

    Convert a microwave power given in W to dBm

    Args:
        power_in_W:   Power in (W)

    """

    power_in_dBm =  10.0*_np.log10(1000 * power_in_W)

    return power_in_dBm



