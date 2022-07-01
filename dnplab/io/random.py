import numpy as np
from ..core.data import DNPData


def fid(points=1024, snr=100.0):
    """Generate an FID dataset for testing

    Args:
        points (int): Number of points in FID
        snr (float): Signal to noise ratio of FID

    Returns:
        DNPData: FID in DNPData object
    """

    t2 = np.r_[0 : 1 : 1j * points]

    values = (
        np.exp(1j * 2 * np.pi * 100.0 * t2) * np.exp(-1 * t2 / 0.10)
        + np.random.randn(points) / snr
    )

    attrs = {"nmr_frequency": 300e6}

    return DNPData(values, ["t2"], [t2])


def ir(points=(1024,), snr=10):
    """Generate an IR dataset for testing"""
    return NotImplemented


def nd(*args):
    """Generate n-dimensional dataset for testing"""
    return NotImplemented
