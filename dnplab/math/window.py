from . import dnpdata

import numpy as np

__all__ = ["exponential"]


def exponential(x, lw, centered=False):
    """Exponential Window Function

    Args:
        x (np.ndarray): Coordinates for exponential window function
        gamma (np.ndarray):
    """

    tau = lw * np.pi

    return np.exp(-1 * x / tau)


def rectangular(x, length, centered=False):
    """Rectangular Window Function"""

    if centered:
        center = (x[-1] - x[0]) / 2
        rect = [abs(center - z) < (0.5 * length) for z in x]

    else:
        rect = [(z < length) for z in x]

    return np.array(rect)


def triangular(x, centered=False):
    """Triangular Window Function"""
    N = len(x)
    n = np.array(range(N))
    if centered:
        tri = 1 - abs((n - N / 2) / (N / 2))
    else:
        tri = -1 * (n - N) / (N)

    return np.array(tri)


def sine(x, centered=False):
    """Sine Window Function"""
    N = len(x)
    n = np.array(range(N))

    if centered:
        s = np.sin(np.pi * n / N)
    else:
        s = np.cos(np.pi * n / N / 2)

    return s
