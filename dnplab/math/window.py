import numpy as np

__all__ = ['exponential']

def exponential(x, lw):
    """Calculate exponential window function

    .. math::
        \mathrm{exponential} =  e^{-2t * \mathrm{linewidth}}

    Args:
        all_data (dnpdata, dict): data container
        dim (str): dimension to window
        lw (int or float): linewidth

    Returns:
        array: exponential window function
    """
    return np.exp(-1. * x * lw)

def gaussian(x, lw):
    """Calculate gaussian window function

    .. math::
        \mathrm{gaussian} = e^{(\sigma * x^{2})}

    Args:
        x (array_like): Vector of points
        lw (float): Standard deviation of gaussian window

    Returns:
        array: gaussian window function
    """
    return np.exp((lw * x) ** 2)

def hann(x):
    """Calculate hann window function

    .. math::
        \mathrm{han} = 0.5 + 0.5\cos(\pi * n / (N-1))

    Args:
        N(int): number of points to return in window function

    Returns:
        array: hann window function
    """

    if type(x) == int:
        N = x
    else:
        N = len(x)

    return 0.5 + 0.5 * np.cos(1.0 * np.pi * np.arange(N) / (N - 1))


def traf(x, lw):
    """Calculate traf window function

    .. math::
        \mathrm{traf}  &=  (f1 * (f1 + f2)) / (f1^{2} + f2^{2}) &

               f1(t)   &=  \exp(-t * \pi * \mathrm{linewidth[0]}) &

               f2(t)   &=  \exp((t - T) * \pi * \mathrm{linewidth[1]}) &


    Args:
        x (array_like): axis for traficant window
        lw (str): linewidth of traficant window

    Returns:
        np.ndarray: traf window function
    """
    data, _ = return_data(all_data)
    T2 = 1.0 / (np.pi * lw)
    t = data.coords[dim]
    T = np.max(t)
    E = np.exp(-1 * t / T2)
    e = np.exp(-1 * (T - t) / T2)
    return E * (E + e) / (E ** 2 + e ** 2)



