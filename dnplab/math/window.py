import numpy as np

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
    return np.exp(-1. * (x - x[0]) * lw)

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
    T2 = 1.0 / (np.pi * lw)
    t = x
    T = np.max(t)
    E = np.exp(-1 * t / T2)
    e = np.exp(-1 * (T - t) / T2)
    return E * (E + e) / (E ** 2 + e ** 2)


def hamming(x):
    """Calculate hamming window function

    .. math::
        \mathrm{hamming} = 0.53836 + 0.46164\cos(\pi * n / (N-1))

    Args:
        dim_size(int): length of array to window

    Returns:
        array: hamming window function
    """
    if type(x) == int:
        N = x
    else:
        N = len(x)

    return 0.53836 + 0.46164 * np.cos(
        1.0 * np.pi * np.arange(N) / (N - 1)
    )

def lorentz_gauss(x, exp_lw, gauss_lw, gaussian_max=0):
    """Calculate lorentz-gauss window function

    .. math::
        \mathrm{lorentz\_gauss} &=  \exp(L -  G^{2}) &

           L(t)    &=  \pi * \mathrm{linewidth[0]} * t &

           G(t)    &=  0.6\pi * \mathrm{linewidth[1]} * (\mathrm{gaussian\_max} * (N - 1) - t) &


    Args:
        all_data (dnpdata, dict): data container
        dim (str): dimension to window
        exp_lw (int or float): exponential linewidth
        gauss_lw (int or float): gaussian linewidth
        gaussian_max (int): location of maximum in gaussian window

    Returns:
        array: gauss_lorentz window function
    """

    expo = np.pi * x * exp_lw
    gaus = 0.6 * np.pi * gauss_lw * (gaussian_max * (N - 1) - x)
    return np.exp(expo - gaus ** 2).reshape(N)


def sin2(x):
    """Calculate sin-squared window function

    .. math::
        \sin^{2}  =  \cos((-0.5\pi * n / (N - 1)) + \pi)^{2}

    Args:
        dim_size(int): length of array to window

    Returns:
        array: sin-squared window function
    """
    if type(x) == int:
        N = x
    else:
        N = len(x)

    return (
        np.cos((-0.5 * np.pi * np.arange(N) / (N - 1)) + np.pi) ** 2
    )