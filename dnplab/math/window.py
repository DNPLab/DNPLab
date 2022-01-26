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
    data, _ = return_data(all_data)
    T2 = 1.0 / (np.pi * lw)
    t = data.coords[dim]
    T = np.max(t)
    E = np.exp(-1 * t / T2)
    e = np.exp(-1 * (T - t) / T2)
    return E * (E + e) / (E ** 2 + e ** 2)



def gaussian_window(all_data, dim, lw):
    """Calculate gaussian window function

    .. math::
        \mathrm{gaussian} = e^{((\mathrm{linewidth}[0] * t) - (\mathrm{linewidth}[1] * t^{2}))}

    Args:
        all_data (dnpdata, dict): data container
        dim (str): dimension to window

    Returns:
        array: gaussian window function
    """
    if (
        not isinstance(lw, list)
        or len(lw) != 2
        or any([isinstance(x, list) for x in lw])
    ):
        raise ValueError("lw must a list with len=2 for the gaussian window")
    else:
        data, _ = return_data(all_data)
        return _np.exp((lw[0] * data.coords[dim]) - (lw[1] * data.coords[dim] ** 2))


def hamming_window(dim_size):
    """Calculate hamming window function

    .. math::
        \mathrm{hamming} = 0.53836 + 0.46164\cos(\pi * n / (N-1))

    Args:
        dim_size(int): length of array to window

    Returns:
        array: hamming window function
    """
    return 0.53836 + 0.46164 * _np.cos(
        1.0 * _np.pi * _np.arange(dim_size) / (dim_size - 1)
    )


def hann_window(dim_size):
    """Calculate hann window function

    .. math::
        \mathrm{han} = 0.5 + 0.5\cos(\pi * n / (N-1))

    Args:
        dim_size(int): length of array to window

    Returns:
        array: hann window function
    """
    return 0.5 + 0.5 * _np.cos(1.0 * _np.pi * _np.arange(dim_size) / (dim_size - 1))


def lorentz_gauss_window(all_data, dim, exp_lw, gauss_lw, gaussian_max=0):
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
    data, _ = return_data(all_data)
    dim_size = data.coords[dim].size
    expo = _np.pi * data.coords[dim] * exp_lw
    gaus = 0.6 * _np.pi * gauss_lw * (gaussian_max * (dim_size - 1) - data.coords[dim])
    return _np.exp(expo - gaus ** 2).reshape(dim_size)


def sin2_window(dim_size):
    """Calculate sin-squared window function

    .. math::
        \sin^{2}  =  \cos((-0.5\pi * n / (N - 1)) + \pi)^{2}

    Args:
        dim_size(int): length of array to window

    Returns:
        array: sin-squared window function
    """
    return (
        _np.cos((-0.5 * _np.pi * _np.arange(dim_size) / (dim_size - 1)) + _np.pi) ** 2
    )


def traf_window(all_data, dim, traf_lw):
    """Calculate traf window function

    .. math::
        \mathrm{traf}  &=  (f1 * (f1 + f2)) / (f1^{2} + f2^{2}) &

               f1(t)   &=  \exp(-t * \pi * \mathrm{linewidth[0]}) &

               f2(t)   &=  \exp((t - T) * \pi * \mathrm{linewidth[1]}) &


    Args:
        all_data (dnpdata, dict): data container
        dim (str): dimension to window
        exp_lw (int or float): exponential linewidth
        gauss_lw (int or float): gaussian linewidth

    Returns:
        array: traf window function
    """
    data, _ = return_data(all_data)
    T2 = 1.0 / (_np.pi * traf_lw)
    t = data.coords[dim]
    T = _np.max(t)
    E = _np.exp(-1 * t / T2)
    e = _np.exp(-1 * (T - t) / T2)
    return E * (E + e) / (E ** 2 + e ** 2)