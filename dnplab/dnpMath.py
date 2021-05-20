import numpy as _np

from . import dnpdata, dnpdata_collection


def return_data(all_data):

    is_workspace = False
    if isinstance(all_data, dnpdata):
        data = all_data.copy()
    elif isinstance(all_data, dict):
        raise ValueError("Type dict is not supported")
    elif isinstance(all_data, dnpdata_collection):
        is_workspace = True
        if all_data.processing_buffer in all_data.keys():
            data = all_data[all_data.processing_buffer]
        else:
            raise ValueError("No data in processing buffer")
    else:
        raise ValueError("Data type not supported")

    return data, is_workspace


def convert_power(watts=False, dBm=False, loss=0):
    """Convert between Watts and dBm

    .. math::
        \mathrm{dBm} =  10 * log((\mathrm{watts} + \mathrm{loss}) * 1000)
        \mathrm{Watts} =  1E^{-3} * 10^{(\mathrm{dBm} + \mathrm{loss}) / 10}

    Args:
        watts (float, int, list, array): microwave power(s) in Watts
        dBm (float, int, list, array): microwave power(s) in dBm
        loss (float, int, list, array): constant loss characteristic of your device, in same units given to convert from

    Returns:
        array: converted microwave power(s)
    """

    if watts and dBm:
        raise ("Give powers in watts OR in dBm, not both together")

    if dBm:
        powers = _np.add(dBm, loss)
        powers = _np.divide(powers, 10)
        powers = _np.power(10, powers)
        powers = _np.multiply(1e-3, powers)

    if watts:
        powers = 10 * _np.log((watts + loss) * 1e3)

    return powers


def exponential_window(all_data, dim, lw):
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
    data, _ = return_data(all_data)
    return _np.exp(-2 * data.coords[dim] * lw)


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


def traf_window(all_data, dim, exp_lw, gauss_lw):
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
    E_t = _np.exp(-1 * data.coords[dim] * _np.pi * exp_lw)
    e_t = _np.exp((data.coords[dim] - max(data.coords[dim])) * _np.pi * gauss_lw)
    return (E_t * (E_t + e_t)) / ((E_t ** 2) + (e_t ** 2))


def t1_function(t, T1, M_0, M_inf):
    """Calculate exponential T1 curve

    .. math::
        f(t) = M_0 - M_{\infty} e^{-t/T_{1}}

    Args:
        t (array): time series
        T_{1} (float): T1 value
        M_{0} (float): see equation
        M_{\infty} (float): see equation

    Returns:
        array: T1 curve
    """

    return M_0 - M_inf * _np.exp(-1.0 * t / T1)


def t2_function(t, M_0, T2, p):
    """Calculate stretched or un-stretched (p=1) exponential T2 curve

    .. math::
        f(t) = M_{0} e^{(-2(t/T_{2})^{p}}

    Args:
        t (array): time series
        M_{0} (float): see equation
        T_{2} (float): T2 value
        p (float): see equation

    Returns:
        array: T2 curve
    """

    return M_0 * _np.exp(-2.0 * (t / T2) ** p)


def monoexp_fit(t, C1, C2, tau):
    """Calculate mono-exponential curve

    .. math::
        f(t) = C1 + C2 e^{-t/tau}

    Args:
        t (array): time series
        C1 (float): see equation
        C2 (float): see equation
        tau (float): see equation

    Returns:
        array: mono-exponential curve
    """

    return C1 + C2 * _np.exp(-1.0 * t / tau)


def biexp_fit(t, C1, C2, tau1, C3, tau2):
    """Calculate bi-exponential curve

    .. math::
        f(t) = C1 + C2 e^{-t/tau1} + C3 e^{-t/tau2}

    Args:
        t (array): time series
        C1 (float): see equation
        C2 (float): see equation
        C3 (float): see equation
        tau1 (float): see equation
        tau2 (float): see equation

    Returns:
        array: bi-exponential curve
    """

    return C1 + C2 * _np.exp(-1.0 * t / tau1) + C3 * _np.exp(-1.0 * t / tau2)


def buildup_function(p, E_max, p_half):
    """Calculate asymptotic buildup curve

    .. math::
        f(p) = E_{max} * p / (p_{1/2} + p)

    Args:
        p (array): power series
        E_{max} (float): maximum enhancement
        p_{1/2} (float): power at half saturation

    Returns:
        array: buildup curve
    """

    return E_max * p / (p_half + p)
