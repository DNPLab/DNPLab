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


def exponential_window(all_data, dim, lw):
    """Calculate exponential window function

    .. math::
        \mathrm{exponential}  &=  \exp(-2t * \mathrm{linewidth}) &

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
        \mathrm{gaussian}  &=  \exp((\mathrm{linewidth[0]} * t) - (\mathrm{linewidth[1]} * t^{2})) &

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
        \mathrm{hamming}  &=  0.53836 + 0.46164\cos(\pi * n/(N-1)) &

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
        \mathrm{han}  &=  0.5 + 0.5\cos(\pi * n/(N-1)) &

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
        \mathrm{sin2}  &=  \cos((-0.5\pi * n/(N - 1)) + \pi)^{2} &

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
        \mathrm{traf}           &=  (f1 * (f1 + f2)) / (f1^{2} + f2^{2}) &

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


def t1_function(t_axis, T1, M_0, M_inf):
    return M_0 - M_inf * _np.exp(-1.0 * t_axis / T1)


def t2_function(x_axis, M_0, T2, p):
    return M_0 * _np.exp(-2.0 * (x_axis / T2) ** p)


def monoexp_fit(x_axis, C1, C2, tau):
    return C1 + C2 * _np.exp(-1.0 * x_axis / tau)


def biexp_fit(x_axis, C1, C2, tau1, C3, tau2):
    return C1 + C2 * _np.exp(-1.0 * x_axis / tau1) + C3 * _np.exp(-1.0 * x_axis / tau2)


def buildup_function(power_array, E_max, power_half):
    return E_max * power_array / (power_half + power_array)


def interpolate_T1(
    E_powers=False,
    T1_powers=False,
    T1_array=False,
    interp_method="linear",
    spin_C=100,
    T10=2.0,
    T100=2.5,
):
    """Returns interpolated T1 data using Eq. 39 of http://dx.doi.org/10.1016/j.pnmrs.2013.06.001 for "linear" or Eq. 22 of https://doi.org/10.1016/bs.mie.2018.09.024 for "second_order"

    Args:
        E_powers: The x-coordinates at which to evaluate.
        T1_powers: The x-coordinates of the data points, must be increasing.
            Otherwise, T1_power is internally sorted.
        T1_array: The y-coordinates of the data points, same length as T1_power.
        interp_method: "second_order" or "linear".
        spin_C: unpaired electron spin concentration in uM.
        T10: T1 measured with unpaired electrons.
        T100: T1 measured without unpaired electrons.

    Returns:
        interp_T1 (np.array): The evaluated values, same shape as E_powers.
    """

    spin_C = spin_C / 1e6

    # 2nd order fit, Franck and Han MIE (Eq. 22) and (Eq. 23)
    if interp_method == "second_order":

        delta_T1_water = T1_array[-1] - T1_array[0]
        T1_water = T100
        macro_C = spin_C

        kHH = (1.0 / T10 - 1.0 / T1_water) / macro_C
        krp = (
            (1.0 / T1_array)
            - (1.0 / (T1_water + delta_T1_water * T1_powers))
            - (kHH * (macro_C))
        ) / (spin_C)

        p = _np.polyfit(T1_powers, krp, 2)
        T1_fit_2order = _np.polyval(p, E_powers)

        interp_T1 = 1.0 / (
            ((spin_C) * T1_fit_2order)
            + (1.0 / (T1_water + delta_T1_water * E_powers))
            + (kHH * (macro_C))
        )

    # linear fit, Franck et al. PNMRS (Eq. 39)
    elif interp_method == "linear":

        linear_t1 = 1.0 / ((1.0 / T1_array) - (1.0 / T10) + (1.0 / T100))

        p = _np.polyfit(T1_powers, linear_t1, 1)
        T1_fit_linear = _np.polyval(p, E_powers)

        interp_T1 = T1_fit_linear / (
            1.0 + (T1_fit_linear / T10) - (T1_fit_linear / T100)
        )

    else:
        raise Exception("invalid interp_method")

    return interp_T1
